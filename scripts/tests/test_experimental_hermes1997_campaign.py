from __future__ import annotations

import contextlib
import copy
import io
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
import run_experimental_campaign as campaign
import phenomenology_reference_data as reference
import run_HERMES_1997_I440904 as wrapper
import check_hermes1997_window_sensitivity as sensitivity

ID = "HERMES_1997_I440904"


class HERMES1997Tests(unittest.TestCase):
    def setUp(self):
        self.measurement = campaign.get_measurement(ID)
        self.snapshot = reference.validate_vendored(ID)
        self.dataset = self.snapshot["datasets"][0]

    def test_table_I_A1_column_and_uncertainties_match_publication(self):
        points = self.dataset["points"]
        self.assertEqual([p["value"] for p in points],
                         [-.111,-.117,-.077,-.126,-.097,-.158,-.078,.146,-.149])
        self.assertEqual([p["stat"] for p in points],
                         [.048,.052,.055,.064,.068,.085,.113,.219,.374])
        self.assertEqual([p["systematic_combined"] for p in points],
                         [.018,.013,.011,.014,.015,.020,.019,.052,.103])
        self.assertEqual([p["x_mean"] for p in points],
                         [.033,.047,.065,.087,.119,.168,.244,.342,.464])
        self.assertEqual([p["q2_mean"] for p in points],
                         [1.22,1.47,1.73,1.99,2.30,2.65,3.07,3.86,5.25])
        self.assertEqual(self.dataset["rivet_path"], f"/{ID}/d01-x01-y01")

    def test_production_budget_and_neutron_only_cards(self):
        jobs = campaign.build_job_matrix(self.measurement,3000000,300000,100,440904)
        self.assertEqual(len(jobs),800)
        self.assertEqual(sum(j["events"] for j in jobs),13200000)
        self.assertEqual(len({j["seed"] for j in jobs}),800)
        self.assertEqual({j["component"] for j in jobs},{"N"})
        self.assertEqual(len({j["stem"] for j in jobs}),8)
        self.assertEqual(set(self.measurement["outputs"]),
                         {"Neutron_W025","Neutron_W050","Neutron_W100"})
        jobs = campaign.build_job_matrix(self.measurement,1000,100,1,10,True)
        self.assertEqual(len(jobs),30)
        for job in jobs:
            card = (ROOT/job["card_source"]).read_text()
            self.assertIn("saverun "+job["stem"]+" EventGenerator",card)
            self.assertIn("TargetParticle /Herwig/Particles/n0",card)
            if job["family"] == "nominal":
                contribution = "PositiveNLO" if job["order"] == "POSNLO" else "NegativeNLO"
                self.assertIn("PowhegMEDISNCPol:Contribution "+contribution,card)

    def test_card_physics_and_generation_envelope(self):
        common=(ROOT/self.measurement["cards"]["directory"]/self.measurement["cards"]["common"]).read_text()
        for setting in ["BeamA /Herwig/Particles/e+", "BeamB /Herwig/Particles/n0",
                        "BeamEMaxA 27.5*GeV", "MaxQ2 15.0*GeV2",
                        "SubProcess:MatrixElements[0] /Herwig/MatrixElements/PowhegMEDISNCPol",
                        "PowhegMEDISNCPol:UsePOWHEGRealSpinVertex Yes",
                        "ShowerHandler:SpinCorrelations Yes"]:
            self.assertIn(setting,common)
        self.assertNotIn("25.51",common)
        self.assertLess(max(self.measurement["physics"]["q2_means_gev2"])*1.1,15.)

    def test_no_second_helium_correction(self):
        family=campaign.campaign_family_specs(self.measurement,False)["nominal"]
        for selection in self.measurement["outputs"]:
            for channel in ("unpolarized","longitudinal"):
                self.assertEqual(campaign._target_component_coefficients(
                    self.measurement,family,channel,selection),{"N":1.0})
        self.assertNotIn("effective_polarizations",self.measurement["physics"])

    def test_yoda_reference_aliases_and_extraction_annotation(self):
        try:
            import yoda
        except ImportError:
            self.skipTest("YODA is not active")
        objects=yoda.read(str(ROOT/self.measurement["analysis"]["reference_yoda"]))
        self.assertEqual(len(objects),3)
        for path in [self.dataset["rivet_path"],*self.dataset["rivet_aliases"]]:
            obj=objects["/REF"+path]
            self.assertEqual(list(obj.xEdges()),self.dataset["bin_edges"])
            self.assertIn("central values neglect A2",obj.annotation("A2G2Assumption"))
            for b,p in zip(obj.bins(),self.dataset["points"]):
                self.assertAlmostEqual(b.val(),p["value"])
                self.assertAlmostEqual(b.errAvg("stat"),p["stat"])
                self.assertAlmostEqual(b.errAvg("published_total"),p["systematic_combined"])
            overlay=campaign._fixed_reference_overlay_points(self.measurement,self.snapshot,Path(path).name)
            self.assertEqual([p["plot_x"] for p in overlay],[p["x_mean"] for p in self.dataset["points"]])

    def test_reference_rejects_wrong_schema_target_energy_or_kinematics(self):
        table=json.loads((ROOT/reference.SOURCES[ID][0]["raw"]).read_text())
        mutations=[
            lambda t:t["headers"][2].update(name="G1"),
            lambda t:t["qualifiers"]["RE"][0].update(value="E+ HE3 --> E+ X"),
            lambda t:t["qualifiers"]["PLAB"][0].update(value="25.51 GeV"),
            lambda t:t["values"][0]["x"][0].update(value="0.018"),
            lambda t:t["values"][0]["x"][1].update(value="5.0"),
        ]
        for mutate in mutations:
            bad=copy.deepcopy(table); mutate(bad)
            with self.assertRaises(reference.ReferenceDataError):
                reference.normalize_hermes1997(bad)

    def test_wrapper_dry_run_no_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            with mock.patch.object(campaign,"CAMPAIGN_ROOT",Path(directory)),contextlib.redirect_stdout(io.StringIO()) as stream:
                status=wrapper.main(["--output-root",directory,"full","--tag","check","--dry-run"])
            self.assertEqual(status,0)
            plan=json.loads(stream.getvalue())
            self.assertEqual(plan["campaign_directory"],str(Path(directory)/ID/"check"))
            self.assertEqual(list(Path(directory).iterdir()),[])

    def test_postprocessing_and_window_report_roundtrip(self):
        try:
            import yoda
        except ImportError:
            self.skipTest("YODA is not active")
        family=campaign.campaign_family_specs(self.measurement,False)["nominal"]
        edges=self.dataset["bin_edges"]
        def load(groups,directory,analysis,name,family_id,component,helicities,orders,*args):
            self.assertEqual(component,"N")
            factor,varfactor=1.,1.
            if name.startswith("SigmaOverD"):
                factor,varfactor=2.,4.
            elif name.startswith("CovarianceProxy"):
                factor,varfactor=2**.5,2.
            return {h:campaign.BinSeries(edges,[(20+(-3 if h in ("PP","MM") else 3))*factor]*9,
                                         [.01*varfactor]*9) for h in helicities}
        with mock.patch.object(campaign,"_load_component_series",side_effect=load):
            results,diagnostics=campaign._load_asymmetry_family_products({},Path('.'),self.measurement,"nominal",family)
        objects=campaign._asymmetry_family_objects(yoda,self.measurement,self.snapshot,results,diagnostics,
                                                  campaign._family_annotations(self.measurement,"nominal",family))
        with tempfile.TemporaryDirectory() as temporary:
            dest=Path(temporary)/'fixture.yoda'
            campaign._write_yoda_objects(yoda,objects,dest)
            loaded=yoda.read(str(dest))
            for selection,output in self.measurement['outputs'].items():
                obj=loaded[output['a1']]
                self.assertAlmostEqual(obj.bin(1).val(),-.3)
                self.assertEqual(obj.annotation('Observable'),'A1n')
                self.assertIn('matching published central-value convention',obj.annotation('MCExtractionAssumption'))
            csv=Path(temporary)/'summary.csv'
            campaign._write_summary_csv(csv,[row for selection,result in results.items()
                                            for row in campaign._summary_rows(selection,result)])
            report=sensitivity.compare(csv)
            self.assertEqual(len(report['records']),18)
            self.assertTrue(all(r['delta_a1']==0 for r in report['records']))


if __name__ == '__main__':
    unittest.main()
