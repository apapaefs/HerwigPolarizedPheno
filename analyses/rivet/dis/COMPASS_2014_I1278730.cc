// -*- C++ -*-
#include "COMPASSSIDIS.hh"
#include "SIDISAzimuthal.hh"
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  /// Unpolarized charged-hadron cos(phi) and cos(2phi) amplitudes.
  class COMPASS_2014_I1278730 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(COMPASS_2014_I1278730);

    struct Histograms {
      Histo1DPtr numerator, denominator, covariancePositive, covarianceNegative;
    };

    struct EventValues {
      std::vector<double> numerator, denominator;
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -13), "PromptMuons");
      _xEdges = {.003,.008,.013,.020,.032,.050,.080,.130};
      _zEdges = {.20,.25,.30,.34,.38,.42,.49,.63,.85};
      _ptEdges = {.10,.20,.27,.33,.39,.46,.55,.64,.77,1.00};
      _x3Edges = {.003,.012,.020,.038,.130};
      _z3Edges = {.20,.25,.32,.40,.55,.70,.85};
      _pt3Edges = {.10,.30,.50,.64,1.00};
      for (const std::string charge : {"hplus", "hminus"}) {
        for (int harmonic : {1, 2}) {
          bookDataset(charge, harmonic, "x", _xEdges.size()-1);
          bookDataset(charge, harmonic, "z", _zEdges.size()-1);
          bookDataset(charge, harmonic, "pt", _ptEdges.size()-1);
          bookDataset(charge, harmonic, "x_z_pt",
                      (_x3Edges.size()-1)*(_z3Edges.size()-1)
                      *(_pt3Edges.size()-1));
        }
      }
      book(_acceptedX, "Accepted_X", _xEdges);
      book(_acceptedQ2, "Accepted_Q2", logspace(40, 1.0, 30.0));
      book(_acceptedY, "Accepted_Y", 35, .2, .9);
      book(_acceptedZ, "Accepted_Z", _zEdges);
      book(_acceptedPt, "Accepted_HadronPt", _ptEdges);
      book(_acceptedPhi, "Accepted_HadronPhi", 32, -M_PI, M_PI);
      _scaled.insert(_scaled.end(), {
        _acceptedX, _acceptedQ2, _acceptedY,
        _acceptedZ, _acceptedPt, _acceptedPhi});
    }

    void analyze(const Event& event) {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptMuons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstruct(event, prompt);
      if (!dis.valid || !COMPASSSIDIS::fixed160GeVBeam(dis)) vetoEvent;
      if (!(dis.Q2 > 1.0 && dis.W2 > 25.0 &&
            dis.x > .003 && dis.x < .13 &&
            dis.y > .2 && dis.y < .9)) vetoEvent;
      if (!(dis.k.angle(dis.q) < .060)) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      std::map<std::string, EventValues> values;
      for (const auto& entry : _histograms) {
        const size_t size = entry.second.numerator->numBins();
        values[entry.first] = {
          std::vector<double>(size, 0.0),
          std::vector<double>(size, 0.0),
        };
      }

      for (const Particle& particle :
           apply<FinalState>(event, "LabFS").particles()) {
        if (!PID::isHadron(particle.pid()) || particle.charge3() == 0) continue;
        const auto hadron = COMPASSSIDIS::hadronKinematics(dis, particle, true);
        if (!hadron.valid || hadron.z < .2 || hadron.z >= .85 ||
            hadron.transverseMomentum < .1 ||
            hadron.transverseMomentum >= 1.0) continue;
        const SIDISAzimuthal::Angle phi =
          SIDISAzimuthal::hadronAzimuth(dis, particle);
        if (!phi.valid) continue;
        const std::string charge = particle.charge3() > 0 ? "hplus" : "hminus";
        const int ix = SIDISAzimuthal::binIndex(dis.x, _xEdges);
        const int iz = SIDISAzimuthal::binIndex(hadron.z, _zEdges);
        const int ipt = SIDISAzimuthal::binIndex(
          hadron.transverseMomentum, _ptEdges);
        const int ix3 = SIDISAzimuthal::binIndex(dis.x, _x3Edges);
        const int iz3 = SIDISAzimuthal::binIndex(hadron.z, _z3Edges);
        const int ipt3 = SIDISAzimuthal::binIndex(
          hadron.transverseMomentum, _pt3Edges);
        for (int harmonic : {1, 2}) {
          const double numerator =
            2.0*SIDISAzimuthal::cosineHarmonic(phi, harmonic);
          const double denominator = SIDISAzimuthal::epsilon(dis.y, harmonic);
          accumulate(values[key(charge, harmonic, "x")], ix,
                     numerator, denominator);
          accumulate(values[key(charge, harmonic, "z")], iz,
                     numerator, denominator);
          accumulate(values[key(charge, harmonic, "pt")], ipt,
                     numerator, denominator);
          const int flat =
            (ix3 < 0 || iz3 < 0 || ipt3 < 0) ? -1
            : (ix3*int(_z3Edges.size()-1) + iz3)
              *int(_pt3Edges.size()-1) + ipt3;
          accumulate(values[key(charge, harmonic, "x_z_pt")], flat,
                     numerator, denominator);
        }
        _acceptedZ->fill(hadron.z);
        _acceptedPt->fill(hadron.transverseMomentum);
        _acceptedPhi->fill(std::atan2(phi.sine, phi.cosine));
      }

      for (const auto& entry : values) {
        Histograms& histograms = _histograms.at(entry.first);
        SIDISAzimuthal::fillEventBins(
          histograms.numerator, entry.second.numerator);
        SIDISAzimuthal::fillEventBins(
          histograms.denominator, entry.second.denominator);
        SIDISAzimuthal::fillSignedCovariance(
          histograms.covariancePositive, histograms.covarianceNegative,
          entry.second.numerator, entry.second.denominator);
      }
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& histogram : _scaled) scale(histogram, factor);
    }

  private:
    static std::string key(
        const std::string& charge, int harmonic, const std::string& projection) {
      return charge + "_cos" + std::to_string(harmonic) + "_" + projection;
    }

    void bookDataset(
        const std::string& charge, int harmonic,
        const std::string& projection, size_t bins) {
      const std::string identifier = key(charge, harmonic, projection);
      Histograms& histograms = _histograms[identifier];
      book(histograms.numerator, "MomentNumerator_" + identifier,
           bins, 0.0, double(bins));
      book(histograms.denominator, "DepolarizationDenominator_" + identifier,
           bins, 0.0, double(bins));
      book(histograms.covariancePositive, "CovariancePositive_" + identifier,
           bins, 0.0, double(bins));
      book(histograms.covarianceNegative, "CovarianceNegative_" + identifier,
           bins, 0.0, double(bins));
      _scaled.insert(_scaled.end(), {
        histograms.numerator, histograms.denominator,
        histograms.covariancePositive, histograms.covarianceNegative});
    }

    static void accumulate(
        EventValues& values, int index,
        double numerator, double denominator) {
      if (index < 0) return;
      values.numerator[size_t(index)] += numerator;
      values.denominator[size_t(index)] += denominator;
    }

    std::vector<double> _xEdges, _zEdges, _ptEdges;
    std::vector<double> _x3Edges, _z3Edges, _pt3Edges;
    std::map<std::string, Histograms> _histograms;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY;
    Histo1DPtr _acceptedZ, _acceptedPt, _acceptedPhi;
  };

  RIVET_DECLARE_PLUGIN(COMPASS_2014_I1278730);
}
