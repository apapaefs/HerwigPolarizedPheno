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

  /// Exact 900-cell HERMES unpolarized cos(phi) and cos(2phi) diagnostics.
  class HERMES_2013_I1111237 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERMES_2013_I1111237);

    struct Histograms {
      Histo1DPtr numerator, denominator, covariancePositive, covarianceNegative;
    };

    struct EventValues {
      std::vector<double> numerator = std::vector<double>(900, 0.0);
      std::vector<double> denominator = std::vector<double>(900, 0.0);
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");
      _xEdges = {.023,.042,.078,.145,.27,.60};
      _yEdges = {.20,.30,.45,.60,.70,.85};
      _zEdges = {.20,.30,.40,.50,.60,.75,1.00};
      _ptEdges = {.05,.20,.35,.50,.70,1.00,1.30};
      for (const std::string species :
           {"hplus","hminus","piplus","piminus","kplus","kminus"})
        for (int harmonic : {1, 2})
          bookDataset(species, harmonic);
      book(_acceptedX, "Accepted_X", _xEdges);
      book(_acceptedY, "Accepted_Y", _yEdges);
      book(_acceptedQ2, "Accepted_Q2", logspace(40, 1.0, 20.0));
      book(_acceptedW2, "Accepted_W2", 40, 10.0, 50.0);
      book(_acceptedZ, "Accepted_Z", _zEdges);
      book(_acceptedPt, "Accepted_HadronPt", _ptEdges);
      book(_acceptedXF, "Accepted_XF", 32, .2, 1.0);
      book(_acceptedPhi, "Accepted_HadronPhi", 32, -M_PI, M_PI);
      _scaled.insert(_scaled.end(), {
        _acceptedX, _acceptedY, _acceptedQ2, _acceptedW2,
        _acceptedZ, _acceptedPt, _acceptedXF, _acceptedPhi});
    }

    void analyze(const Event& event) {
      const Particles prompt =
        apply<PromptFinalState>(event, "PromptPositrons").particles();
      const COMPASSSIDIS::DISKinematics dis =
        COMPASSInclusiveDIS::reconstructForLeptons(event, prompt, {-11});
      if (!dis.valid || !COMPASSSIDIS::fixedBeamEnergy(dis, 27.6)) vetoEvent;
      if (!(dis.Q2 > 1.0 && dis.W2 > 10.0 &&
            dis.x >= .023 && dis.x < .60 &&
            dis.y >= .20 && dis.y < .85)) vetoEvent;
      const int ix = SIDISAzimuthal::binIndex(dis.x, _xEdges);
      const int iy = SIDISAzimuthal::binIndex(dis.y, _yEdges);
      if (ix < 0 || iy < 0) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedY->fill(dis.y);
      _acceptedQ2->fill(dis.Q2);
      _acceptedW2->fill(dis.W2);
      std::map<std::string, EventValues> values;
      for (const auto& entry : _histograms) values[entry.first] = EventValues();

      for (const Particle& particle :
           apply<FinalState>(event, "LabFS").particles()) {
        if (!PID::isHadron(particle.pid()) || particle.charge3() == 0) continue;
        const int absPid = std::abs(particle.pid());
        const bool pion = absPid == 211;
        const bool kaon = absPid == 321;
        const auto hadron = COMPASSSIDIS::hadronKinematics(
          dis, particle, !pion && !kaon);
        if (!hadron.valid) continue;
        const double xF = SIDISAzimuthal::feynmanX(dis, particle);
        if (!(xF > .2)) continue;
        const int iz = SIDISAzimuthal::binIndex(hadron.z, _zEdges);
        const int ipt = SIDISAzimuthal::binIndex(
          hadron.transverseMomentum, _ptEdges);
        if (iz < 0 || ipt < 0) continue;
        const double momentum = hadron.momentum;
        const bool unidentifiedAcceptance = momentum > 1.0 && momentum < 15.0;
        const bool pionAcceptance = pion && momentum > 1.0 && momentum < 15.0;
        const bool kaonAcceptance = kaon && momentum > 2.0 && momentum < 15.0;
        if (!unidentifiedAcceptance && !pionAcceptance && !kaonAcceptance)
          continue;
        const SIDISAzimuthal::Angle phi =
          SIDISAzimuthal::hadronAzimuth(dis, particle);
        if (!phi.valid) continue;
        const int flat = ((ix*5 + iy)*6 + iz)*6 + ipt;
        std::vector<std::string> species;
        if (unidentifiedAcceptance)
          species.push_back(particle.charge3() > 0 ? "hplus" : "hminus");
        if (pionAcceptance)
          species.push_back(particle.pid() > 0 ? "piplus" : "piminus");
        if (kaonAcceptance)
          species.push_back(particle.pid() > 0 ? "kplus" : "kminus");
        for (const std::string& name : species) {
          for (int harmonic : {1, 2}) {
            EventValues& target = values[key(name, harmonic)];
            target.numerator[size_t(flat)] +=
              SIDISAzimuthal::cosineHarmonic(phi, harmonic);
            target.denominator[size_t(flat)] += 1.0;
          }
        }
        _acceptedZ->fill(hadron.z);
        _acceptedPt->fill(hadron.transverseMomentum);
        _acceptedXF->fill(xF);
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
    static std::string key(const std::string& species, int harmonic) {
      return species + "_cos" + std::to_string(harmonic);
    }

    void bookDataset(const std::string& species, int harmonic) {
      const std::string identifier = key(species, harmonic);
      Histograms& histograms = _histograms[identifier];
      book(histograms.numerator, "MomentNumerator_" + identifier,
           900, 0.0, 900.0);
      book(histograms.denominator, "MomentDenominator_" + identifier,
           900, 0.0, 900.0);
      book(histograms.covariancePositive, "CovariancePositive_" + identifier,
           900, 0.0, 900.0);
      book(histograms.covarianceNegative, "CovarianceNegative_" + identifier,
           900, 0.0, 900.0);
      _scaled.insert(_scaled.end(), {
        histograms.numerator, histograms.denominator,
        histograms.covariancePositive, histograms.covarianceNegative});
    }

    std::vector<double> _xEdges, _yEdges, _zEdges, _ptEdges;
    std::map<std::string, Histograms> _histograms;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _acceptedX, _acceptedY, _acceptedQ2, _acceptedW2;
    Histo1DPtr _acceptedZ, _acceptedPt, _acceptedXF, _acceptedPhi;
  };

  RIVET_DECLARE_PLUGIN(HERMES_2013_I1111237);
}
