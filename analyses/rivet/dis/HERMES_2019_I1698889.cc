// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Tools/Beams.hh"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace Rivet {

  /// Identified-hadron longitudinal double-spin asymmetries in HERMES SIDIS.
  class HERMES_2019_I1698889 : public Analysis {
  public:
    RIVET_DEFAULT_ANALYSIS_CTOR(HERMES_2019_I1698889);

    struct DISKinematics {
      bool valid = false;
      int targetPid = 0;
      FourMomentum k, kp, target, q;
      double Q2 = -1.0;
      double x = -1.0;
      double y = -1.0;
      double W2 = -1.0;
      double theta = -1.0;
    };

    struct HadronKinematics {
      bool valid = false;
      double z = -1.0;
      double pt = -1.0;
      double xF = -1.0;
      double cosPhi = 0.0;
      double momentum = -1.0;
    };

    struct SpeciesHistograms {
      Histo1DPtr x, xz, xpt, xzpt;
      Histo1DPtr moment, momentDenominator;
      Histo1DPtr momentCovariancePositive, momentCovarianceNegative;
    };

    struct EventSpecies {
      std::vector<double> x = std::vector<double>(9, 0.0);
      std::vector<double> xz = std::vector<double>(21, 0.0);
      std::vector<double> xpt = std::vector<double>(18, 0.0);
      std::vector<double> xzpt = std::vector<double>(81, 0.0);
      std::vector<double> moment = std::vector<double>(30, 0.0);
      std::vector<double> momentDenominator = std::vector<double>(30, 0.0);
    };

    void init() {
      declare(FinalState(), "LabFS");
      declare(PromptFinalState(Cuts::pid == -11), "PromptPositrons");
      _xEdges = {0.023,0.040,0.055,0.075,0.100,0.140,0.200,0.300,0.400,0.600};
      _coarseXEdges = {0.023,0.055,0.100,0.600};
      _threeXEdges = {0.023,0.040,0.055,0.075,0.100,0.140,0.200,0.300,0.400,0.600};
      _zEdges = {0.100,0.200,0.300,0.400,0.500,0.600,0.700,0.800};
      _ptEdges = {0.0,0.15,0.30,0.40,0.50,0.60,2.0};
      _threeZEdges = {0.20,0.35,0.50,0.80};
      _threePtEdges = {0.0,0.30,0.50,2.0};

      bookSpecies("proton_piplus");
      bookSpecies("proton_piminus");
      for (const std::string species :
           {"piplus","piminus","kplus","kminus"}) {
        bookSpecies("deuteron_" + species);
      }

      // The unidentified-hadron yields are used only for the published
      // deuteron charge-difference observable.
      book(_hDeuteronHPlusX, "Yield_deuteron_hplus_x", _xEdges);
      book(_hDeuteronHMinusX, "Yield_deuteron_hminus_x", _xEdges);
      _scaled.push_back(_hDeuteronHPlusX);
      _scaled.push_back(_hDeuteronHMinusX);

      book(_covProtonPi, "CovarianceProxy_proton_pi_charge_difference", _xEdges);
      book(_covDeuteronPi, "CovarianceProxy_deuteron_pi_charge_difference", _xEdges);
      book(_covDeuteronK, "CovarianceProxy_deuteron_k_charge_difference", _xEdges);
      book(_covDeuteronH, "CovarianceProxy_deuteron_h_charge_difference", _xEdges);
      _scaled.insert(_scaled.end(), {
        _covProtonPi, _covDeuteronPi, _covDeuteronK, _covDeuteronH});

      book(_acceptedX, "Accepted_X", _xEdges);
      book(_acceptedQ2, "Accepted_Q2", logspace(40, 1.0, 30.0));
      book(_acceptedY, "Accepted_Y", 34, 0.0, 0.85);
      book(_acceptedW2, "Accepted_W2", 40, 10.0, 60.0);
      book(_acceptedTheta, "Accepted_Theta", 36, 0.04, 0.22);
      book(_acceptedZ, "Accepted_Z", 35, 0.10, 0.80);
      book(_acceptedPt, "Accepted_PhT", 40, 0.0, 2.0);
      book(_acceptedXF, "Accepted_XF", 36, 0.10, 1.0);
      _scaled.insert(_scaled.end(), {
        _acceptedX, _acceptedQ2, _acceptedY, _acceptedW2, _acceptedTheta,
        _acceptedZ, _acceptedPt, _acceptedXF});
    }

    void analyze(const Event& event) {
      const Particles finalState = apply<FinalState>(event, "LabFS").particles();
      const Particles promptPositrons =
        apply<PromptFinalState>(event, "PromptPositrons").particles();
      const DISKinematics dis = reconstruct(event, promptPositrons);
      if (!dis.valid) vetoEvent;
      if (dis.Q2 <= 1.0 || dis.W2 <= 10.0 || dis.y >= 0.85) vetoEvent;
      if (dis.theta <= 0.04 || dis.theta >= 0.22) vetoEvent;
      if (binIndex(dis.x, _xEdges) < 0) vetoEvent;

      _acceptedX->fill(dis.x);
      _acceptedQ2->fill(dis.Q2);
      _acceptedY->fill(dis.y);
      _acceptedW2->fill(dis.W2);
      _acceptedTheta->fill(dis.theta);

      std::map<std::string, double> xCounts;
      std::map<std::string, EventSpecies> eventSpecies;
      for (const Particle& particle : finalState) {
        const int absPid = std::abs(particle.pid());
        const bool identified = absPid == 211 || absPid == 321;
        const bool unidentifiedChargedHadron =
          PID::isHadron(particle.pid()) && particle.charge3() != 0;
        if (!identified && !unidentifiedChargedHadron) continue;

        const HadronKinematics hadron = hadronKinematics(dis, particle);
        if (!hadron.valid || hadron.xF <= 0.10 ||
            hadron.z <= 0.10 || hadron.z >= 0.80) continue;

        _acceptedZ->fill(hadron.z);
        _acceptedPt->fill(hadron.pt);
        _acceptedXF->fill(hadron.xF);

        // A proton event supplies both the physical hydrogen prediction and
        // the proton component entering the deuteron impulse approximation.
        if (dis.targetPid == 2212 && absPid == 211 &&
            hadron.momentum > 4.0 && hadron.momentum < 13.8) {
          const std::string prefix =
            particle.pid() > 0 ? "proton_piplus" : "proton_piminus";
          accumulateSpecies(eventSpecies[prefix], dis.x, hadron);
          if (hadron.z > 0.20) xCounts[prefix] += 1.0;
        }

        // Both proton and neutron target samples contribute to deuterium,
        // with the common deuterium acceptance applied independently.
        if (identified && hadron.momentum > 2.0 && hadron.momentum < 15.0) {
          std::string species;
          if (particle.pid() ==  211) species = "piplus";
          if (particle.pid() == -211) species = "piminus";
          if (particle.pid() ==  321) species = "kplus";
          if (particle.pid() == -321) species = "kminus";
          const std::string prefix = "deuteron_" + species;
          accumulateSpecies(eventSpecies[prefix], dis.x, hadron);
          if (hadron.z > 0.20) xCounts[prefix] += 1.0;
        }

        // The charge-unidentified result uses the relaxed p_h > 0.5 GeV
        // requirement stated in the publication.
        if (unidentifiedChargedHadron && hadron.momentum > 0.5 &&
            hadron.momentum < 15.0 && hadron.z > 0.20) {
          if (particle.charge3() > 0) {
            xCounts["deuteron_hplus"] += 1.0;
          } else {
            xCounts["deuteron_hminus"] += 1.0;
          }
        }
      }

      for (const auto& entry : eventSpecies)
        flushSpecies(entry.first, entry.second);
      if (xCounts["deuteron_hplus"] > 0.0)
        _hDeuteronHPlusX->fill(dis.x, xCounts["deuteron_hplus"]);
      if (xCounts["deuteron_hminus"] > 0.0)
        _hDeuteronHMinusX->fill(dis.x, xCounts["deuteron_hminus"]);

      fillChargeCovariance(
        dis.x, xCounts["proton_piplus"], xCounts["proton_piminus"], _covProtonPi);
      fillChargeCovariance(
        dis.x, xCounts["deuteron_piplus"], xCounts["deuteron_piminus"], _covDeuteronPi);
      fillChargeCovariance(
        dis.x, xCounts["deuteron_kplus"], xCounts["deuteron_kminus"], _covDeuteronK);
      fillChargeCovariance(
        dis.x, xCounts["deuteron_hplus"], xCounts["deuteron_hminus"], _covDeuteronH);
    }

    void finalize() {
      if (sumW() == 0.0) return;
      const double factor = crossSection()/picobarn/sumW();
      for (Histo1DPtr& hist : _scaled) scale(hist, factor);
    }

  private:
    static Particle scatteredLepton(
        const Particles& finalState, const Particle& incoming) {
      Particles candidates;
      for (const Particle& particle : finalState) {
        if (particle.pid() == incoming.pid()) candidates.push_back(particle);
      }
      if (candidates.empty()) return Particle();
      return *std::max_element(
        candidates.begin(), candidates.end(),
        [](const Particle& first, const Particle& second) {
          return first.E() < second.E();
        });
    }

    static DISKinematics reconstruct(
        const Event& event, const Particles& finalState) {
      DISKinematics out;
      const ParticlePair incoming = Rivet::beams(event);
      Particle lepton, target;
      if (PID::isLepton(incoming.first.pid()) &&
          PID::isHadron(incoming.second.pid())) {
        lepton = incoming.first;
        target = incoming.second;
      } else if (PID::isLepton(incoming.second.pid()) &&
                 PID::isHadron(incoming.first.pid())) {
        lepton = incoming.second;
        target = incoming.first;
      } else {
        return out;
      }
      if (lepton.pid() != -11 ||
          (target.pid() != 2212 && target.pid() != 2112)) return out;
      const Particle outgoing = scatteredLepton(finalState, lepton);
      if (outgoing.pid() == PID::ANY) return out;

      out.k = lepton.momentum();
      out.kp = outgoing.momentum();
      out.target = target.momentum();
      out.q = out.k-out.kp;
      const double pdotq = out.target*out.q;
      const double pdotk = out.target*out.k;
      if (pdotq <= 0.0 || pdotk <= 0.0) return out;
      out.Q2 = -out.q.mass2()/GeV2;
      out.x = -out.q.mass2()/(2.0*pdotq);
      out.y = pdotq/pdotk;
      out.W2 = (out.target+out.q).mass2()/GeV2;
      out.theta = out.k.angle(out.kp);
      out.targetPid = target.pid();
      out.valid =
        std::isfinite(out.Q2) && std::isfinite(out.x) &&
        std::isfinite(out.y) && std::isfinite(out.W2) &&
        std::isfinite(out.theta);
      return out;
    }

    static HadronKinematics hadronKinematics(
        const DISKinematics& dis, const Particle& particle) {
      HadronKinematics out;
      const FourMomentum momentum = particle.momentum();
      const double pdotq = dis.target*dis.q;
      if (pdotq <= 0.0 || dis.q.p3().mod() == 0.0) return out;
      out.z = (dis.target*momentum)/pdotq;
      const Vector3 qhat = dis.q.p3().unit();
      const double longitudinal = momentum.p3().dot(qhat);
      const double pt2 = momentum.p3().mod2()-longitudinal*longitudinal;
      out.pt = std::sqrt(std::max(0.0, pt2))/GeV;
      out.momentum = momentum.p3().mod()/GeV;

      LorentzTransform boost;
      const FourMomentum hadronicSystem = dis.target+dis.q;
      boost.setBetaVec(-hadronicSystem.betaVec());
      const FourMomentum qcm = boost.transform(dis.q);
      const FourMomentum hcm = boost.transform(momentum);
      if (qcm.p3().mod() == 0.0 || hadronicSystem.mass() <= 0.0) return out;
      out.xF = 2.0*hcm.p3().dot(qcm.p3().unit())/hadronicSystem.mass();

      const Vector3 leptonNormal = dis.k.p3().cross(dis.kp.p3());
      const Vector3 hadronNormal = dis.q.p3().cross(momentum.p3());
      if (leptonNormal.mod() > 0.0 && hadronNormal.mod() > 0.0) {
        out.cosPhi = std::max(
          -1.0, std::min(1.0, leptonNormal.unit().dot(hadronNormal.unit())));
      }
      out.valid =
        std::isfinite(out.z) && std::isfinite(out.pt) &&
        std::isfinite(out.xF) && std::isfinite(out.cosPhi);
      return out;
    }

    static int binIndex(double value, const std::vector<double>& edges) {
      if (value < edges.front() || value >= edges.back()) return -1;
      const auto found = std::upper_bound(edges.begin(), edges.end(), value);
      return int(found-edges.begin())-1;
    }

    void bookSpecies(const std::string& prefix) {
      SpeciesHistograms& hist = _species[prefix];
      book(hist.x, "Yield_" + prefix + "_x", _xEdges);
      book(hist.xz, "Yield_" + prefix + "_xz", 21, 0.0, 21.0);
      book(hist.xpt, "Yield_" + prefix + "_xpt", 18, 0.0, 18.0);
      book(hist.xzpt, "Yield_" + prefix + "_xzpt", 81, 0.0, 81.0);
      book(hist.moment, "CosPhiNumerator_" + prefix, 30, 0.0, 30.0);
      book(hist.momentDenominator, "CosPhiDenominator_" + prefix, 30, 0.0, 30.0);
      book(hist.momentCovariancePositive,
           "CosPhiCovariancePositive_" + prefix, 30, 0.0, 30.0);
      book(hist.momentCovarianceNegative,
           "CosPhiCovarianceNegative_" + prefix, 30, 0.0, 30.0);
      _scaled.insert(_scaled.end(), {
        hist.x, hist.xz, hist.xpt, hist.xzpt, hist.moment,
        hist.momentDenominator, hist.momentCovariancePositive,
        hist.momentCovarianceNegative});
    }

    void accumulateSpecies(
        EventSpecies& counts, double x, const HadronKinematics& hadron) const {
      const int ix = binIndex(x, _xEdges);
      const int icx = binIndex(x, _coarseXEdges);
      if (ix < 0 || icx < 0) return;

      if (hadron.z > 0.20) counts.x[ix] += 1.0;
      const int iz = binIndex(hadron.z, _zEdges);
      if (iz >= 0) counts.xz[icx*7+iz] += 1.0;
      if (hadron.z > 0.20) {
        const int ipt = binIndex(hadron.pt, _ptEdges);
        if (ipt >= 0) counts.xpt[icx*6+ipt] += 1.0;
        const int ix3 = binIndex(x, _threeXEdges);
        const int iz3 = binIndex(hadron.z, _threeZEdges);
        const int ipt3 = binIndex(hadron.pt, _threePtEdges);
        if (ix3 >= 0 && iz3 >= 0 && ipt3 >= 0)
          counts.xzpt[ix3*9+iz3*3+ipt3] += 1.0;
      }

      if (hadron.z > 0.20) {
        for (const int momentBin :
             azimuthalMomentBins(x, hadron.z, hadron.pt)) {
          counts.moment[momentBin] += 2.0*hadron.cosPhi;
          counts.momentDenominator[momentBin] += 1.0;
        }
      }
    }

    void flushSpecies(
        const std::string& prefix, const EventSpecies& counts) {
      SpeciesHistograms& hist = _species.at(prefix);
      fillEventBins(hist.x, counts.x, false, _xEdges);
      fillEventBins(hist.xz, counts.xz);
      fillEventBins(hist.xpt, counts.xpt);
      fillEventBins(hist.xzpt, counts.xzpt);
      fillEventBins(hist.moment, counts.moment);
      fillEventBins(hist.momentDenominator, counts.momentDenominator);
      for (size_t index = 0; index < counts.moment.size(); ++index) {
        const double covariance =
          counts.moment[index]*counts.momentDenominator[index];
        if (covariance > 0.0)
          hist.momentCovariancePositive->fill(
            index+0.5, std::sqrt(covariance));
        else if (covariance < 0.0)
          hist.momentCovarianceNegative->fill(
            index+0.5, std::sqrt(-covariance));
      }
    }

    static void fillEventBins(
        const Histo1DPtr& hist, const std::vector<double>& values,
        bool flat = true, const std::vector<double>& edges = {}) {
      for (size_t index = 0; index < values.size(); ++index) {
        if (values[index] == 0.0) continue;
        const double coordinate = flat
          ? index+0.5 : 0.5*(edges[index]+edges[index+1]);
        hist->fill(coordinate, values[index]);
      }
    }

    static std::vector<int> azimuthalMomentBins(
        double x, double z, double pt) {
      const std::vector<double> x5 = {0.023,0.040,0.055,0.075,0.140,0.600};
      const std::vector<double> x2 = {0.023,0.100,0.600};
      const std::vector<double> z2 = {0.20,0.40,0.80};
      const std::vector<double> z5 = {0.20,0.32,0.44,0.56,0.68,0.80};
      const std::vector<double> pt5 = {0.0,0.30,0.40,0.50,0.60,2.0};
      const int ix5 = binIndex(x, x5);
      const int ix2 = binIndex(x, x2);
      const int iz2 = binIndex(z, z2);
      const int iz5 = binIndex(z, z5);
      const int ipt5 = binIndex(pt, pt5);
      std::vector<int> output;
      if (ix5 >= 0 && iz2 >= 0) output.push_back(iz2*5+ix5);
      if (ix2 >= 0 && iz5 >= 0) output.push_back(10+ix2*5+iz5);
      if (ix2 >= 0 && ipt5 >= 0) output.push_back(20+ix2*5+ipt5);
      return output;
    }

    static void fillChargeCovariance(
        double x, double positive, double negative, const Histo1DPtr& hist) {
      if (positive <= 0.0 || negative <= 0.0) return;
      hist->fill(x, std::sqrt(positive*negative));
    }

    std::vector<double> _xEdges, _coarseXEdges, _threeXEdges;
    std::vector<double> _zEdges, _ptEdges;
    std::vector<double> _threeZEdges, _threePtEdges;
    std::map<std::string, SpeciesHistograms> _species;
    std::vector<Histo1DPtr> _scaled;
    Histo1DPtr _hDeuteronHPlusX, _hDeuteronHMinusX;
    Histo1DPtr _covProtonPi, _covDeuteronPi, _covDeuteronK, _covDeuteronH;
    Histo1DPtr _acceptedX, _acceptedQ2, _acceptedY, _acceptedW2;
    Histo1DPtr _acceptedTheta, _acceptedZ, _acceptedPt, _acceptedXF;
  };

  RIVET_DECLARE_PLUGIN(HERMES_2019_I1698889);

}
