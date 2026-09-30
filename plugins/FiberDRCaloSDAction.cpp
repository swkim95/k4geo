// DD4hep Framework include files
#include <cmath>

#include "DD4hep/Segmentations.h"
#include "DD4hep/detail/SegmentationsInterna.h"
#include "DDG4/Geant4Mapping.h"
#include "DDG4/Geant4Random.h"
#include "DDG4/Geant4SensDetAction.inl"

#include "G4OpticalPhoton.hh"
// k4geo Framework include files

#include "DRCaloFastSimModel.h"
#include "FiberDRCaloSDAction.h"
#include "detectorSegmentations/GridDRcalo_k4geo.h"

// Geant4 include files
#include "G4HCofThisEvent.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpProcessSubType.hh"
#include "G4OpticalPhoton.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Poisson.hh"
#include "G4PhysicalConstants.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "G4TouchableHistory.hh"

#include <vector>

/// Namespace for the AIDA detector description toolkit
namespace dd4hep {

/// Namespace for the Geant4 based simulation part of the AIDA detector description toolkit
namespace sim {

  /*
   *  Geant4SensitiveAction<DRcaloSDAction> sensitive detector for the Dual-readout calorimeter
   *
   *  \author  Sungwon Kim, Sanghyun Ko
   *  \version 1.0
   *  \ingroup DD4HEP_SIMULATION
   */

  struct DRCData {
  public:
    bool skipScint = true;
    DRCFiberModel fastfiber;

    G4int fWavBin;
    G4int fTimeBin;
    G4float fWavlenStart;
    G4float fWavlenEnd;
    G4float fTimeStart;
    G4float fTimeEnd;
    G4float fWavlenStep;
    G4float fTimeStep;

    // --- Stage 4 (o3 chi recovery plan, 2026-09-26): Birks quenching (S ch) + SiPM
    // PDE + yellow filter (both channels), applied at the Sim level instead of the
    // Digi level. All three gated by declareProperty switches, DEFAULT OFF so any
    // steering file that does not explicitly set them reproduces the exact
    // pre-existing (noOpt baseline) behavior -- see FiberDRCaloSDAction.cpp's top
    // comment / the plan's Stage 4 section for the physics motivation and the
    // Stage 2 (Digi-level) finding that most of the chi improvement traces to
    // SiPM PDE + filter narrowing the (h/e)_C anomaly, not Birks.
    bool applyBirks = false;    // S ch only; Cerenkov has no Birks quenching (no scintillation)
    bool applySiPMEff = false;  // both channels: C ch per-photon PDE Bernoulli test,
                                 // S ch folded into the effective factor below
    bool applyFilter = false;   // S ch only (yellow filter blocks short-wavelength light
                                 // reaching the SiPM); folded into the effective factor
    double kBirks = 0.126;      // mm/MeV, DR_Polystyrene XML value (generic polystyrene
                                 // literature constant, NOT RD52-measured -- see plan)
    double scintYieldPerKeV = 13.9; // Kuraray SCSF-78 fiber, matches XML's SCINTILLATIONYIELD

    // Wavelength grids + response curves, all copied verbatim from
    // example/condor/common/run_digi_reco.py's SimulateSiPMwithEdep/
    // SimulateSiPMwithOpticalPhoton property values (kept in the source's own
    // descending-wavelength order; converted to ascending-eV once in
    // initEffectiveFactors()). Cross-referencing those exact numbers means Stage 4
    // (Sim-level) and Stage 5 (re-Digitized Stage-4 output, eff=1) use internally
    // consistent physics.
    static inline const std::vector<double> fWavlen_S_nm = { // 31 points, S ch (Edep-mode digitizer)
        900., 850., 800., 750., 725., 700., 675., 650., 625., 600.,
        590., 580., 570., 560., 550., 540., 530., 520., 510., 500.,
        490., 480., 470., 460., 450., 440., 430., 420., 400., 350., 300.
    };
    static inline const std::vector<double> fScintSpectrum = { // Kuraray SCSF-78 emission spectrum
        0., 0., 0., 0., 0., 0., 0., 0.0003, 0.0008, 0.0032,
        0.0057, 0.0084, 0.0153, 0.0234, 0.0343, 0.0604, 0.0927, 0.1398, 0.2105, 0.2903,
        0.4122, 0.5518, 0.7086, 0.8678, 1., 0.8676, 0.2311, 0.0033, 0.0012, 0., 0.
    };
    static inline const std::vector<double> fFilterEff = { // Kodak Wratten No.9 yellow filter
        0.903, 0.903, 0.903, 0.903, 0.903, 0.903, 0.902, 0.901, 0.898, 0.895,
        0.893, 0.891, 0.888, 0.883, 0.87, 0.838, 0.76, 0.62, 0.488, 0.345,
        0.207, 0.083, 0.018, 0., 0., 0., 0., 0., 0., 0., 0.
    };
    static inline const std::vector<double> fPDE_S = { // Hamamatsu S14160-1310PS
        0.02, 0.025, 0.045, 0.06, 0.0675, 0.075, 0.0925, 0.11, 0.125, 0.14,
        0.146, 0.152, 0.158, 0.164, 0.17, 0.173, 0.176, 0.178, 0.179, 0.18,
        0.181, 0.182, 0.183, 0.184, 0.18, 0.173, 0.166, 0.158, 0.15, 0.12, 0.05
    };
    static inline const std::vector<double> fWavlen_C_nm = { // 32 points, C ch (OpticalPhoton-mode digitizer)
        900., 850., 800., 750., 725., 700., 675., 650., 625., 600.,
        590., 580., 570., 560., 550., 540., 530., 520., 510., 500.,
        490., 480., 470., 460., 450., 440., 430., 420., 400., 350., 300., 280.
    };
    static inline const std::vector<double> fPDE_C = { // Hamamatsu S14160-1315PS
        0.03, 0.05, 0.07, 0.1, 0.13, 0.14, 0.16, 0.19, 0.21, 0.23,
        0.238, 0.246, 0.254, 0.262, 0.27, 0.28, 0.29, 0.3, 0.31, 0.32,
        0.322, 0.324, 0.326, 0.328, 0.33, 0.33, 0.32, 0.30, 0.27, 0.22, 0.12, 0.
    };

    std::vector<double> fEnergy_S_eV; // fWavlen_S_nm converted to eV, ascending
    std::vector<double> fEnergy_C_eV; // fWavlen_C_nm converted to eV, ascending
    double fEffS = 1.0;         // precomputed spectrum-weighted <filterEff * PDE_S> for the S ch
    bool fEffSComputed = false; // initEffectiveFactors() is lazy (see call site in process()):
                                 // declareProperty() in the C++ ctor only registers the default;
                                 // any steering-file override of applyFilter/applySiPMEff is applied
                                 // later (before the event loop starts, but after construction), so
                                 // computing fEffS eagerly in the ctor would silently ignore it.

    // Linear interpolation of a tabulated (ascending-energy) response curve at
    // photon energy `en` -- same technique as the pre-124aeddd o3 SD action's
    // getFilterEff() (git show 124aeddd~1:plugins/FiberDRCaloSDAction.cpp).
    static double interpolate(const std::vector<double>& xEnergyEV, const std::vector<double>& y, double en) {
      if (en <= xEnergyEV.front())
        return y.front();
      if (en >= xEnergyEV.back())
        return y.back();
      for (std::size_t i = 1; i < xEnergyEV.size(); ++i) {
        if (en <= xEnergyEV[i]) {
          double x1 = xEnergyEV[i - 1], x2 = xEnergyEV[i];
          double y1 = y[i - 1], y2 = y[i];
          return y1 + (y2 - y1) / (x2 - x1) * (en - x1);
        }
      }
      return y.back();
    }

    // Birks quenching: dE_quenched = dE / (1 + kB * dE/dx), dE/dx computed naively
    // as (this step's deposit) / (this step's length) -- same formula as the o2
    // reference (plugins/DRTubesSglHpr.hh's ApplyBirks(), algebraically identical).
    //
    // CAP on dE/dx (2026-09-27, Stage 6 postmortem): a first version without this
    // cap gave a ~96% reduction in total S-ch energy for a 40 GeV e- shower --
    // implying an EVENT-AVERAGE dE/dx of ~20000 MeV/cm, far beyond any real charged
    // particle in polystyrene (even a near-Bragg-peak alpha is ~100s of MeV/cm).
    // This is a well-known Geant4 Birks-implementation pitfall: G4's stepping
    // algorithm (multiple scattering, boundary crossings in this detector's thin,
    // 0.485mm-radius fiber core) produces many steps whose length is an artifact of
    // step-size control, not of the true local energy-loss rate, so a raw
    // edep/stepLength ratio can spike to unphysical values for those steps even
    // though edep itself is unremarkable. The physically-correct fix is Geant4's
    // G4EmSaturation helper (uses the track's restricted dE/dx from the EM physics
    // tables, not a raw per-step ratio) -- not used here for time reasons. Instead,
    // dE/dx is capped at kMaxDEdxMeVPerMm (chosen as a generous but bounded ceiling,
    // comparable to a very slow proton/alpha near its Bragg peak in plastic, i.e.
    // already an extreme case for any electron/positron/relativistic hadron
    // fragment), which caps the worst-case single-step quenching at
    // 1/(1+kBirks*kMaxDEdxMeVPerMm) while leaving the vast majority of steps
    // (genuinely at or below this dE/dx) completely unaffected.
    static constexpr double kMaxDEdxMeVPerMm = 50.0; // ~5000 MeV/cm
    static double applyBirksQuenching(double dE_GeV, double stepLength_G4mm, double kBirks) {
      double stepLength_mm = stepLength_G4mm / CLHEP::mm;
      if (stepLength_mm <= 0.)
        return dE_GeV;
      double dEdx_MeVPerMm = std::min((dE_GeV * 1.e3) / stepLength_mm, kMaxDEdxMeVPerMm); // GeV -> MeV
      return dE_GeV / (1. + kBirks * dEdx_MeVPerMm);
    }

    // Precompute fEnergy_S_eV/fEnergy_C_eV and fEffS. Called once from the
    // Geant4SensitiveAction constructor, after the declareProperty defaults (or
    // steering-file overrides) for applyFilter/applySiPMEff are already set.
    void initEffectiveFactors() {
      auto nmToEv = [](double nm) { return 1239.84193 / nm; }; // hc in eV*nm
      fEnergy_S_eV.clear();
      for (double nm : fWavlen_S_nm)
        fEnergy_S_eV.push_back(nmToEv(nm)); // descending nm -> ascending eV
      fEnergy_C_eV.clear();
      for (double nm : fWavlen_C_nm)
        fEnergy_C_eV.push_back(nmToEv(nm));

      // Discrete spectrum-weighted average (NOT a true wavelength-bin-width-weighted
      // integral -- the source arrays' bins are non-uniform, denser near the
      // scintillation peak; treated as an approximation, documented in the plan) of
      // filterEff * PDE_S over the Kuraray SCSF-78 emission spectrum:
      //   effS = sum(scintSpectrum_i * filterEff_i * PDE_S_i) / sum(scintSpectrum_i)
      double num = 0.0, den = 0.0;
      for (std::size_t i = 0; i < fScintSpectrum.size(); ++i) {
        double filt = applyFilter ? fFilterEff[i] : 1.0;
        double pde = applySiPMEff ? fPDE_S[i] : 1.0;
        num += fScintSpectrum[i] * filt * pde;
        den += fScintSpectrum[i];
      }
      fEffS = (den > 0.0) ? (num / den) : 1.0;
      fEffSComputed = true;
    }

  public:
    G4double wavToE(G4double wav) const { return CLHEP::h_Planck * CLHEP::c_light / wav; }

    int findWavBin(G4double en) const {
      int i = 0;
      for (; i < fWavBin + 1; i++) {
        if (en < wavToE((fWavlenStart - static_cast<float>(i) * fWavlenStep) * CLHEP::nm))
          break;
      }

      return fWavBin + 1 - i;
    }

    int findTimeBin(G4double stepTime) const {
      int i = 0;
      for (; i < fTimeBin + 1; i++) {
        if (stepTime < ((fTimeStart + static_cast<float>(i) * fTimeStep) * CLHEP::ns))
          break;
      }

      return i;
    }

    // default constructor
    DRCData() : fWavBin(120), fTimeBin(2000), fWavlenStart(900.), fWavlenEnd(300.), fTimeStart(-100.), fTimeEnd(100.) {
      fWavlenStep = (fWavlenStart - fWavlenEnd) / (float)fWavBin;
      fTimeStep = (fTimeEnd - fTimeStart) / (float)fTimeBin;
    }
  }; // struct DRCData

  template <>
  Geant4SensitiveAction<DRCData>::Geant4SensitiveAction(Geant4Context* ctxt, const std::string& nam, DetElement det,
                                                        Detector& desc)
      : Geant4Sensitive(ctxt, nam, det, desc), m_collectionName(), m_collectionID(0) {
    declareProperty("skipScint", m_userData.skipScint = true);
    declareProperty("ReadoutName", m_readoutName);
    declareProperty("CollectionName", m_collectionName);
    declareProperty("timeStart", m_userData.fTimeStart); // in ns
    declareProperty("timeEnd", m_userData.fTimeEnd);     // in ns
    // delegate to user's responsibility to ensure timeStep divides (timeEnd - timeStart)
    declareProperty("timeStep", m_userData.fTimeStep); // in ns
    declareProperty("timeBin", m_userData.fTimeBin);   // in ns
    declareProperty("verbose", m_userData.fastfiber.fVerbose);
    declareProperty("safety", m_userData.fastfiber.fSafety);
    // Stage 4 (o3 chi recovery plan): all default OFF, so any steering file that
    // does not set these reproduces the exact pre-existing behavior.
    declareProperty("applyBirks", m_userData.applyBirks);
    declareProperty("applySiPMEff", m_userData.applySiPMEff);
    declareProperty("applyFilter", m_userData.applyFilter);
    declareProperty("kBirks", m_userData.kBirks);
    declareProperty("scintYieldPerKeV", m_userData.scintYieldPerKeV);
    InstanceCount::increment(this);

    m_hitCreationMode = HitCreationFlags::DETAILED_MODE; // always store step position
  }

  /// Define collections created by this sensitive action object
  template <>
  void Geant4SensitiveAction<DRCData>::defineCollections() {
    std::string readout_name = m_sensitive.readout().name();
    m_collectionID = defineCollection<Geant4DRCalorimeter::Hit>(m_sensitive.readout().name());
    std::cout << "defineCollection Geant4DRCalorimeter readout_name   : " << readout_name << std::endl;
    std::cout << "defineCollection Geant4DRCalorimeter m_collectionID : " << m_collectionID << std::endl;

    if (m_userData.skipScint) {
      defineCollection<Geant4Calorimeter::Hit>(std::string(m_sensitive.readout().name()) + "_scint");
      std::cout << "defineCollection Geant4Calorimeter readout_name   : " << readout_name + "_scint" << std::endl;
      std::cout << "defineCollection Geant4Calorimeter m_collectionID : " << m_collectionID + 1 << std::endl;
    }
  }

  /// Method for generating hit(s) using the information of G4Step object.
  template <>
  G4bool Geant4SensitiveAction<DRCData>::process(G4Step const* step, G4TouchableHistory*) {
    // Lazy one-time init (Stage 4): see fEffSComputed's comment for why this can't
    // run in the ctor -- by the time the first Step is processed, any steering-file
    // override of applyFilter/applySiPMEff is guaranteed to already be in effect.
    if (!m_userData.fEffSComputed)
      m_userData.initEffectiveFactors();

    // optical photons traveling through the cladding is not interesting
    // but consume lots of computation power
    // let's kill optical photons inside the cladding whose status is not StepTooSmall
    if (step->GetPreStepPoint() &&
        step->GetPreStepPoint()->GetPhysicalVolume()->GetLogicalVolume()->GetNoDaughters() == 1 &&
        step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
      // 1e-9 is the default tolerance
      // for the warnings, see
      // https://geant4-forum.web.cern.ch/t/error-occurs-when-an-optical-photon-hits-the-edge-of-a-cubic-scintillator/8748
      // we're not interested in the optical photon position
      // but only the number of photons (and timing)
      // so we assume it's fine to ignore the warnings
      if (step->GetStepLength() > 1e-8 * CLHEP::mm)
        step->GetTrack()->SetTrackStatus(fStopAndKill);
    }
    // the fiber itself is the SD
    // we skip the scintillation process and only account for the Cherenkov process
    // remember to turn off the scintillation process!
    if (m_userData.skipScint) {
      // we need to move the touchable to the tower to retrieve volID
      auto* touchable = const_cast<G4VTouchable*>(step->GetPostStepPoint()->GetTouchable());
      const auto* logicalVol = touchable->GetVolume()->GetLogicalVolume();

      // we're not interested in the world or assembly volume
      if (touchable->GetHistoryDepth() < 2)
        return false;

      // we're only interested in the fiber core for optical photons
      if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition() &&
          logicalVol->GetNoDaughters() != 0)
        return false;

      // now let's make the touchable points to the tower
      // world -> assembly -> tower
      touchable->MoveUpHistory(touchable->GetHistoryDepth() - 2);

      // now the touchable is the tower
      // get local position and volumeID
      // dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
      // dd4hep::VolumeID volID = volMgr.volumeID(touchable);
      // note the above method started to throw warnings since the PR DD4hep#1390
      // so we switch to the copy number

      // trick: we retrieve the tower's volID by using the fact that
      // in the DRconstructor.cpp, the first 32bits of the volID
      // contain the information relevant to the tower
      // and this 32bit integer becomes the tower's copy number
      auto copyNo = touchable->GetCopyNumber();

      // ideally the below operation should be GridDRcalo_k4geo::convertFirst32to64
      // but this is trivial, so let's avoid doing dynamic cast
      dd4hep::VolumeID volID = static_cast<dd4hep::VolumeID>(copyNo);
      G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector local = touchable->GetHistory()->GetTopTransform().TransformPoint(global);

      // convert G4 position to dd4hep position
      dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter,
                           local.y() * dd4hep::millimeter / CLHEP::millimeter,
                           local.z() * dd4hep::millimeter / CLHEP::millimeter);
      dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                            global.y() * dd4hep::millimeter / CLHEP::millimeter,
                            global.z() * dd4hep::millimeter / CLHEP::millimeter);

      // Guard against the y-bit cellID overflow crash (KNOWN_ISSUES_o3_DRC.md #3):
      // a shower secondary can escape this tower's front face into the solenoid;
      // MoveUpHistory() above still walks up to *some* tower ancestor and reports a
      // volID for it, but the step's actual local position (loc) then falls well
      // outside that tower's real x/y extent. Encoding such a position would make
      // GridDRcalo_k4geo::cellID() compute an out-of-range x/y index and crash
      // Geant4Kernel ("BitFieldElement 'y': out of range") when the decoder tries to
      // set it. Check the tower's known front-face bounds first and skip the step
      // (no hit) instead of crashing.
      {
        auto* gridSeg = dynamic_cast<dd4hep::DDSegmentation::GridDRcalo_k4geo*>(m_segmentation->segmentation);
        if (gridSeg) {
          int numx = gridSeg->numX(volID);
          int numy = gridSeg->numY(volID);
          double gridSize = gridSeg->getGridSize();
          double halfX = 0.5 * numx * gridSize;
          double halfY = 0.5 * numy * gridSize;
          if (std::abs(loc.x()) > halfX || std::abs(loc.y()) > halfY)
            return false;
        }
      }

      // retrieve cellID and distinguish Cherenkov & scintillation fibers
      auto cID = m_segmentation->cellID(loc, glob, volID);
      auto ceren = static_cast<dd4hep::DDSegmentation::VolumeID>(m_segmentation->decoder()->get(cID, "c"));
      bool IsCeren = static_cast<bool>(ceren);

      Geant4HitCollection* coll = collection(m_collectionID);

      if (IsCeren) {
        // Cherenkov fiber
        // skip anything else than optical photons
        if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
          return false;

        const auto* track = step->GetTrack();

        // reset when moving to the next track
        if (m_userData.fastfiber.mDataCurrent.trackID != track->GetTrackID())
          m_userData.fastfiber.reset();

        // need repetitive total internal reflection
        if (!m_userData.fastfiber.check_trigger(track))
          return false;

        // absorption
        if (m_userData.fastfiber.fKill) {
          step->GetTrack()->SetTrackStatus(fStopAndKill);

          return false;
        }

        // backward transportation
        if (m_userData.fastfiber.mTransportUnit < 0.) {
          step->GetTrack()->SetTrackStatus(fStopAndKill);

          return false;
        }

        // for timing measurement
        double timeUnit = m_userData.fastfiber.mDataCurrent.globalTime - m_userData.fastfiber.mDataPrevious.globalTime;
        double timeShift = timeUnit * m_userData.fastfiber.mNtransport;
        G4double energy = step->GetTrack()->GetTotalEnergy();

        // Stage 4 (o3 chi recovery plan): SiPM PDE (Hamamatsu S14160-1315PS) as a
        // per-photon Bernoulli survival test, restoring (with the updated, C-ch-
        // specific curve from run_digi_reco.py) the logic removed in commit
        // 124aeddd. No yellow filter here -- it only shapes the (broadband)
        // scintillation spectrum, not Cerenkov light.
        if (m_userData.applySiPMEff) {
          double pde = DRCData::interpolate(m_userData.fEnergy_C_eV, m_userData.fPDE_C, energy);
          Geant4Event& evt = context()->event();
          Geant4Random& rnd = evt.random();
          if (rnd.uniform(0., 1.) > pde) {
            step->GetTrack()->SetTrackStatus(fStopAndKill);
            return false;
          }
        }

        // default hit (optical photon count)
        Geant4DRCalorimeter::Hit* drHit = coll->findByKey<Geant4DRCalorimeter::Hit>(cID);

        if (!drHit) {
          drHit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
          drHit->cellID = cID;
          drHit->position = m_segmentation->position(cID) * CLHEP::mm / dd4hep::mm; // segmentation gives dd4hep unit
          drHit->SetSiPMnum(cID);
          drHit->SetTimeStart(m_userData.fTimeStart);
          drHit->SetTimeEnd(m_userData.fTimeEnd);
          drHit->SetWavlenMax(m_userData.fWavlenStart);
          drHit->SetWavlenMin(m_userData.fWavlenEnd);
          coll->add(cID, drHit);
        }

        // everything should be in the G4 unit
        // (approximate) timing at the end of the fiber
        G4double hitTime = step->GetPostStepPoint()->GetGlobalTime() + timeShift;

        drHit->photonCount();
        int wavBin = m_userData.findWavBin(energy);
        drHit->CountWavlenSpectrum(wavBin);
        int timeBin = m_userData.findTimeBin(hitTime);
        drHit->CountTimeStruct(timeBin);

        // finally kill optical photon
        step->GetTrack()->SetTrackStatus(fStopAndKill);

        return true;
      } else {
        // scintillation calo hit

        // assume nPhoton_scint >> nPhoton_cherenkov
        // kill optical photons (from Cherenkov process)
        if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
          step->GetTrack()->SetTrackStatus(fStopAndKill);

          return false;
        }

        // copy-paste of the dd4hep scintillation calorimeter SD
        //
        // IMPORTANT (2026-09-28 S-ch resolution regression fix): the second argument
        // here is dd4hep's OWN internal Birks Law switch (Geant4HitData::
        // extractContribution(step, ApplyBirksLaw) -> Geant4StepHandler::
        // doApplyBirksLaw() -> birkAttenuation() -> G4EmSaturation::
        // VisibleEnergyDeposition()), COMPLETELY INDEPENDENT of this file's own
        // m_userData.applyBirks flag (which defaults to false and gates the
        // separate applyBirksQuenching() call below). DR_Polystyrene's XML material
        // definition sets a BirksConstant property, so passing `true` here silently
        // quenched EVERY step's deposit via G4EmSaturation regardless of
        // m_userData.applyBirks, in every noOpt run, the entire time. This uncapped,
        // per-step, nonlinear saturation in a thin (0.485mm-radius) fiber core adds
        // substantial extra event-to-event variance to the S ch (dE/dx-based) sum
        // without touching the C ch (photon-counting, never goes through this code
        // path) -- confirmed to be the root cause of the S-ch-worse-than-C-ch
        // resolution regression (see chat discussion; a July 2025 local, never-
        // upstreamed fix flipped this same argument to `false`). Hardcoded to
        // `false` so dd4hep's internal Birks Law is never silently applied here;
        // m_userData.applyBirks (with its own dE/dx-capped applyBirksQuenching()
        // below) remains the only intentional Birks control for this SD action.
        //
        // Also part of the same 2026-09-28 fix: skip non-ionizing (neutral, e.g.
        // gamma) tracks and zero-length steps entirely for the S ch. These carry no
        // real local dE/dx signal but still land here (nothing upstream filters
        // them out), and letting them through creates spurious low/zero-length
        // contributions that inflate event-to-event variance of the S ch sum
        // without changing its mean much -- exactly a resolution-only regression.
        // Restores the July 2025 local (never-upstreamed) fix.
        //
        // RE-ENABLED (2026-09-29, EM/hadronic resolution campaign): the Calib-B vs
        // Calib-B-nofilter cross-check (TEST_HISTORY Parts I/J) confirmed this filter
        // is the SOLE cause of chi_fiber moving from ~0.394 (filter OFF, matches the
        // 2026-09-26 0.39/0.47 reference) to ~0.439 (filter ON) -- both are physical
        // (positive, sane S/C distributions), and the user judged the filter-ON case
        // more physical (removes non-ionizing/zero-length spurious S-ch contributions,
        // see the comment block above) going forward. This is now the baseline for
        // the EM/hadronic resolution measurement campaign and beyond.
        if (step->GetTrack()->GetDefinition()->GetPDGCharge() == 0 || step->GetStepLength() == 0.) {
          return false;
        }

        // NOTE (2026-09-28): a core-only S-ch variant (matching o2's DRTubesSDAction.cpp,
        // which unconditionally skips cladding steps before its channel-specific logic)
        // was tested here and REJECTED -- it made both S-ch resolution and p0-p3
        // uniformity measurably worse at every point (see KNOWN_ISSUES_o3_DRC.md, section
        // on the (h/e)_S investigation). Reverted back to the original, unrestricted
        // core+cladding S-ch sum below. Do not re-add the core-only filter without new
        // evidence -- that specific change has already been tried and ruled out.

        Geant4HitCollection* coll_scint = collection(m_collectionID + 1);
        Geant4Calorimeter::Hit* caloHit = coll_scint->findByKey<Geant4Calorimeter::Hit>(cID);
        // 2026-09-29 (Calib-B campaign): dd4hep-internal Birks (G4EmSaturation) set back
        // to `true`, i.e. the same Birks implementation as the 2026-09-26 0.39/0.47
        // reference. The Birks ON/OFF matrix (TEST_HISTORY Phase 1-4) showed this only
        // rescales the S ch by ~3% with no resolution change, so the long comment above
        // claiming it is the root cause of the S-ch regression is SUPERSEDED (kept for
        // history). Note: m_userData.applyBirks (custom path) has a MeV/GeV unit bug
        // (see KNOWN_ISSUES) and must stay false.
        HitContribution contrib = Geant4Calorimeter::Hit::extractContribution(step, true);

        if (!caloHit) {
          caloHit = new Geant4Calorimeter::Hit(glob);
          caloHit->cellID = cID;
          caloHit->position = m_segmentation->position(cID) * CLHEP::mm / dd4hep::mm; // segmentation gives dd4hep unit
          coll_scint->add(cID, caloHit);
        }

        caloHit->truth.emplace_back(contrib);

        // Stage 4 (o3 chi recovery plan): Birks quenching + effective SiPM PDE *
        // yellow filter factor. noOpt mode never tracks individual scintillation
        // photons (skipScint=true skips G4Scintillation entirely), so PDE/filter
        // can't be applied per-photon here the way it is for the C ch above --
        // instead, this step's Birks-corrected deposit is converted to a MEAN
        // photon count (deposit * scintYieldPerKeV * fEffS) and Poisson-sampled;
        // sum of independent Poisson draws with means mu_i ~ Poisson(sum mu_i), so
        // per-step sampling here is statistically equivalent to sampling once from
        // the hit's total accumulated mean. When ALL THREE flags are off (default),
        // this reduces to exactly the original `energyDeposit += contrib.deposit`
        // (GeV) -- units only change (GeV -> p.e. count) once applySiPMEff or
        // applyFilter is enabled.
        //
        // Stage 6 (o3 chi recovery plan, 2026-09-27): applyBirks alone (both
        // applySiPMEff and applyFilter off) is a DELIBERATELY separate branch that
        // keeps energyDeposit in GeV (just Birks-quenched), instead of converting to
        // a p.e. count. Rationale: this lets the Sim output be digitized by the
        // STOCK run_digi_reco.py unchanged -- SimulateSiPMwithEdep's own internal
        // photon generation (scintYield/scintSpectrum/filterEfficiency/
        // sipmEfficiency) expects a genuine GeV energy deposit as input, exactly
        // what it always received. Stage 5 (applyBirks+applySiPMEff+applyFilter at
        // Sim level, then bypassing Digi's own PDE/filter/yield to avoid double-
        // counting) instead fed Digi's SiPM electronics an already-p.e.-count input
        // it was never designed for, and the electronics' per-cell dark-count noise
        // floor (roughly constant per hit, ~162 counts/gate at DCR=120kHz over the
        // ~16900-pixel SiPM) ended up dominating the now-much-smaller real signal --
        // see the chat discussion for the full derivation. Birks-only avoids that
        // failure mode entirely by never touching the photon-count scale at Sim
        // level; SiPM PDE/filter are applied exactly once, in Digi, on their native
        // GeV input.
        if (m_userData.applyBirks && !m_userData.applySiPMEff && !m_userData.applyFilter) {
          caloHit->energyDeposit += DRCData::applyBirksQuenching(contrib.deposit, step->GetStepLength(),
                                                                  m_userData.kBirks);
        } else if (m_userData.applyBirks || m_userData.applySiPMEff || m_userData.applyFilter) {
          double dE_quenched_GeV = contrib.deposit;
          if (m_userData.applyBirks) {
            dE_quenched_GeV = DRCData::applyBirksQuenching(contrib.deposit, step->GetStepLength(),
                                                            m_userData.kBirks);
          }
          double nPhotonsMean =
              dE_quenched_GeV * 1.e6 * m_userData.scintYieldPerKeV * m_userData.fEffS; // GeV -> keV
          caloHit->energyDeposit += static_cast<double>(G4Poisson(nPhotonsMean));
        } else {
          caloHit->energyDeposit += contrib.deposit;
        }

        Geant4StepHandler h(step);
        mark(h.track);

        return true;
      } // !IsCeren
    } else {
      // no skipping optical photon propagation
      // SiPM wafers are SD in this case
      if (step->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition())
        return false;

      typedef Geant4DRCalorimeter::Hit Hit;

      Geant4HitCollection* coll = collection(m_collectionID);

      auto* touchable = step->GetPostStepPoint()->GetTouchable();
      dd4hep::sim::Geant4VolumeManager volMgr = dd4hep::sim::Geant4Mapping::instance().volumeManager();
      dd4hep::VolumeID volID = volMgr.volumeID(touchable);
      G4ThreeVector global = step->GetPostStepPoint()->GetPosition();
      G4ThreeVector local = touchable->GetHistory()->GetTopTransform().TransformPoint(global);

      dd4hep::Position loc(local.x() * dd4hep::millimeter / CLHEP::millimeter,
                           local.y() * dd4hep::millimeter / CLHEP::millimeter,
                           local.z() * dd4hep::millimeter / CLHEP::millimeter);
      dd4hep::Position glob(global.x() * dd4hep::millimeter / CLHEP::millimeter,
                            global.y() * dd4hep::millimeter / CLHEP::millimeter,
                            global.z() * dd4hep::millimeter / CLHEP::millimeter);

      auto cID = m_segmentation->cellID(loc, glob, volID); // This returns cID corresponding to SiPM Wafer
      Hit* hit = coll->findByKey<Hit>(cID);

      G4double hitTime = step->GetPostStepPoint()->GetGlobalTime();
      G4double energy = step->GetTrack()->GetTotalEnergy();

      if (!hit) {
        hit = new Geant4DRCalorimeter::Hit(m_userData.fWavlenStep, m_userData.fTimeStep);
        hit->cellID = cID;
        hit->position = m_segmentation->position(cID) * CLHEP::mm / dd4hep::mm; // segmentation gives dd4hep unit
        hit->SetSiPMnum(cID);
        hit->SetTimeStart(m_userData.fTimeStart);
        hit->SetTimeEnd(m_userData.fTimeEnd);
        hit->SetWavlenMax(m_userData.fWavlenStart);
        hit->SetWavlenMin(m_userData.fWavlenEnd);
        coll->add(cID, hit);
      }

      hit->photonCount();
      int wavBin = m_userData.findWavBin(energy);
      hit->CountWavlenSpectrum(wavBin);
      int timeBin = m_userData.findTimeBin(hitTime);
      hit->CountTimeStruct(timeBin);

      return true;
    } // !skipScint
  } // Geant4SensitiveAction::process

  typedef Geant4SensitiveAction<DRCData> DRCaloSDAction;
} // namespace sim
} // namespace dd4hep

#include "DDG4/Factories.h"
DECLARE_GEANT4SENSITIVE(DRCaloSDAction)
