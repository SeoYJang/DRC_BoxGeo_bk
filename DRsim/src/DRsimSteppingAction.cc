#include "DRsimSteppingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"

DRsimSteppingAction::DRsimSteppingAction(DRsimEventAction* eventAction)
: G4UserSteppingAction(), fEventAction(eventAction)
{}

DRsimSteppingAction::~DRsimSteppingAction() {}

void DRsimSteppingAction::UserSteppingAction(const G4Step* step) {
  if (step->GetTrack()->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) return;
  /*if (step->GetTrack()->GetDefinition() == G4Electron::ElectronDefinition() || step->GetTrack()->GetDefinition() == G4Positron::PositronDefinition()) 
  {
    if (step->GetPreStepPoint()->GetTotalEnergy() < 0.7 * MeV) {
      step->GetTrack()->SetTrackStatus(fStopAndKill);
      return;
    }
  }*/

  G4Track* track = step->GetTrack();
  G4ParticleDefinition* particle = track->GetDefinition();
  G4int pdgID = particle->GetPDGEncoding();

  G4StepPoint* presteppoint = step->GetPreStepPoint();
  G4StepPoint* poststeppoint = step->GetPostStepPoint();
  G4LogicalVolume* preVol = presteppoint->GetPhysicalVolume()->GetLogicalVolume();
  G4TouchableHandle theTouchable = presteppoint->GetTouchableHandle();

  if (poststeppoint->GetStepStatus() == fWorldBoundary) {
    fLeak.E = track->GetTotalEnergy();
    fLeak.px = track->GetMomentum().x();
    fLeak.py = track->GetMomentum().y();
    fLeak.pz = track->GetMomentum().z();
    fLeak.vx = presteppoint->GetPosition().x();
    fLeak.vy = presteppoint->GetPosition().y();
    fLeak.vz = presteppoint->GetPosition().z();
    fLeak.vt = presteppoint->GetGlobalTime();
    fLeak.pdgId = track->GetDefinition()->GetPDGEncoding();

    fEventAction->fillLeaks(fLeak);
  }

  G4String matName = preVol->GetMaterial()->GetName();

  if ( matName=="G4_Galactic" || matName=="Air" ) return;

  G4VPhysicalVolume* motherTower = GetMotherTower(theTouchable);
  fEdep.ModuleNum = GetModuleNum(motherTower->GetName());

  G4double pdgCharge = particle->GetPDGCharge();

  fEdep.Edep = step->GetTotalEnergyDeposit();
  fEdep.EdepEle = (std::abs(pdgID)==11) ? fEdep.Edep : 0.;
  fEdep.EdepGamma = (std::abs(pdgID)==22) ? fEdep.Edep : 0.;
  fEdep.EdepCharged = ( std::round(std::abs(pdgCharge)) != 0. ) ? fEdep.Edep : 0.;
  /*fEdep.vx = presteppoint->GetPosition().x();
  fEdep.vy = presteppoint->GetPosition().y();
  fEdep.vz = presteppoint->GetPosition().z();

  fEdep.EdepModule = (matName=="Copper") ? fEdep.Edep : 0.;
  fEdep.EdepCeren = (matName=="FluorinatedPolymer") ? fEdep.Edep : 0.;
  fEdep.EdepScint = (matName=="Polystyrene") ? fEdep.Edep : 0.;

  if ( matName=="PMMA") {
    fEdep.EdepScint = (preVol->GetNoDaughters()>0) ? fEdep.Edep : 0.;
    fEdep.EdepCeren = (preVol->GetNoDaughters()>0) ? 0. : fEdep.Edep;
  }*/

  if ( fEdep.Edep > 0. ) {
    fEventAction->fillEdeps(fEdep);
    //G4cout<<fEdep.ModuleNum<<G4endl;
  //  fEventAction->fillEdepPoss(fEdep);
  }

  return;
}
