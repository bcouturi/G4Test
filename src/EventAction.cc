//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "EventAction.hh"
#include "SensitiveDetectorHit.hh"
#include "Analysis.hh"
#include "RunAction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::EventAction(RunAction* runAction):
fSensitiveDetector_ID(-1),
fVerboseLevel(0),
fRunAction(runAction),
Edep(0.)
{
	G4cout << "### EventAction instantiated ###" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

EventAction::~EventAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::BeginOfEventAction(const G4Event*)
{
	Edep = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void EventAction::EndOfEventAction(const G4Event* aEvent)
{   
	//instantiating The Sensitive Detector Manager
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
        
    //Hit Detection System
	if (fSensitiveDetector_ID == -1) {
    	G4String SensitiveDetectorName;
        if (SDman->FindSensitiveDetector(SensitiveDetectorName="det",0)) {
            fSensitiveDetector_ID = SDman->GetCollectionID(SensitiveDetectorName="det/collection");
        }
    }
    
    SensitiveDetectorHitsCollection* fSensitiveDetectorHC = 0;
    G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
    
    if (HCE) {
        if (fSensitiveDetector_ID != -1) {
            G4VHitsCollection* aHC = HCE->GetHC(fSensitiveDetector_ID);
            fSensitiveDetectorHC = (SensitiveDetectorHitsCollection*)(aHC);
        }
    }

    //instantiating The Analysis Manager
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
	G4double eHit = 0.;

    //filling the scoring ntuple
    if (fSensitiveDetectorHC) {
        int vNumberOfHit = fSensitiveDetectorHC->entries();
        for (int i=0; i<vNumberOfHit; i++) {
            SensitiveDetectorHit* aHit = (*fSensitiveDetectorHC)[i];
			eHit = aHit->GetEnergy();
			Edep += eHit;
        }
		if (Edep > 0.) {
		   	analysisManager->FillNtupleDColumn(0,0,Edep/CLHEP::keV);
			analysisManager->AddNtupleRow(0);   
		}
    }
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

