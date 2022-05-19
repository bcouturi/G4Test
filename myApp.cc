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


#include "G4MTRunManager.hh"
#include "G4RunManager.hh"

#include "G4ScoringManager.hh"
#include "G4UImanager.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "UserActionInitialization.hh"
#include "Analysis.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "QGSP_BERT.hh"
#include "QGSP_BIC_HP.hh"
#include "G4EmStandardPhysics_option3.hh"
#include "G4VUserPhysicsList.hh"
#include "G4PhysListFactory.hh"
#include "G4OpticalPhysics.hh"

#include "Randomize.hh"
#include <sys/time.h>
#include "G4Timer.hh"



int main(int argc,char** argv)
{
	//Get current time
    G4Timer* theTimer = new G4Timer();
    theTimer->Start();
    
    //Choose the Random engine
	G4Random::setTheEngine(new CLHEP::RanecuEngine);
	time_t systime = time(NULL);
    G4int seed = (long)(systime*G4UniformRand());
    G4Random::setTheSeed(seed);

    //Construct the default run manager
#undef G4MULTITHREADED
#ifdef G4MULTITHREADED
    G4MTRunManager* runManager = new G4MTRunManager;
    int vNumberOfThreads = 1;
    if (argc>2) {
        vNumberOfThreads = atoi(argv[2]);
    }
    if (vNumberOfThreads > 0) {
        runManager->SetNumberOfThreads(vNumberOfThreads);
    }
    G4cout << "### MT MODE ON " << runManager->GetNumberOfThreads() << " ###" << G4endl;
#else
    G4RunManager* runManager = new G4RunManager;
    G4cout << "### MT MODE OFF ###" << G4endl;
#endif
    
    //Activate UI-command base scorer
    // G4ScoringManager* scManager = G4ScoringManager::GetScoringManager();
    // scManager->SetVerboseLevel(0);
    
    //Set mandatory initialization classes    
    runManager->SetUserInitialization(new DetectorConstruction);

	G4bool IWantModularPL = true;
	if (IWantModularPL) {
		G4PhysListFactory* physListFactory = new G4PhysListFactory();
		G4VModularPhysicsList* physics = physListFactory->GetReferencePhysList("QGSP_BIC_HP");
		physics->SetVerboseLevel(1);
		physics->RegisterPhysics(new G4OpticalPhysics());
		physics->ReplacePhysics(new G4EmStandardPhysics_option3());
		runManager->SetUserInitialization(physics);
	} else {
		runManager->SetUserInitialization(new PhysicsList);
	}

    runManager->SetUserInitialization(new UserActionInitialization());
 
    //Get the pointer to the User Interface manager
    G4UImanager* UI = G4UImanager::GetUIpointer();  
    
    if (argc!=1) {
        //Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UI->ApplyCommand(command+fileName);
    }
    
    else {
        //Define visualization and UI terminal for interactive mode
        G4VisManager* visManager = new G4VisExecutive;
        visManager->Initialize();

        G4UIExecutive* ui = new G4UIExecutive(argc,argv);
		UI->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
	        
		delete visManager;
    }
    
    //Job termination
    delete runManager;
    
    G4AnalysisManager* man = G4AnalysisManager::Instance();
    man->CloseFile();
    
    theTimer->Stop();
    G4cout << "Execution terminated" << G4endl;
    G4cout << (*theTimer) << G4endl;
    delete theTimer;
    
    return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

