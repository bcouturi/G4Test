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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "RunAction.hh"
#include "RunActionMessenger.hh"
#include "Analysis.hh"
#include "Run.hh"

#include "G4SDManager.hh"

#include "G4AccumulableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction():
G4UserRunAction(),
fIsFileOpened(false)
{
	G4RunManager::GetRunManager()->SetPrintProgress(1000);
    
    fMessenger = new RunActionMessenger(this);

    fFileName = "myOutput";
    
    G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    analysisManager->SetDefaultFileType("root"); //default is root 
    analysisManager->SetFileName(fFileName);
    analysisManager->SetNtupleMerging(true);
    analysisManager->SetVerboseLevel(0);
    
    //Creating the scoring ntuple
    analysisManager->CreateNtuple("ntuple","scoring");
    analysisManager->CreateNtupleDColumn("Edep");
    analysisManager->FinishNtuple();
      	
  	G4cout << "### RunAction instantiated ###" << G4endl;    	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction() {
   	//write results
	if (fIsFileOpened) {
    	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance(); 
      	analysisManager->Write();
      	//analysisManager->CloseFile();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun() {return new Run;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
	//merging is only for root
    if (fFileName.find(".csv") != std::string::npos) {
    	G4AnalysisManager::Instance()->SetNtupleMerging(false);
    }
    
	//open the output file         
  	if (!fIsFileOpened) {
    	G4AnalysisManager::Instance()->OpenFile(fFileName);
      	fIsFileOpened = true;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{  	   
/*    
    G4int nofEvents = run->GetNumberOfEvent();
  	if (nofEvents == 0) return;
  
  	//results from the Scorer
  	const Run* aRun = static_cast<const Run*>(run);
  	G4int nbGoodEvents = aRun->GetNbGoodEvents(); 
  	G4double sumDose   = aRun->GetSumDose();                           
        
  	//print
  	if (IsMaster()) {
    	G4cout
     	<< G4endl
     	<< "--------------------End of Global Run-----------------------"
     	<< G4endl
     	<< " The run had " << nofEvents << " events";
  	} else {
    	G4cout
     	<< G4endl
     	<< "--------------------End of Local Run------------------------"
     	<< G4endl
     	<< " The run had " << nofEvents << " events";
  	}      
  	
  	G4cout
    << "; Nb of 'good' events: " << nbGoodEvents << G4endl 
    << G4endl
    << " Total dose in solid: " << sumDose/gray << " Gray" << G4endl  
    << "------------------------------------------------------------" << G4endl 
    << G4endl;


	//Write Results from Scorer on a text file  
    if (IsMaster()) {
    	std::ofstream fFileOut;
        fFileOut.open("SimulationSummary.dat", std::ofstream::out | std::ofstream::app);
        fFileOut << "seed: " << G4Random::getTheSeed() 
        		 << " | " << sumDose/gray << " Gray"    
        	 	 << std::endl;
        fFileOut.close();
  	}
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetFileName(G4String filename) 
{
	if (filename != "") fFileName = filename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

