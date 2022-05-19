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
// $Id: Run.cc 68058 2013-03-13 14:47:43Z gcosmo $
// edited by gpaterno (june 2016)
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Run.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::Run():
G4Run(), 
fCollID_dose(-1),
fGoodEvents(0),
fSumDose(0.)
{
	G4cout << "### Run instantiated ###" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Run::~Run(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::RecordEvent(const G4Event* event)
{
    //Dose hits collection ID 
	if (fCollID_dose < 0) {
   		fCollID_dose = G4SDManager::GetSDMpointer()->GetCollectionID("multisd/dose");
   		G4cout << "fCollID: " << fCollID_dose << G4endl;   
   	}
                
  	//get all Hits Collections of this Event
  	G4HCofThisEvent* HCE = event->GetHCofThisEvent();
  	if (!HCE) return;
    
	//get the hits Collection, which is basically a map: (hit,dose*), with ID = fCollID_dose
  	G4THitsMap<G4double>* evtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID_dose));
    
	//map iterator
  	std::map<G4int,G4double*>::iterator itr;
  
  	//Dose deposit 
  	G4double dose = 0.;
               
  	for (itr = evtMap->GetMap()->begin(); itr != evtMap->GetMap()->end(); itr++) {
    	dose = *(itr->second);
  		fSumDose += dose;
  		if (dose > 0.) fGoodEvents ++;   //gpaterno
   	}

	//call base method
  	G4Run::RecordEvent(event);      
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Run::Merge(const G4Run* aRun)
{
	//accumulate variables from each instance (one for thread) of Run
	const Run* localRun = static_cast<const Run*>(aRun);
  	fGoodEvents += localRun->fGoodEvents;
  	fSumDose    += localRun->fSumDose;

	//call base method (to pass accumulated variables to the master thread)
  	G4Run::Merge(aRun); 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

