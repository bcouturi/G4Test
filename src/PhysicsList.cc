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

#include "PhysicsList.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4ProcessManager.hh"
#include "G4Threading.hh"

//Particles
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

//EM Physics Lists
#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmLowEPPhysics.hh"
#include "G4EmStandardPhysics_option4.hh"

//EM options
#include "G4EmSaturation.hh"
#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"

//Hadronic and Extra Physics Lists
#include "G4EmExtraPhysics.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
#include "G4HadronElasticPhysicsHP.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"

//Optical processes
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"

//Decays
#include "G4Decay.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4PhysicsListHelper.hh"
#include "G4RadioactiveDecay.hh"
#include "G4NuclideTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::PhysicsList(): G4VUserPhysicsList()
{
	G4cout << "### PhysicsList instantiated ###" << G4endl; 

	G4LossTableManager::Instance();
	
	//Hadronic and Extra Physics option
	IWantHadExtraPhysics = false;
	
	//set default cuts value
	defaultCutValue = 0.001*mm;
	  	    
    //set verbosity
  	SetVerboseLevel(1);
  	
  	//particle list
  	fParticleList = new G4DecayPhysics(verboseLevel);

  	//EM Physics
  	fEmPhysicsList = new G4EmStandardPhysics_option4(verboseLevel);  
  	
  	//Hadronic and Extra Physics	
	fEmExtraPhysicsList = new G4EmExtraPhysics(verboseLevel);         
    fHadPhysicsList = new G4HadronPhysicsQGSP_BIC_HP(verboseLevel);
    fHadElPhysicsList = new G4HadronElasticPhysicsHP(verboseLevel);  
    fIonPhysicsList = new G4IonPhysics(verboseLevel);        
	fStoppingPhysicsList = new G4StoppingPhysics(verboseLevel);  
	
/*  		
	//add new units for radioactive decays
  	const G4double minute = 60*second;
  	const G4double hour   = 60*minute;
  	const G4double day    = 24*hour;
  	const G4double year   = 365*day;
  	new G4UnitDefinition("minute", "min", "Time", minute);
  	new G4UnitDefinition("hour",   "h",   "Time", hour);
  	new G4UnitDefinition("day",    "d",   "Time", day);
  	new G4UnitDefinition("year",   "y",   "Time", year);
*/

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PhysicsList::~PhysicsList() 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructParticle()
{
	fParticleList->ConstructParticle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructProcess()
{	
	//transportation
 	AddTransportation();
 	
 	//EM physics
 	fEmPhysicsList->ConstructProcess();
 	
 	//Hadronic and Extra Physics
 	if (IWantHadExtraPhysics) {
 		fEmExtraPhysicsList->ConstructProcess();
    	fHadPhysicsList->ConstructProcess();
		fHadElPhysicsList->ConstructProcess();
    	fIonPhysicsList->ConstructProcess();
    	fStoppingPhysicsList->ConstructProcess();
    }

	//Atomic deexcitation
	G4VAtomDeexcitation* de = new G4UAtomicDeexcitation();
	de->SetFluo(true);	//To activate deexcitation processes and fluorescence
	de->SetAuger(true);	//To activate Auger effect if deexcitation is activated
	de->SetAugerCascade(true);	//To activate Auger Cascade if deexcitation is activated
	de->SetPIXE(true);  //To activate Particle Induced X-Ray Emission (PIXE)
	G4LossTableManager::Instance()->SetAtomDeexcitation(de);
 	
  	//decay
  	ConstructDecay(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PhysicsList::ConstructDecay()
{
  	//Add Decay Process
	G4Decay* theDecayProcess = new G4Decay();
  	GetParticleIterator()->reset();
  	while((*GetParticleIterator())()){
    	G4ParticleDefinition* particle = GetParticleIterator()->value();
    	G4ProcessManager* pmanager = particle->GetProcessManager();
    	if (theDecayProcess->IsApplicable(*particle)) {
      		pmanager->AddProcess(theDecayProcess);
      		//set ordering for PostStepDoIt and AtRestDoIt
      		pmanager->SetProcessOrdering(theDecayProcess, idxPostStep);
      		pmanager->SetProcessOrdering(theDecayProcess, idxAtRest);
    	}
  	}

	//Radioactive decay processes
	G4RadioactiveDecay* radioactiveDecay = new G4RadioactiveDecay();
  	radioactiveDecay->SetARM(true);	//Atomic Rearangement
  	G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();  
  	ph->RegisterProcess(radioactiveDecay, G4GenericIon::GenericIon());
  	G4NuclideTable::GetInstance()->SetThresholdOfHalfLife(0.1*picosecond);
  	G4NuclideTable::GetInstance()->SetLevelTolerance(1.0*eV);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

