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


#include "PrimaryGeneratorAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4AutoLock.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"

#include "Randomize.hh"

#include "G4GeneralParticleSource.hh"
#include "G4ParticleGun.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"



PrimaryGeneratorAction::PrimaryGeneratorAction() : fGPS{new G4GeneralParticleSource()}, fGun{new G4ParticleGun()}
{
	G4cout << "### PrimaryGeneratorAction instantiated ###" << G4endl;

	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle("pi+");
	fGun->SetParticleDefinition(particle);
	fGun->SetParticleMomentumDirection(G4ThreeVector(0., 0., 1.));
	fGun->SetParticleEnergy(1 * GeV);
	fGun->SetParticlePosition(G4ThreeVector(0., 0., -75 * cm));
}



PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
	delete fGPS;
	delete fGun;
}



void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	if (IWantGPS)
	{
		fGPS->GeneratePrimaryVertex(anEvent);
	}
	else
	{
		G4double cosTheta = 1 - 2 * G4UniformRand();// / 100;
		G4double sinTheta2 = 1. - cosTheta * cosTheta;
		if (sinTheta2 < 0.)
			sinTheta2 = 0.;
		G4double sinTheta = std::sqrt(sinTheta2);
		G4double phi = CLHEP::twopi * G4UniformRand();
		fGun->SetParticleMomentumDirection(G4ThreeVector(sinTheta * std::cos(phi),
														 sinTheta * std::sin(phi), cosTheta));
		fGun->GeneratePrimaryVertex(anEvent);
	}
}


