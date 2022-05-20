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

#include "DetectorConstruction.hh"
#include "SensitiveDetector.hh"

#include "G4RunManager.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"

#include "G4ExtrudedSolid.hh"
#include "G4Polyhedra.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LogicalVolumeStore.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "G4VSensitiveDetector.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4PSDoseDeposit.hh"
#include "G4PSFlatSurfaceFlux.hh"
#include "G4PSEnergyDeposit.hh"
#include "G4SDParticleFilter.hh"
#include "G4VSDFilter.hh"
#include "G4SDKineticEnergyFilter.hh"

#include "G4FieldManager.hh"
#include "G4UniformMagField.hh"
#include "G4TransportationManager.hh"
#include "CLHEP/Units/SystemOfUnits.h"

DetectorConstruction::DetectorConstruction()
{
	checkOverlaps = false;
	G4cout << "### DetectorConstruction instantiated ###" << G4endl;
}

DetectorConstruction::~DetectorConstruction() {}

G4VPhysicalVolume *DetectorConstruction::Construct()
{

	// Build the world !
	G4Material *G4_Air = G4NistManager::Instance()->FindOrBuildMaterial("G4_AIR");
	G4double worldSize = 4. * m;
	G4Box *worldSolid = new G4Box("World",
								  worldSize / 2.,
								  worldSize / 2.,
								  worldSize / 2.);
	lv_world = new G4LogicalVolume(worldSolid,
								   G4_Air,
								   "World");
	lv_world->SetVisAttributes(G4VisAttributes::GetInvisible());
	G4PVPlacement *pv_world = new G4PVPlacement(0,
												G4ThreeVector(),
												lv_world,
												"World",
												0,
												false,
												0);
	placements["world"] = pv_world;

	// prepare the beam pipe, a beryllium GTubs with bvacuum inside
	double pipe_dz = worldSize / 2 - 0.1 * mm;
	// beam pipe
	G4Material *G4_Be = G4NistManager::Instance()->FindOrBuildMaterial("G4_Be");
	G4VSolid *sd_bpipe = new G4Tubs("beampipe", 0., 5 * cm, pipe_dz, 0., twopi * rad);
	G4LogicalVolume *lv_bpipe = new G4LogicalVolume(sd_bpipe, G4_Be, "lb_bpipe", 0, 0, 0);
	G4VisAttributes *visatt = new G4VisAttributes(G4Colour(1, 0, 0));
	lv_bpipe->SetVisAttributes(visatt);
	G4VPhysicalVolume *pv_bpipe = new G4PVPlacement(0, G4ThreeVector(0., 0., 0), "pv_bpipe", lv_bpipe, pv_world, false, 0);
	placements["beam_pipe"] = pv_bpipe;

	// Vaccuum in the beam pipe
	G4Material *G4_vacuum = G4NistManager::Instance()->FindOrBuildMaterial("G4_Galactic");
	G4VSolid *sd_vacuum = new G4Tubs("vacuum", 0., 4.5 * cm, pipe_dz, 0., twopi * rad);
	G4LogicalVolume *lv_vacuum = new G4LogicalVolume(sd_vacuum, G4_vacuum, "lv_vacuum", 0, 0, 0);
	G4VisAttributes *visatttransp = new G4VisAttributes(G4Colour(0, 0, 1));
	visatttransp->SetVisibility(visatttransp);
	lv_vacuum->SetVisAttributes(visatttransp);
	G4VPhysicalVolume *pv_vacuum = new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), "pv_vacuum", lv_vacuum, pv_bpipe, false, 0);
	placements["vacuum"] = pv_vacuum;

	// Now prepare the sensors
	G4Material *G4_Si = G4NistManager::Instance()->FindOrBuildMaterial("G4_Si");

	G4Box *sd_sensor_box = new G4Box("sensor",
									 40 * cm,
									 40 * cm,
									 2 * mm);
	G4Tubs *si_sensor_hole = new G4Tubs("hole",
										0 * cm,
										6 * cm,
										10 * mm, 0., twopi);
	G4SubtractionSolid *sd_sensor = new G4SubtractionSolid("sensor", sd_sensor_box, si_sensor_hole);
	lv_sensor = new G4LogicalVolume(sd_sensor,
									G4_Si,
									"lv_sensor");

	// First group of sensors
	for (int i = 0; i < 10; i++)
	{
		double z = -100.0 + 5.0 * i;
		G4VPhysicalVolume *p = new G4PVPlacement(nullptr,
												 G4ThreeVector(0., 0., z * cm),
												 lv_sensor,
												 "pv_sensor",
												 lv_world,
												 false,
												 i);
		sensor_placements["A" + std::to_string(i)] = p;
	}

	// Second group of sensors
	for (int i = 0; i < 10; i++)
	{
		double z = 100.0 + 5.0 * i;
		G4VPhysicalVolume *p = new G4PVPlacement(nullptr,
												 G4ThreeVector(0., 0., z * cm),
												 lv_sensor,
												 "pv_sensor",
												 lv_world,
												 false,
												 100 + i);
		sensor_placements["B" + std::to_string(i)] = p;
	}

	G4Box *sd_hsensor = new G4Box("hsensor",
								  80 * cm,
								  80 * cm,
								  2 * mm);

	lv_hsensor = new G4LogicalVolume(sd_hsensor,
									 G4_Si,
									 "lv_hsensor");
	G4RotationMatrix *myRotationMatrix = new G4RotationMatrix();
	myRotationMatrix->rotateX(-90. * CLHEP::deg);
	G4VPhysicalVolume *h0 = new G4PVPlacement(myRotationMatrix,
											  G4ThreeVector(0., 80. * cm, -100 * cm),
											  lv_hsensor,
											  "pv_hsensor",
											  lv_world,
											  false,
											  0);
	sensor_placements["H0"] = h0;

	G4VPhysicalVolume *h1 = new G4PVPlacement(myRotationMatrix,
											  G4ThreeVector(0., -80. * cm, -100 * cm),
											  lv_hsensor,
											  "pv_hsensor",
											  lv_world,
											  false,
											  1);
	sensor_placements["H1"] = h1;
	return pv_world;
}

void DetectorConstruction::ConstructSDandField()
{

	// Setting the field
	G4FieldManager *globalFieldMgr =
		G4TransportationManager::GetTransportationManager()->GetFieldManager();
	G4MagneticField *magField = new G4UniformMagField(G4ThreeVector(0.5 * tesla, 0., 0.));
	globalFieldMgr->SetDetectorField(magField);
	globalFieldMgr->CreateChordFinder(magField);

	// Sensitive Detector
	G4VSensitiveDetector *vDetector = new SensitiveDetector("det");
	G4SDManager::GetSDMpointer()->AddNewDetector(vDetector);
	lv_sensor->SetSensitiveDetector(vDetector);
	lv_hsensor->SetSensitiveDetector(vDetector);
	
	// Dose Scoring
	G4MultiFunctionalDetector *multisd = new G4MultiFunctionalDetector("multisd");
	G4VPrimitiveScorer *dosedet = new G4PSDoseDeposit("dose");
	multisd->RegisterPrimitive(dosedet);
	G4SDManager::GetSDMpointer()->AddNewDetector(multisd);
	SetSensitiveDetector("lv_sensor", multisd);

	// if I use the following expression
	// myPolyhedraLogical->SetSensitiveDetector(multisd);
	// the previous SD does not work anymore (I can associate only one SD to a logical volume)

	/*
		G4SDParticleFilter* particleFilter = new G4SDParticleFilter("particleFilter");
		particleFilter->add("e-");
		dosedet->SetFilter(particleFilter);
	*/

	/*
		G4SDKineticEnergyFilter* energyFilter = new G4SDKineticEnergyFilter("particleFilter");
		energyFilter->SetKineticEnergy(10.*keV,100.*keV);
		dosedet->SetFilter(energyFilter);
	 */
}
