#include "SHPrimaryGeneratorAction.hh"

#include <iostream>

#include "G4HEPEvtInterface.hh"
#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "CLHEP/Random/RandFlat.h"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

using namespace CLHEP;
using namespace std;

//	0-7 are old energies(for EM)
//	8-12 are new energies(for HAD)
//
const int numEnergies = 13;
const int numPrims = 3;
G4double primEnergies[numEnergies] = {
	1.*GeV, 2.*GeV, 4.*GeV, 8.*GeV, 16.*GeV, 32.*GeV, 50.*GeV, 60.*GeV,
	20*GeV, 50*GeV, 100*GeV, 200*GeV, 300*GeV
};
G4String primNames[numPrims] = {"e-", "pi-", "mu-"};



//	Constructor
//
SHPrimaryGeneratorAction::SHPrimaryGeneratorAction(RunParams params, 
		HGCReadoutModule *readout)
	: runParams(params)
{
	_readout = readout;
	_rand.SetSeed(runParams.seed);

	//	Define a generator
	//	Random dir(up to dims of Endcap)
	//	Random Energy(up to 200GeV)
	//
	primName = primNames[runParams.iPrim];
//	G4PrimarypParticle *primParticle = new G4PrimaryParitcle()
//	primName = "d_quark";
	primEnergy = _rand.Uniform(1, 200)*GeV;

	particleGun = new G4ParticleGun(1);
	particleGun->SetParticleEnergy(primEnergy);
	G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
	G4ParticleDefinition *particle = particleTable->FindParticle(primName);
//	G4PrimaryParticle *primParticle = new G4PrimaryParitcle(particle);
	

	particleGun->SetParticleDefinition(particle);

	primPos = G4ThreeVector(0*mm, 0*mm, 0.*mm);
//	primDir = G4ThreeVector(1.*m, 0., 0);
	
	G4double dirx = _rand.Uniform(-1600, 1600)*mm;
	G4double diry = _rand.Uniform(-1600, 1600)*mm;
	G4double dirz = 3150*mm;
	primDir = G4ThreeVector(dirx, diry, dirz);
//	primDir = G4ThreeVector(0, 0, -3150*mm);
	particleGun->SetParticlePosition(primPos);
	particleGun->SetParticleMomentumDirection(primDir);
	particleGun->SetParticleEnergy(primEnergy);
}

//	Destructor
//
SHPrimaryGeneratorAction::~SHPrimaryGeneratorAction()
{
	delete particleGun;
}

//	Generate Primary Event
//
void SHPrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
	//	All gun's settings have been set up in Constructor
	//
	//particleGun->SetParticleEnergy(10.*keV);
	G4cout << particleGun->GetParticleEnergy()/GeV << "  " 
		<< particleGun->GetParticleDefinition()->GetParticleName()
		<< G4endl; 

	primEnergy = _rand.Uniform(1, 200)*GeV;
	particleGun->SetParticleEnergy(primEnergy);

	G4double dirx = _rand.Uniform(-1600., 1600.)*mm;
	G4double diry = _rand.Uniform(-1600., 1600.)*mm;
	G4double dirz = 3150*mm;
	primDir = G4ThreeVector(dirx, diry, dirz);
//	primDir = G4ThreeVector(0, 0, -3150*mm);
	particleGun->SetParticleMomentumDirection(primDir);

/*
	G4PrimaryParticle *primParticle = new G4PrimaryParticle(1);
	G4PrimaryVertex *primVertex = new G4PrimaryVertex(G4ThreeVector(0, 500*mm, 
				0), 0);
	primParticle->SetMomentum(0, 0, 1*GeV);
	primVertex->SetPrimary(primParticle);
	anEvent->AddPrimaryVertex(primVertex);
*/

	particleGun->GeneratePrimaryVertex(anEvent);
}
