


#ifndef SHPRIMARYGENERATORACTION_H
#define SHPRIMARYGENERATORACTION_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4ThreeVector.hh"
#include "G4ios.hh"
#include "G4Event.hh"
#include "G4VPrimaryGenerator.hh"

#include "TRandom.h"

#include "SHDefs.hh"
#include "HGCReadoutModule.hh"

class SHPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
	public:
		//	Constructor and Destructor
		SHPrimaryGeneratorAction(RunParams params, TTree *tree);
		SHPrimaryGeneratorAction(RunParams, HGCReadoutModule*);
		virtual ~SHPrimaryGeneratorAction();

		//	Standard func	
		void GeneratePrimaries(G4Event *anEvent);

		//	To keep things simplier
		G4ParticleGun *particleGun;
		G4String primName;
		G4double primEnergy;
		G4ThreeVector primPos;
		G4ThreeVector primDir;
		G4int verbosityLevel;

		RunParams runParams;
		TTree *shTree;
		HGCReadoutModule *_readout;
		TRandom _rand;
};

#endif
