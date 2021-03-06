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
/// \file analysis/AnaEx01/AnaEx01.cc
/// \brief Main program of the analysis/AnaEx01 example
//
//
// $Id: AnaEx01.cc 73919 2013-09-17 07:38:47Z gcosmo $
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

#include "Randomize.hh"

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
    // Choose the Random engine
    G4Random::setTheEngine(new CLHEP::RanecuEngine);
	 // G4long is at least 32 bits. 
    G4long seed = time(NULL);
    G4Random::setTheSeed(seed);

    // Construct the default run manager
#ifdef G4MULTITHREADED
	 G4int nThreads = 2;
	 if(argc == 3) {
		 nThreads = strtol(argv[2], nullptr, 10);
	 }
	 G4cout<<"RUNNING MULTITHREADED WITH "<<nThreads<<" THREADS"<<G4endl;
	 G4MTRunManager* runManager = new G4MTRunManager;
	 runManager->SetNumberOfThreads(nThreads);
#else
	 G4cout<<"NOT RUNNING MULTITHREADED"<<G4endl;
	 G4RunManager* runManager = new G4RunManager;
#endif

	 // Set mandatory initialization classes
	 DetectorConstruction* detector = new DetectorConstruction;
	 runManager->SetUserInitialization(detector);
	 runManager->SetUserInitialization(new PhysicsList);
	 runManager->SetUserInitialization(new ActionInitialization(detector));

	 // We don't initialize the G4 kernel at run time so the physics list can be changed!

	 // Get the pointer to the User Interface manager
	 G4UImanager* UImanager = G4UImanager::GetUIpointer();

#ifdef G4VIS_USE
	 G4VisManager* visManager = new G4VisExecutive;
	 visManager->Initialize();
#endif

	 if(argc != 1) { // batch mode
		 G4String command = "/control/execute ";
		 G4String fileName = argv[1];
		 UImanager->ApplyCommand(command+fileName);
	 } else { // interactive mode : define visualization and UI terminal
#ifdef G4UI_USE
		 G4UIExecutive* ui = new G4UIExecutive(argc, argv);
#endif
		 UImanager->ApplyCommand("/control/execute run.mac");

#ifdef G4UI_USE
		 ui->SessionStart();

		 delete ui;
#endif
	 }

#ifdef G4VIS_USE
		 delete visManager;
#endif

	 // Job termination
	 delete runManager;

	 return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
