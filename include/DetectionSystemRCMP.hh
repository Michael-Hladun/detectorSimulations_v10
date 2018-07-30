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
// $Id: DetectionSystemRCMP.hh,v 1.0 2018-06-12 $
// GEANT4 tag $Name: geant4-10-04-p02 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DETECTIONSYSTEMRCUBE_HH
#define DETECTIONSYSTEMRCUBE_HH

#include "globals.hh"

class DetectionSystemRCMP
{
public:
    DetectionSystemRCMP();
    ~DetectionSystemRCMP();

    // Assembly volumes
    G4AssemblyVolume* fAssembly_Si;
    G4AssemblyVolume* fAssembly_PCB;
    
private:
    // Logical volumes
    G4LogicalVolume* fRCMPSilicon_log;
    G4LogicalVolume* fRCMPPCB_log;
    G4LogicalVolume* fRCMPPCB_withcut_log;
public:
    // Solid parameters
    G4double fPCB_thickness;
    G4double fPCB_length;
    G4double fPCB_width;

    G4double fSi_thickness;
    G4double fSi_length;

    G4Colour fPCB_colour;
    G4Colour fSi_colour;

// Called from DetectorConstruction.cc
public:
    G4int Build();
    G4int PlaceDetector(G4LogicalVolume* expHallLog, G4double beam_opening);

// Called from within DetectionSystemRCMP.cc, hence these functions are private, and not public
private:
    // Construction methods
    G4int SetUpFace();
};

#endif
