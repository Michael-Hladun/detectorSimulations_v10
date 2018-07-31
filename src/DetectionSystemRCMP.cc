#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"

#include "G4Material.hh"

#include "globals.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4AssemblyVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "DetectionSystemRCMP.hh"

#include "G4SystemOfUnits.hh" // new version geant4.10 requires units


DetectionSystemRCMP::DetectionSystemRCMP() :
	// Logical Volumes
	fRCMPSilicon_log(0),
	fRCMPPCB_log(0),
	fRCMPPCB_withcut_log(0)
{
	// Solid parameters - Based off Micron Semiconductor Limited A-3893
	fPCB_thickness = (5. / 2.)*mm;
	fPCB_length = (67.975 / 2.)*mm;
	fPCB_width = (67.975 / 2.)*mm;

	fSi_thickness = (1.5 / 2.)*mm;
	fSi_length = (63.96 / 2.)*mm;

	fPCB_colour = G4Colour(0.5, 0.0, 0.5);
	fSi_colour = G4Colour(1.0, 1.0, 1.0);
}

DetectionSystemRCMP::~DetectionSystemRCMP() {
	// Logical volumes
	delete fRCMPSilicon_log;
	delete fRCMPPCB_log;
	delete fRCMPPCB_withcut_log;
}

G4int DetectionSystemRCMP::Build() {
	// Build assembly volumes
	fAssembly_PCB = new G4AssemblyVolume();
	fAssembly_Si = new G4AssemblyVolume();
	// Add a face to the cube.
	SetUpFace();	
	
	return 1;
}//end ::Build



G4int DetectionSystemRCMP::PlaceDetector(G4LogicalVolume* expHallLog, G4double beam_opening) {
	// beam_opening is the width the beam enters:
	if(beam_opening > 60*mm){
		std::cout<<std::endl<<"Beam opening specified is too large. DSSD's hitting vacuum chamber."<<std::endl<<std::endl;
		return 0;
	}
	if(beam_opening < 0.*mm){
		std::cout<<std::endl<<"Beam opening must be >=0."<<std::endl<<std::endl;
		return 0;
	}

	G4double sin45 = std::sin(45.*3.14159265/180.);

	// specify front DSSD offsets from the middle bisection:
	G4double up_offset = 2.0*mm;
	up_offset = (up_offset * sin45 / std::sqrt(2.))*mm;

	// specify back DSSD offsets from the middle bisection:
	G4double down_offset = 2.0*mm;
	down_offset = -(down_offset * sin45 / std::sqrt(2.))*mm;
	
	// Distance of DSSD centers from origin:
	G4double r_vec = fPCB_length + fPCB_thickness;
 
	// Pivoting DSSD's on the bisection axis (with 2mm offset):
	// Since the trigonometry of this problem is not trivial, I loop through angles of DSSD's,
	// making translations after each rotation, in order to keep the pivot points stationary
	// while being able to specify different beam openings above.
	G4double Xmove_up, Xmove_over, Xbeam_opening = 0.0*mm;
	G4double Xtheta = 0.0;
	if(beam_opening > 0.0*mm){
		do{
			Xtheta += 0.001*deg;
			Xmove_up = std::sin(Xtheta*3.14159265/180.) * fPCB_length;
			Xbeam_opening = 8.0 * std::sin((Xtheta/2.0)*3.14159265/180.) * fPCB_length;
		}while(Xbeam_opening < beam_opening);

		// adjustment needed to make Xtheta = 45deg at std::sqrt(2.0 * 63.975 * 63.975), i.e. the hypotenuse of the active Si.
		Xtheta = Xtheta * 45.0 / 38.8714;
		
		Xmove_over = (fPCB_length * (1.0 - std::cos(Xtheta*3.14159265/180.)))*mm;
		// std::cout<<std::endl<<"side move: "<<Xmove_over<<std::endl<<std::endl;
	}
	std::cout<<std::endl<<"theta: "<<Xtheta<<std::endl<<std::endl;
	// std::cout<<std::endl<<"beam opening diff: "<<Xbeam_opening-beam_opening<<std::endl<<std::endl;
	std::cout<<std::endl<<"beam opening: "<<Xbeam_opening<<std::endl<<std::endl;
		

	// Moving the DSSD's in any direction:
	//	{x, y, z}
	G4double trans[6][3] = {
		{up_offset - Xmove_over, 			0.*mm, 		r_vec + up_offset + Xmove_up},// 0 - upstream
		{r_vec + up_offset + Xmove_up,	 	0.*mm, 		up_offset - Xmove_over},// 1 - upstream
		{0.*mm, 							r_vec, 		0.*mm},// 2 - top
		{0.*mm, 							-r_vec, 	0.*mm},// 3 - bottom
		{down_offset, 						0.*mm, 		-r_vec + down_offset},// 4 - downstream
		{-r_vec + down_offset,				0.*mm, 		down_offset}// 5 - downstream
	};


	// Rotating the cube faces:
		//{x, y, z}
	G4double rot[6][3] = {
		{0.*deg, 	(-Xtheta)*deg, 			90.*deg},
		{90.0*deg, 	(+Xtheta+270.0)*deg, 	90.0*deg},
		{0.*deg, 	270.*deg, 				0.*deg},
		{0.*deg, 	270.*deg, 				0.*deg},
		{0.*deg, 	180.*deg, 				90.*deg},
		{90.0*deg, 	90.0*deg, 				90.0*deg}
	};

	// Variables which change for each face:
	G4ThreeVector movePCB;
	G4ThreeVector moveSi;
	G4RotationMatrix* rotate;// Generic rotation object for cube faces.

	// MakeImprint variables. Y rotation to line up beam gap with incoming beam.
	G4ThreeVector Ta;
	G4RotationMatrix* Ra = new G4RotationMatrix;
	Ra->rotateY(315.0*deg);

	// Loop over 6 faces of the RCMP cube:
		// i==0,1 - Front and Back (z displacement) faces
		// i==2,3 - Top and Bottom (y-displacement) faces
		// i==4,5 - Left and Right (x displacement) faces
	for(G4int i=0; i<6; i++){

		movePCB = G4ThreeVector(trans[i][0], trans[i][1], trans[i][2]);
		moveSi =  G4ThreeVector(trans[i][0], trans[i][1], trans[i][2]);

		rotate = new G4RotationMatrix(rot[i][0], rot[i][1], rot[i][2]);

		fAssembly_Si->AddPlacedVolume(fRCMPSilicon_log, moveSi, rotate);
		fAssembly_PCB->AddPlacedVolume(fRCMPPCB_log, movePCB, rotate);

		// Put blank translations and rotations in MakeImprint:
		fAssembly_Si->MakeImprint(expHallLog, Ta, Ra, i);
		fAssembly_PCB->MakeImprint(expHallLog, Ta, Ra, i);
	}

	return 1;
}// end ::PlaceDetector


// Creates the foundation for a face of the RCMP cube.
G4int DetectionSystemRCMP::SetUpFace() {

	// Materials
	// PCB - FR-4 Substrate material (similar to quartz)
	// Z, A, and rho taken from: http://personalpages.to.infn.it/~tosello/EngMeet/ITSmat/SDD/SDD_G10FR4.html
	G4Material* mat_PCB = new G4Material("FR4_substrate", 9.4328, 18.9415*g/mole, 1.80*g/cm3);
	G4Material* mat_Si = G4Material::GetMaterial("Silicon");

	if(!mat_PCB) {
		G4cout << " ----> Material  FR-4 PCB Substrate not made, cannot build the detector! " << G4endl;
		return 0;
	}//else{std::cout<<"good\n";}
	if(!mat_Si) {
		G4cout << " ----> Material  Silicon not found, cannot build the detector! " << G4endl;
		return 0;
	}//else{std::cout<<"good\n";}
	

	// Solid Volumes
		// Solid state silicon detector:
	G4Box* si_face = new G4Box("siliconFace", fSi_length, fSi_length, fSi_thickness);

		// Create a box for cutting the PCB:
	G4Box* si_face_for_cut = new G4Box("siliconFace", fSi_length, fSi_length, fPCB_thickness+0.5*mm);

		// Solid volume which will have it's center cut out for the silicon face:
	G4Box* PCB_Face = new G4Box("PCB", fPCB_length, fPCB_width, fPCB_thickness);

		// PCB with cut to produce 2mm border around active Si:
	G4SubtractionSolid* PCB_withcut = new G4SubtractionSolid("PCB_WithCut", PCB_Face, si_face_for_cut);


	// Vis Attributes
	G4VisAttributes* si_face_vis_att = new G4VisAttributes(fSi_colour);
	si_face_vis_att->SetVisibility(true);
	si_face_vis_att->SetForceSolid(true);
	

	G4VisAttributes* PCB_face_vis_att = new G4VisAttributes(fPCB_colour);
	PCB_face_vis_att->SetVisibility(true);
	PCB_face_vis_att->SetForceSolid(true);


	// Logical Volumes
		// PCB with cut added as logical volume:
	if(fRCMPPCB_log == NULL) {
		fRCMPPCB_log = new G4LogicalVolume(PCB_withcut, mat_PCB, "fRCMPPCB_log");
		fRCMPPCB_log->SetVisAttributes(PCB_face_vis_att);
		//std::cout<<"HELLO fRCMPPCB_log\n";
	}
	if(fRCMPSilicon_log == NULL) {
		fRCMPSilicon_log = new G4LogicalVolume(si_face, mat_Si, "fRCMPSilicon_log");
		fRCMPSilicon_log->SetVisAttributes(si_face_vis_att);
		//std::cout<<"HELLO fRCMPSilicon_log\n";
	}

	return 1;
}// end ::end SetUpFace
