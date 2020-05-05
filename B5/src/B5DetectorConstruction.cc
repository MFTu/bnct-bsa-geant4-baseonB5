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
/// \file B5DetectorConstruction.cc
/// \brief Implementation of the B5DetectorConstruction class

#include "B5DetectorConstruction.hh"
#include "B5MagneticField.hh"
#include "B5CellParameterisation.hh"
#include "B5HodoscopeSD.hh"
#include "B5DriftChamberSD.hh"
#include "B5EmCalorimeterSD.hh"
#include "B5HadCalorimeterSD.hh"

#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4Mag_UsualEqRhs.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4MaterialTable.hh"
#include "G4NistManager.hh"

#include "G4VSolid.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4PVReplica.hh"
#include "G4UserLimits.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

//added headers
#include "G4SubtractionSolid.hh"
#include "G4Cons.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal B5MagneticField* B5DetectorConstruction::fMagneticField = 0;
G4ThreadLocal G4FieldManager* B5DetectorConstruction::fFieldMgr = 0;
    
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::B5DetectorConstruction()
: G4VUserDetectorConstruction(), 
  fMessenger(nullptr),
  fHodoscope1Logical(nullptr), fHodoscope2Logical(nullptr),
  fWirePlane1Logical(nullptr), fWirePlane2Logical(nullptr),
  fCellLogical(nullptr), fHadCalScintiLogical(nullptr),
  fMagneticLogical(nullptr),
  fVisAttributes(),
  fArmAngle(30.*deg), fArmRotation(nullptr), fSecondArmPhys(nullptr)

  // added for BNCT BSA
  // logicBSAeDetector(nullptr)

{
  fArmRotation = new G4RotationMatrix();
  fArmRotation->rotateY(fArmAngle);
  
  // define commands for this class
  DefineCommands();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B5DetectorConstruction::~B5DetectorConstruction()
{
  delete fArmRotation;
  delete fMessenger;
  
  for (auto visAttributes: fVisAttributes) {
    delete visAttributes;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B5DetectorConstruction::Construct()
{
  // Construct materials
  ConstructMaterials();
  auto air = G4Material::GetMaterial("G4_AIR");
  //auto argonGas = G4Material::GetMaterial("B5_Ar");
  auto argonGas = G4Material::GetMaterial("G4_Ar");
  auto scintillator = G4Material::GetMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  auto csI = G4Material::GetMaterial("G4_CESIUM_IODIDE");
  auto lead = G4Material::GetMaterial("G4_Pb");


  // added materials
  auto vacuum = G4Material::GetMaterial("G4_Galactic");
  auto shifter_mat = G4Material::GetMaterial("MatConTiF3");
  auto fnfilter_mat = G4Material::GetMaterial("G4_Ni");
  auto grshielding_mat = G4Material::GetMaterial("G4_Bi");
  auto thneuabs_mat = G4Material::GetMaterial("G4_Cd");
  auto rside_reflector = G4Material::GetMaterial("G4_Pb");
  auto delimiter_mat = G4Material::GetMaterial("G4_POLYETHYLENE");
  auto tar_mat = G4Material::GetMaterial("G4_Li");
  auto BSA_exit_detector_mat = G4Material::GetMaterial("G4_POLYSTYRENE");





  
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  // geometries --------------------------------------------------------------
  // experimental hall (world volume)
  auto worldSolid 
    = new G4Box("worldBox",10.*m,3.*m,10.*m);
  auto worldLogical
    //= new G4LogicalVolume(worldSolid,air,"worldLogical");
    = new G4LogicalVolume(worldSolid,vacuum,"worldLogical");
  auto worldPhysical
    = new G4PVPlacement(0,G4ThreeVector(),worldLogical,"worldPhysical",0,
                        false,0,checkOverlaps);
  
  // Tube with Local Magnetic field
  
  auto magneticSolid 
    = new G4Tubs("magneticTubs",0.,1.*m,1.*m,0.,360.*deg);

  fMagneticLogical
    = new G4LogicalVolume(magneticSolid, air, "magneticLogical");

  // placement of Tube
  
  G4RotationMatrix* fieldRot = new G4RotationMatrix();
  fieldRot->rotateX(90.*deg);
  /*
  new G4PVPlacement(fieldRot,G4ThreeVector(),fMagneticLogical,
                    "magneticPhysical",worldLogical,
                    false,0,checkOverlaps);
  */

  // set step limit in tube with magnetic field  
  G4UserLimits* userLimits = new G4UserLimits(1*m);
  fMagneticLogical->SetUserLimits(userLimits);
  
  // first arm
  auto firstArmSolid 
    = new G4Box("firstArmBox",1.5*m,1.*m,3.*m);
  auto firstArmLogical
    = new G4LogicalVolume(firstArmSolid,air,"firstArmLogical");
  /*  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,-5.*m),firstArmLogical,
                    "firstArmPhysical",worldLogical,
                    false,0,checkOverlaps);
  */

  
  // second arm
  auto secondArmSolid 
    = new G4Box("secondArmBox",2.*m,2.*m,3.5*m);
  auto secondArmLogical
    = new G4LogicalVolume(secondArmSolid,air,"secondArmLogical");
  auto x = -5.*m * std::sin(fArmAngle);
  auto z = 5.*m * std::cos(fArmAngle);

  /*
  fSecondArmPhys
    = new G4PVPlacement(fArmRotation,G4ThreeVector(x,0.,z),secondArmLogical,
                        "fSecondArmPhys",worldLogical,
                        false,0,checkOverlaps);
  */


  // hodoscopes in first arm
  auto hodoscope1Solid 
    = new G4Box("hodoscope1Box",5.*cm,20.*cm,0.5*cm);
  fHodoscope1Logical
    = new G4LogicalVolume(hodoscope1Solid,scintillator,"hodoscope1Logical");

  /*
  for (auto i=0;i<kNofHodoscopes1;i++) {
      G4double x1 = (i-kNofHodoscopes1/2)*10.*cm;
      new G4PVPlacement(0,G4ThreeVector(x1,0.,-1.5*m),fHodoscope1Logical,
                        "hodoscope1Physical",firstArmLogical,
                        false,i,checkOverlaps);
  }
  */

  
  // drift chambers in first arm
  auto chamber1Solid 
    = new G4Box("chamber1Box",1.*m,30.*cm,1.*cm);
  auto chamber1Logical
    = new G4LogicalVolume(chamber1Solid,argonGas,"chamber1Logical");

  /*
  for (auto i=0;i<kNofChambers;i++) {
    G4double z1 = (i-kNofChambers/2)*0.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z1),chamber1Logical,
                      "chamber1Physical",firstArmLogical,
                      false,i,checkOverlaps);
  }
  */

  
  // "virtual" wire plane
  auto wirePlane1Solid 
    = new G4Box("wirePlane1Box",1.*m,30.*cm,0.1*mm);
  fWirePlane1Logical
    = new G4LogicalVolume(wirePlane1Solid,argonGas,"wirePlane1Logical");

  /*  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWirePlane1Logical,
                    "wirePlane1Physical",chamber1Logical,
                    false,0,checkOverlaps);
  */


  // hodoscopes in second arm
  auto hodoscope2Solid 
    = new G4Box("hodoscope2Box",5.*cm,20.*cm,0.5*cm);
  fHodoscope2Logical
    = new G4LogicalVolume(hodoscope2Solid,scintillator,"hodoscope2Logical");

  /*
  for (auto i=0;i<kNofHodoscopes2;i++) {
      G4double x2 = (i-kNofHodoscopes2/2)*10.*cm;
      new G4PVPlacement(0,G4ThreeVector(x2,0.,0.),fHodoscope2Logical,
                        "hodoscope2Physical",secondArmLogical,
                        false,i,checkOverlaps);
  }
  */

  
  // drift chambers in second arm
  auto chamber2Solid 
    = new G4Box("chamber2Box",1.5*m,30.*cm,1.*cm);
  auto chamber2Logical
    = new G4LogicalVolume(chamber2Solid,argonGas,"chamber2Logical");
  
  /*
  for (auto i=0;i<kNofChambers;i++) {
    G4double z2 = (i-kNofChambers/2)*0.5*m - 1.5*m;
    new G4PVPlacement(0,G4ThreeVector(0.,0.,z2),chamber2Logical,
                      "chamber2Physical",secondArmLogical,
                      false,i,checkOverlaps);
  }
  */

  
  // "virtual" wire plane
  auto wirePlane2Solid 
    = new G4Box("wirePlane2Box",1.5*m,30.*cm,0.1*mm);
  fWirePlane2Logical
    = new G4LogicalVolume(wirePlane2Solid,argonGas,"wirePlane2Logical");

  /*  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,0.),fWirePlane2Logical,
                    "wirePlane2Physical",chamber2Logical,
                    false,0,checkOverlaps);
  */


  // CsI calorimeter
  auto emCalorimeterSolid 
    = new G4Box("EMcalorimeterBox",1.5*m,30.*cm,15.*cm);
  auto emCalorimeterLogical
    = new G4LogicalVolume(emCalorimeterSolid,csI,"EMcalorimeterLogical");

  /*  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*m),emCalorimeterLogical,
                    "EMcalorimeterPhysical",secondArmLogical,
                    false,0,checkOverlaps);
  */


  // EMcalorimeter cells
  /*
  auto cellSolid 
    = new G4Box("cellBox",7.5*cm,7.5*cm,15.*cm);
  fCellLogical
    = new G4LogicalVolume(cellSolid,csI,"cellLogical");
  G4VPVParameterisation* cellParam = new B5CellParameterisation();
  new G4PVParameterised("cellPhysical",fCellLogical,emCalorimeterLogical,
                        kXAxis,kNofEmCells,cellParam);
  */


  // hadron calorimeter
  auto hadCalorimeterSolid
    = new G4Box("HadCalorimeterBox",1.5*m,30.*cm,50.*cm);
  auto hadCalorimeterLogical
    = new G4LogicalVolume(hadCalorimeterSolid,lead,"HadCalorimeterLogical");

  /*  
  new G4PVPlacement(0,G4ThreeVector(0.,0.,3.*m),hadCalorimeterLogical,
                    "HadCalorimeterPhysical",secondArmLogical,
                    false,0,checkOverlaps);
  */


  // hadron calorimeter column
  auto HadCalColumnSolid
    = new G4Box("HadCalColumnBox",15.*cm,30.*cm,50.*cm);
  auto HadCalColumnLogical
    = new G4LogicalVolume(HadCalColumnSolid,lead,"HadCalColumnLogical");
  new G4PVReplica("HadCalColumnPhysical",HadCalColumnLogical,
                  hadCalorimeterLogical,kXAxis,kNofHadColumns,30.*cm);
  
  // hadron calorimeter cell
  auto HadCalCellSolid
    = new G4Box("HadCalCellBox",15.*cm,15.*cm,50.*cm);
  auto HadCalCellLogical
    = new G4LogicalVolume(HadCalCellSolid,lead,"HadCalCellLogical");
  new G4PVReplica("HadCalCellPhysical",HadCalCellLogical,
                  HadCalColumnLogical,kYAxis,kNofHadRows,30.*cm);
  
  // hadron calorimeter layers
  auto HadCalLayerSolid
    = new G4Box("HadCalLayerBox",15.*cm,15.*cm,2.5*cm);
  auto HadCalLayerLogical
    = new G4LogicalVolume(HadCalLayerSolid,lead,"HadCalLayerLogical");
  new G4PVReplica("HadCalLayerPhysical",HadCalLayerLogical,
                  HadCalCellLogical,kZAxis,kNofHadCells,5.*cm);
  
  // scintillator plates
  auto HadCalScintiSolid
    = new G4Box("HadCalScintiBox",15.*cm,15.*cm,0.5*cm);
  fHadCalScintiLogical
    = new G4LogicalVolume(HadCalScintiSolid,scintillator,
                          "HadCalScintiLogical");

  /*                        
  new G4PVPlacement(0,G4ThreeVector(0.,0.,2.*cm),fHadCalScintiLogical,
                    "HadCalScintiPhysical",HadCalLayerLogical,
                    false,0,checkOverlaps);
  */


  
  // visualization attributes ------------------------------------------------
  
  // all attributes for original B5 design
  auto visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  worldLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.9));   // LightGray
  fMagneticLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  visAttributes->SetVisibility(false);
  firstArmLogical->SetVisAttributes(visAttributes);
  secondArmLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.8888,0.0,0.0));
  fHodoscope1Logical->SetVisAttributes(visAttributes);
  fHodoscope2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0,1.0,0.0));
  chamber1Logical->SetVisAttributes(visAttributes);
  chamber2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0,0.8888,0.0));
  visAttributes->SetVisibility(false);
  fWirePlane1Logical->SetVisAttributes(visAttributes);
  fWirePlane2Logical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.8888,0.8888,0.0));
  visAttributes->SetVisibility(false);
  emCalorimeterLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  /*
  visAttributes = new G4VisAttributes(G4Colour(0.9,0.9,0.0));
  fCellLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  */


  visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  hadCalorimeterLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);
  
  visAttributes = new G4VisAttributes(G4Colour(0.0, 0.0, 0.9));
  visAttributes->SetVisibility(false);
  HadCalColumnLogical->SetVisAttributes(visAttributes);
  HadCalCellLogical->SetVisAttributes(visAttributes);
  HadCalLayerLogical->SetVisAttributes(visAttributes);
  fHadCalScintiLogical->SetVisAttributes(visAttributes);
  fVisAttributes.push_back(visAttributes);



  //******************************************************************
  // added geometries for BNCT BSA


  // build spectrum shifter material TiF3

  G4double target_thick  = 0.001*cm;

  G4ThreeVector pos_shifter = G4ThreeVector(0, 0, (14/2)*cm);
  //beam(10cm diameter, 82cm long) plus target(2cm lithium)Tub
  G4double BeaminnerRadius = 0.*cm; G4double BeamouterRadius = 10/2.*cm; G4double Beamhz = (84+target_thick)/2.*cm;
  G4double BeamstartAngle = 0.*deg; G4double BeamspanningAngle = 360.*deg;
  G4Tubs* BeamPlusTar
     = new G4Tubs("BeamPlusTar",
               BeaminnerRadius,
               BeamouterRadius,
               Beamhz,
               BeamstartAngle,   
               BeamspanningAngle);


  //positon and orientation of BeamTarget "relative to shiter cone" 
  G4RotationMatrix* yBeamRot = new G4RotationMatrix; // Rotates X and Z axes only 
  yBeamRot->rotateY(-M_PI/4.*rad); // Rotates 45 degrees 
  // yBeamRot->rotateY(0.*rad); // Rotates 45 degrees 
   // G4ThreeVector zBeamTrans(0, 0, 0);   
   //G4ThreeVector zBeamTrans(((((-84-2)/2.)+1.)/sqrt(2))*cm, 0, (-((14+82+2)/2.)+1.)*cm); 
   
  G4ThreeVector zxBeamTrans(((-((14+84+target_thick)/2.)+(target_thick)/2.)/(sqrt(2)))*cm, 0, (-((84+target_thick)/2.)+(target_thick)/2.)*cm);
        
  // Conical section for shifter       
  G4double shifterwosubeam_rmina =  0.*cm, shifterwosubeam_rmaxb = (6+43/3.)*cm;
  G4double shifterwosubeam_rminb =  0.*cm, shifterwosubeam_rmaxa = 25.*cm;
  G4double shifterwosubeam_hz = (14/2)*cm;
  G4double shifterwosubeam_phimin = 0.*deg, shifterwosubeam_phimax = 360.*deg;
  G4Cons* solidShifterwosubeam =    
    new G4Cons("Shifterwosubeam", 
    shifterwosubeam_rmina, shifterwosubeam_rmaxa, shifterwosubeam_rminb, shifterwosubeam_rmaxb, 
    shifterwosubeam_hz, shifterwosubeam_phimin, shifterwosubeam_phimax);

  //Conical section for shifter subtract beamplustatget
  G4SubtractionSolid* solidShifter =
    new G4SubtractionSolid("solidShifter", solidShifterwosubeam, BeamPlusTar, yBeamRot, zxBeamTrans);


                      
  G4LogicalVolume* logicShifter =                         
    new G4LogicalVolume(solidShifter,         //its solid
                        shifter_mat,          //its material
                        "Shifter");           //its name
  logicShifter->SetVisAttributes(new G4VisAttributes(G4Colour::Grey()));
               
  new G4PVPlacement(0,                       //no rotation
                    pos_shifter,                    //at position
                    logicShifter,             //its logical volume
                    "Shifter",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  // Fast_Neutron_Filter
  //  
  //G4Material* fnfilter_mat = nist->FindOrBuildMaterial("G4_Ni"); //Ni
  G4ThreeVector pos_fnfilter = G4ThreeVector(0, 0, (14+30/2)*cm);
        
  // Conical section for fast neutron filter       
  G4double fnfilter_rmina =  0.*cm, fnfilter_rmaxb = (6+13/3.)*cm;
  G4double fnfilter_rminb =  0.*cm, fnfilter_rmaxa = (6+43/3.)*cm;
  G4double fnfilter_hz = (30/2)*cm;
  G4double fnfilter_phimin = 0.*deg, fnfilter_phimax = 360.*deg;
  G4Cons* solidFnfilter =    
    new G4Cons("Fnfilter", 
    fnfilter_rmina, fnfilter_rmaxa, fnfilter_rminb, fnfilter_rmaxb, fnfilter_hz,
    fnfilter_phimin, fnfilter_phimax);
                      
  G4LogicalVolume* logicFnfilter =                         
    new G4LogicalVolume(solidFnfilter,         //its solid
                        fnfilter_mat,          //its material
                        "Fnfilter");           //its name
  logicFnfilter->SetVisAttributes(new G4VisAttributes(G4Colour::Yellow()));
               
  new G4PVPlacement(0,                       //no rotation
                    pos_fnfilter,                    //at position
                    logicFnfilter,             //its logical volume
                    "Fnfilter",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  // Gamma_Ray_Sheilding
  //  
  //G4Material* grshielding_mat = nist->FindOrBuildMaterial("G4_Bi"); //Bi
  G4ThreeVector pos_grshielding = G4ThreeVector(0, 0, (14+30+3.5/2)*cm);
        
  // Conical section for Gamma_Ray_Sheilding       
  G4double grshielding_rmina =  0.*cm, grshielding_rmaxa = (6+13/3.)*cm;
  G4double grshielding_rminb =  0.*cm, grshielding_rmaxb = (6+9.5/3.)*cm;
  G4double grshielding_hz = (3.5/2)*cm;
  G4double grshielding_phimin = 0.*deg, grshielding_phimax = 360.*deg;
  G4Cons* solidGrshielding =    
    new G4Cons("Grshielding", 
    grshielding_rmina, grshielding_rmaxa, grshielding_rminb, grshielding_rmaxb, grshielding_hz,
    grshielding_phimin, grshielding_phimax);
                      
  G4LogicalVolume* logicGrshielding =                         
    new G4LogicalVolume(solidGrshielding,         //its solid
                        grshielding_mat,          //its material
                        "Grshielding");           //its name
  logicGrshielding->SetVisAttributes(new G4VisAttributes(G4Colour(255/255.,183/255.,169/255.)));

  new G4PVPlacement(0,                       //no rotation
                    pos_grshielding,                    //at position
                    logicGrshielding,             //its logical volume
                    "Grshielding",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking




  // thermal_neutron_absorber
  //  
  //G4Material* thneuabs_mat = nist->FindOrBuildMaterial("G4_Cd"); //Cd
  G4ThreeVector pos_thneuabs = G4ThreeVector(0, 0, (14+30+3.5+0.1/2)*cm);
        
  // Conical section for thermal_neutron_absorber      
  G4double thneuabs_rmina =  0.*cm, thneuabs_rmaxb = (6+9.4/3.)*cm;
  G4double thneuabs_rminb =  0.*cm, thneuabs_rmaxa = (6+9.5/3.)*cm;
  G4double thneuabs_hz = (0.1/2)*cm;
  G4double thneuabs_phimin = 0.*deg, thneuabs_phimax = 360.*deg;
  G4Cons* solidThneuabs =    
    new G4Cons("Thneuabs", 
    thneuabs_rmina, thneuabs_rmaxa, thneuabs_rminb, thneuabs_rmaxb, thneuabs_hz,
    thneuabs_phimin, thneuabs_phimax);
                      
  G4LogicalVolume* logicThneuabs =                         
    new G4LogicalVolume(solidThneuabs,         //its solid
                        thneuabs_mat,          //its material
                        "Thneuabs");           //its name
  logicThneuabs->SetVisAttributes(new G4VisAttributes(G4Colour::Red()));

               
  new G4PVPlacement(0,                       //no rotation
                    pos_thneuabs,                    //at position
                    logicThneuabs,             //its logical volume
                    "Thneuabs",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  
  // full cone minus poly-Li(delimiter)
  //  
  G4ThreeVector pos_fullcone_minus_polyLi = G4ThreeVector(0, 0, ((14+30+3.5+0.1+9.4-2.0)/2)*cm);
        
  // Conical section minus poly-Li(delimiter)       
  G4double fullconemp_rmina =  0.*cm, fullconemp_rmaxa = (50/2.)*cm;
  G4double fullconemp_rminb =  0.*cm, fullconemp_rmaxb = (6+2/3)*cm;
  G4double fullconemp_hz = ((14+30+3.5+0.1+9.4-2.0)/2)*cm;
  G4double fullconemp_phimin = 0.*deg, fullconemp_phimax = 360.*deg;
  G4Cons* solidFullconemp =    
    new G4Cons("Fullconemp", 
    fullconemp_rmina, fullconemp_rmaxa, fullconemp_rminb, fullconemp_rmaxb, fullconemp_hz,
    fullconemp_phimin, fullconemp_phimax);

  // right side Pd(reflector) minus poly-Li(delimiter)
  //  
  //G4Material* rside_reflector = nist->FindOrBuildMaterial("G4_Pb"); //Pb
  G4ThreeVector rside_refle_minus_polyLi = G4ThreeVector(0, 0, ((14+30+3.5+0.1+9.4-2.0)/2)*cm);
        
  // right side box minus poly-Li(delimiter) length
  G4double rsbmp_sizeXY = 70.*cm;
  G4double rsbmp_sizeZ  = (14+30+3.5+0.1+9.4-2.0)*cm;      

  G4Box* solidRsbmp =    
   new G4Box("RsideEnvelopemp",                    //its name
       0.5*rsbmp_sizeXY, 0.5*rsbmp_sizeXY, 0.5*rsbmp_sizeZ); //its size

  G4SubtractionSolid* RsbsubtconeminusPL =
   new G4SubtractionSolid("RsbsubtconeminusPL", solidRsbmp, solidFullconemp);


  G4LogicalVolume* logicRsbsubtconeminusPL =                         
    new G4LogicalVolume(RsbsubtconeminusPL,         //its solid
                        rside_reflector,          //its material
                        "RsbsubtconeminusPL");               //its name
  logicRsbsubtconeminusPL->SetVisAttributes(new G4VisAttributes(G4Colour(204/255.,0/255.,204/255.)));
              
  new G4PVPlacement(0,                       //no rotation
                    rside_refle_minus_polyLi,                    //at position
                    logicRsbsubtconeminusPL,             //its logical volume
                    "RsbsubtconeminusPL",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  // last section poly-Li(delimiter) (polyethylene)
  //  
  //G4Material* delimiter_mat = nist->FindOrBuildMaterial("G4_POLYETHYLENE"); //POLYETHYLENE
  G4ThreeVector pos_delimiter = G4ThreeVector(0, 0, (14+30+3.5+0.1+9.4-2.0+(2.0)/2)*cm);
        
  // poly-Li(delimiter) length
  G4double delimiter_sizeXY = 70.*cm;
  G4double delimiter_sizeZ  = 2.*cm;   

  // last section cone 
  //G4ThreeVector pos_fullcone_minus_polyLi = G4ThreeVector(0, 0, ((14+30+3.5+0.1+9.4-2.0)/2)*cm);
      
  G4double lasectioncone_rmina =  0.*cm, lasectioncone_rmaxb = (12/2.)*cm;
  G4double lasectioncone_rminb =  0.*cm, lasectioncone_rmaxa = (6+2/3)*cm;
  G4double lasectioncone_hz = (2/2.)*cm;
  G4double lasectioncone_phimin = 0.*deg, lasectioncone_phimax = 360.*deg;
  G4Cons* solidLasectioncone =    
    new G4Cons("Lasectioncone", 
    lasectioncone_rmina, lasectioncone_rmaxa, lasectioncone_rminb, lasectioncone_rmaxb, lasectioncone_hz,
    lasectioncone_phimin, lasectioncone_phimax);  


  G4Box* solidlastlid =    
   new G4Box("lastlid",                    //its name
       0.5*delimiter_sizeXY, 0.5*delimiter_sizeXY, 0.5*delimiter_sizeZ); //its size

  G4SubtractionSolid* DelimiterLid =
   new G4SubtractionSolid("DelimiterLid", solidlastlid, solidLasectioncone);


  G4LogicalVolume* logicDelimiterLid =                         
    new G4LogicalVolume(DelimiterLid,         //its solid
                        delimiter_mat,          //its material
                        "DelimiterLid");               //its name
  logicDelimiterLid->SetVisAttributes(new G4VisAttributes(G4Colour(255/255.,0/255.,127/255.)));
              
  new G4PVPlacement(0,                       //no rotation
                    pos_delimiter,                    //at position
                    logicDelimiterLid,             //its logical volume
                    "DelimiterLid",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



  // left side reflector
  
  G4ThreeVector lside_refle_minus_beam = G4ThreeVector(0, 0, (-(14+30+3.5+0.1+9.4)/2)*cm);
  //positon and orientation of BeamTarget relative to left box 
  G4RotationMatrix* LsideyBeamRot = new G4RotationMatrix; // Rotates X and Z axes only 
  LsideyBeamRot->rotateY(-(M_PI/4.)*rad); // Rotates 45 degrees 
    
  G4ThreeVector LsidezxBeamTrans(-(42./(sqrt(2)))*cm, 0, -1.5*cm);    

  // left side box 
  G4double leftsb_sizeXY = 70.*cm;
  G4double leftsb_sizeZ  = (14+30+3.5+0.1+9.4)*cm;      

  G4Box* solidLeftsb =    
   new G4Box("Leftsb",                    //its name
       0.5*leftsb_sizeXY, 0.5*leftsb_sizeXY, 0.5*leftsb_sizeZ); //its size

  //Left side box subtract beamplustatget
  G4SubtractionSolid* LeftsbminusBeam =
    new G4SubtractionSolid("LeftsbminusBeam", solidLeftsb, BeamPlusTar, LsideyBeamRot, LsidezxBeamTrans);
     


  G4LogicalVolume* logicLeftsbminusBeam =                         
    new G4LogicalVolume(LeftsbminusBeam,         //its solid
                        rside_reflector,          //its material
                        "LeftsbminusBeam");               //its name
  logicLeftsbminusBeam->SetVisAttributes(new G4VisAttributes(G4Colour(204/255.,0/255.,204/255.)));
              
  new G4PVPlacement(0,                       //no rotation
                    lside_refle_minus_beam,                    //at position
                    logicLeftsbminusBeam,             //its logical volume
                    "LeftsbminusBeam",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking




  //proton to neutron conversion target

  //G4Material* tar_mat = nist->FindOrBuildMaterial("G4_Li"); //Li
  G4ThreeVector pos_tar = G4ThreeVector(0, 0, 0);

  G4double TargetinnerRadius = 0.*cm; G4double TargetouterRadius = 10/2.*cm; 
  G4double Targethz = (target_thick)/2.*cm;
  G4double TargetstartAngle = 0.*deg; G4double TargetspanningAngle = 360.*deg;
  G4Tubs* solidPNConversionTarget
     = new G4Tubs("PNConversionTarget",
               TargetinnerRadius,
               TargetouterRadius,
               Targethz,
               TargetstartAngle,
               TargetspanningAngle);   


                      
  G4LogicalVolume* logicPNConversionTarget =                         
    new G4LogicalVolume(solidPNConversionTarget,         //its solid
                        tar_mat,          //its material
                        "PNConversionTarget");           //its name
  logicPNConversionTarget->SetVisAttributes(new G4VisAttributes(G4Colour::Green()));


  G4RotationMatrix* yTarRot = new G4RotationMatrix; // Rotates X and Z axes only 
   yTarRot->rotateY(-M_PI/4.*rad); // Rotates 45 degrees 
               
  new G4PVPlacement(yTarRot,                       //rotation
                    pos_tar,                    //at position
                    logicPNConversionTarget,             //its logical volume
                    "PNConversionTarget",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

             



  //BSA exit detector

  G4double BSAeDetector_innerRadius = 0.*cm; G4double BSAeDetector_outerRadius = 12/2.*cm; 
  G4double BSAeDetector_hz = 10/2.*cm;
  G4double BSAeDetector_startAngle = 0.*deg; G4double BSAeDetector_spanningAngle = 360.*deg;

  //G4Material* BSA_exit_detector_mat = nist->FindOrBuildMaterial("G4_POLYSTYRENE"); //POLYSTYRENE
  G4ThreeVector pos_BSA_exit_detector = G4ThreeVector(0, 0, ((14+30+3.5+0.1+9.4+((10.)*6.8))/2)*cm); //??


/*
  G4Tubs* solidBSAeDetector
     = new G4Tubs("BSAeDetector",
               BSAeDetector_innerRadius,
               BSAeDetector_outerRadius,
               BSAeDetector_hz,
               BSAeDetector_startAngle,
               BSAeDetector_spanningAngle);
 
*/



/*                      
  G4LogicalVolume* logicBSAeDetector =                         
    new G4LogicalVolume(solidBSAeDetector,         //its solid
                        BSA_exit_detector_mat,          //its material
                        "BSAeDetector");           //its name
  logicBSAeDetector->SetVisAttributes(new G4VisAttributes(G4Colour:: Cyan()));

  //G4RotationMatrix* yTarRot = new G4RotationMatrix; // Rotates X and Z axes only 
  // yTarRot->rotateY(-M_PI/4.*rad); // Rotates 45 degrees 
               
  new G4PVPlacement(0,                       //no rotation
                    pos_BSA_exit_detector,                    //at position
                    logicBSAeDetector,             //its logical volume
                    "BSAeDetector",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
*/

// (important!!!!) use auto instead of G4Tubs* 
  auto solidBSAeDetector
     = new G4Tubs("BSAeDetector",
               BSAeDetector_innerRadius,
               BSAeDetector_outerRadius,
               BSAeDetector_hz,
               BSAeDetector_startAngle,
               BSAeDetector_spanningAngle);  

// (important!!!!) with no G4LogicalVolume*
  fCellLogical =                         
    new G4LogicalVolume(solidBSAeDetector,         //its solid
                        BSA_exit_detector_mat,          //its material
                        "BSAeDetector");           //its name
  fCellLogical->SetVisAttributes(new G4VisAttributes(G4Colour:: Cyan()));

  //G4RotationMatrix* yTarRot = new G4RotationMatrix; // Rotates X and Z axes only 
  // yTarRot->rotateY(-M_PI/4.*rad); // Rotates 45 degrees 
               
  new G4PVPlacement(0,                       //no rotation
                    pos_BSA_exit_detector,                    //at position
                    fCellLogical,             //its logical volume
                    "BSAeDetector",                //its name
                    worldLogical,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking



//    auto cellSolid 
//    = new G4Box("cellBox",7.5*cm,7.5*cm,15.*cm);
//  fCellLogical
//    = new G4LogicalVolume(cellSolid,csI,"cellLogical");
//  G4VPVParameterisation* cellParam = new B5CellParameterisation();
//  new G4PVParameterised("cellPhysical",fCellLogical,emCalorimeterLogical,
//                        kXAxis,kNofEmCells,cellParam);                  



  
  // return the world physical volume ----------------------------------------
  
  return worldPhysical;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructSDandField()
{
  // sensitive detectors -----------------------------------------------------
  auto sdManager = G4SDManager::GetSDMpointer();
  G4String SDname;
  
  // all attributes for original B5 design

  auto hodoscope1 = new B5HodoscopeSD(SDname="/hodoscope1");
  sdManager->AddNewDetector(hodoscope1);
  fHodoscope1Logical->SetSensitiveDetector(hodoscope1);

  auto hodoscope2 = new B5HodoscopeSD(SDname="/hodoscope2");
  sdManager->AddNewDetector(hodoscope2);
  fHodoscope2Logical->SetSensitiveDetector(hodoscope2);
  
  auto chamber1 = new B5DriftChamberSD(SDname="/chamber1");
  sdManager->AddNewDetector(chamber1);
  fWirePlane1Logical->SetSensitiveDetector(chamber1);

  auto chamber2 = new B5DriftChamberSD(SDname="/chamber2");
  sdManager->AddNewDetector(chamber2);
  fWirePlane2Logical->SetSensitiveDetector(chamber2);


  auto hadCalorimeter = new B5HadCalorimeterSD(SDname="/HadCalorimeter");
  sdManager->AddNewDetector(hadCalorimeter);
  fHadCalScintiLogical->SetSensitiveDetector(hadCalorimeter);



  
  // use this for BNCT BSA
  auto emCalorimeter = new B5EmCalorimeterSD(SDname="/EMcalorimeter");
  sdManager->AddNewDetector(emCalorimeter);
  fCellLogical->SetSensitiveDetector(emCalorimeter);

  




  // added for BNCT BSA
  //auto neutronmeter = new B5HadCalorimeterSD(SDname="/neutronmeter"); 
  //sdManager->AddNewDetector(neutronmeter);
  //logicBSAeDetector->SetSensitiveDetector(neutronmeter);



  // magnetic field ----------------------------------------------------------
  fMagneticField = new B5MagneticField();
  fFieldMgr = new G4FieldManager();
  fFieldMgr->SetDetectorField(fMagneticField);
  fFieldMgr->CreateChordFinder(fMagneticField);
  G4bool forceToAllDaughters = true;
  fMagneticLogical->SetFieldManager(fFieldMgr, forceToAllDaughters);
}    

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B5DetectorConstruction::ConstructMaterials()
{
  auto nistManager = G4NistManager::Instance();

  // Air 
  nistManager->FindOrBuildMaterial("G4_AIR");
  
  // Argon gas
  nistManager->FindOrBuildMaterial("G4_Ar");
  // With a density different from the one defined in NIST
  // G4double density = 1.782e-03*g/cm3; 
  // nistManager->BuildMaterialWithNewDensity("B5_Ar","G4_Ar",density);
  // !! cases segmentation fault

  // Scintillator
  // (PolyVinylToluene, C_9H_10)
  nistManager->FindOrBuildMaterial("G4_PLASTIC_SC_VINYLTOLUENE");
  
  // CsI
  nistManager->FindOrBuildMaterial("G4_CESIUM_IODIDE");
  
  // Lead
  nistManager->FindOrBuildMaterial("G4_Pb");
  
  // Vacuum "Galactic"
  // nistManager->FindOrBuildMaterial("G4_Galactic");

  // Vacuum "Air with low density"
  // auto air = G4Material::GetMaterial("G4_AIR");
  // G4double density = 1.0e-5*air->GetDensity();
  // nistManager
  //   ->BuildMaterialWithNewDensity("Air_lowDensity", "G4_AIR", density);



  // add materials for BNCT BSA 
  // added nist material

  //G4Material* vacuum = 
  nistManager->FindOrBuildMaterial("G4_Galactic");

  nistManager->FindOrBuildMaterial("G4_Ni");
  nistManager->FindOrBuildMaterial("G4_Bi");
  nistManager->FindOrBuildMaterial("G4_Cd");
  nistManager->FindOrBuildMaterial("G4_Pb");
  nistManager->FindOrBuildMaterial("G4_POLYETHYLENE");
  nistManager->FindOrBuildMaterial("G4_Li");
  nistManager->FindOrBuildMaterial("4_POLYSTYRENE");


   // build spectrum shifter material TiF3
  G4int ncomponents, natoms;
  G4double TiF3density = 3.4*g/cm3;
  G4Element* elemTi = nistManager->FindOrBuildElement("Ti");
  G4Element* elemF = nistManager->FindOrBuildElement("F");
  G4Material* shifter_mat = new G4Material("MatConTiF3", TiF3density, ncomponents = 2);
  shifter_mat -> AddElement(elemTi, natoms = 1);
  shifter_mat -> AddElement(elemF, natoms = 3);


  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// all attributes for original B5 design

void B5DetectorConstruction::SetArmAngle(G4double val)
{
  if (!fSecondArmPhys) {
      G4cerr << "Detector has not yet been constructed." << G4endl;
      return;
  }
  
  fArmAngle = val;
  *fArmRotation = G4RotationMatrix();  // make it unit vector
  fArmRotation->rotateY(fArmAngle);
  auto x = -5.*m * std::sin(fArmAngle);
  auto z = 5.*m * std::cos(fArmAngle);
  fSecondArmPhys->SetTranslation(G4ThreeVector(x,0.,z));
  
  // tell G4RunManager that we change the geometry
  G4RunManager::GetRunManager()->GeometryHasBeenModified();
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
// all attributes for original B5 design

void B5DetectorConstruction::DefineCommands()
{

  // Define /B5/detector command directory using generic messenger class
  fMessenger = new G4GenericMessenger(this, 
                                      "/B5/detector/", 
                                      "Detector control");


  // armAngle command
  auto& armAngleCmd
    = fMessenger->DeclareMethodWithUnit("armAngle","deg",
                                &B5DetectorConstruction::SetArmAngle, 
                                "Set rotation angle of the second arm.");
  armAngleCmd.SetParameterName("angle", true);
  armAngleCmd.SetRange("angle>=0. && angle<180.");
  armAngleCmd.SetDefaultValue("30.");

}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
