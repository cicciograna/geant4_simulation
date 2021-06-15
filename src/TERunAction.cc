/// \file TERunAction.cc
/// \brief Implementation of the TERunAction class

#include "TERunAction.hh"

#include "TEAnalysis.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TERunAction::TERunAction()
 : G4UserRunAction()
{
# if defined(ISO_CURVES_GENERATOR) || defined(ENERGIES_RANGE_ISO_CURVES_GENERATOR) ||defined(FIXED_PARTICLES)
  // set printing event number per each 100 events
  G4RunManager::GetRunManager()->SetPrintProgress(1000);

  // Create analysis manager
  // The choice of analysis technology is done via selectin of a namespace
  // in B4Analysis.hh
  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;

  // Create directories
  analysisManager->SetVerboseLevel(0);
  analysisManager->SetNtupleMerging(true);

  // Creating ntuple
  //
  analysisManager->CreateNtuple("tends", "Track info");
  analysisManager->CreateNtupleDColumn("ionn"         );
  analysisManager->CreateNtupleDColumn("mass"         );
  analysisManager->CreateNtupleDColumn("charge"       );
  analysisManager->CreateNtupleDColumn("tof"          );
  analysisManager->CreateNtupleDColumn("init_x"       );
  analysisManager->CreateNtupleDColumn("init_y"       );
  analysisManager->CreateNtupleDColumn("init_z"       );
  analysisManager->CreateNtupleDColumn("init_azm"     );
  analysisManager->CreateNtupleDColumn("init_elev"    );
  analysisManager->CreateNtupleDColumn("init_dircos"  );
  analysisManager->CreateNtupleDColumn("init_vx"      );
  analysisManager->CreateNtupleDColumn("init_vy"      );
  analysisManager->CreateNtupleDColumn("init_vz"      );
  analysisManager->CreateNtupleDColumn("init_vtot"    );
  analysisManager->CreateNtupleDColumn("init_vrad"    );
  analysisManager->CreateNtupleDColumn("init_vlong"   );
  analysisManager->CreateNtupleDColumn("init_px"      );
  analysisManager->CreateNtupleDColumn("init_py"      );
  analysisManager->CreateNtupleDColumn("init_pz"      );
  analysisManager->CreateNtupleDColumn("init_ptot"    );
  analysisManager->CreateNtupleDColumn("init_prad"    );
  analysisManager->CreateNtupleDColumn("init_plong"   );
  analysisManager->CreateNtupleDColumn("init_ke"      );
  analysisManager->CreateNtupleDColumn("splat_x"      );
  analysisManager->CreateNtupleDColumn("splat_y"      );
  analysisManager->CreateNtupleDColumn("splat_z"      );
  analysisManager->CreateNtupleDColumn("splat_azm"    );
  analysisManager->CreateNtupleDColumn("splat_elev"   );
  analysisManager->CreateNtupleDColumn("splat_dircos" );
  analysisManager->CreateNtupleDColumn("splat_vx"     );
  analysisManager->CreateNtupleDColumn("splat_vy"     );
  analysisManager->CreateNtupleDColumn("splat_vz"     );
  analysisManager->CreateNtupleDColumn("splat_vtot"   );
  analysisManager->CreateNtupleDColumn("splat_vrad"   );
  analysisManager->CreateNtupleDColumn("splat_vlong"  );
  analysisManager->CreateNtupleDColumn("splat_px"     );
  analysisManager->CreateNtupleDColumn("splat_py"     );
  analysisManager->CreateNtupleDColumn("splat_pz"     );
  analysisManager->CreateNtupleDColumn("splat_ptot"   );
  analysisManager->CreateNtupleDColumn("splat_prad"   );
  analysisManager->CreateNtupleDColumn("splat_plong"  );
  analysisManager->CreateNtupleDColumn("splat_ke"     );
  analysisManager->CreateNtupleDColumn("splat_radius" );
  analysisManager->CreateNtupleDColumn("run_id"       );
  analysisManager->CreateNtupleIColumn("ontarget"     );
  analysisManager->CreateNtupleIColumn("onmcpplane"   );
  analysisManager->FinishNtuple();

# endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TERunAction::~TERunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TERunAction::BeginOfRunAction(const G4Run*)
{
# if defined(ISO_CURVES_GENERATOR) || defined(ENERGIES_RANGE_ISO_CURVES_GENERATOR) || defined(FIXED_PARTICLES)
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);

    // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  # ifdef ISO_CURVES_GENERATOR
//  G4String fileName = "../iso_curves_generator_theta1_phi1_KE1_150_1eV_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
  G4String fileName = "../iso_curves_generator_theta1_phi1_KE1_150_1eV_MuFe_WithSleeves_WithPMTSh_450mm_pa20210111.root";
//  G4String fileName = "../iso_curves_generator_theta1_phi1_KE1_150_1eV_MuFe_WithSleeves_NoPMTSh_pa20210111.root";

  # elif defined(ENERGIES_RANGE_ISO_CURVES_GENERATOR)
//  G4String fileName = "../iso_curves_generator_theta0.5_phi1_StandardErange_5eV_MuFe_WithSleeves_NoPMTSh_pa20210111.root";
//  G4String fileName = "../iso_curves_generator_theta0.5_phi1_StandardErange_5eV_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
  G4String fileName = "../iso_curves_generator_theta0.5_phi1_StandardErange5eV_dE0.25eV_MuFe_WithSleeves_WithPMTSh_pa20210111.root";

  # elif defined(FIXED_PARTICLES)
//  G4String fileName = "../output_10k_FixedElev30_GaussSrc_20201029.root"; // Default
//  G4String fileName = "../output_DefaultFieldMap_CrossSrc_z_ne414720_ds0.1_dth5.0_dph5.0.root"; // Default
//  G4String fileName = "../output_HiEnergyAugers_PointEnergies_GaussSrc_20201104.root"; // Default
//  G4String fileName = "../output_ShakeOffs_GaussSrc_20201122.root"; // Default
//  G4String fileName = "../output_100kEVENTS_realProbabilities_GaussSrc_20201202.root"; // Real probabilities, Gaussian source
//  G4String fileName = "../output_NoSleevesWithShield_CrossSrc_z_ne414720_ds0.1_dth5.0_dph5.0.root"; // NoSleeves_WithShield

//  G4String fileName = "../output_RunID1_fixTheta60_minE20-maxE130_100k_FeFe_NoSleeves_NoPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID2_fixTheta60_minE20-maxE130_100k_FeFe_NoSleeves_WithPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID3_fixTheta60_minE20-maxE130_100k_MuFe_NoSleeves_NoPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID4_fixTheta60_minE20-maxE130_100k_MuFe_NoSleeves_WithPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID5_fixTheta60_minE20-maxE130_100k_MuFe_WithSleeves_NoPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID6_fixTheta60_minE20-maxE130_100k_MuFe_WithSleeves_WithPMTSh_PointSrc_20201214.root";

//  G4String fileName = "../output_RunID1_fixTheta30_standardE_100k_FeFe_NoSleeves_NoPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID2_fixTheta30_standardE_100k_FeFe_NoSleeves_WithPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID3_fixTheta30_standardE_100k_MuFe_NoSleeves_NoPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID4_fixTheta30_standardE_100k_MuFe_NoSleeves_WithPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID5_fixTheta30_standardE_100k_MuFe_WithSleeves_NoPMTSh_PointSrc_20201214.root";
//  G4String fileName = "../output_RunID6_fixTheta30_standardE_100k_MuFe_WithSleeves_WithPMTSh_PointSrc_20201214.root";

//  G4String fileName = "../output_RunID18100_GaussSrc_standardE_100k_FeFe_NoSleeves_NoPMTSh_GaussSrc_20210106.root";
//  G4String fileName = "../output_RunID28100_GaussSrc_standardE_100k_FeFe_NoSleeves_WithPMTSh_GaussSrc_20210106.root";
//  G4String fileName = "../output_RunID38100_GaussSrc_standardE_100k_MuFe_NoSleeves_NoPMTSh_GaussSrc_20210106.root";
//  G4String fileName = "../output_RunID338100_GaussSrc_standardE_100k_MuFe_NoSleeves_NoPMTSh_GaussSrc_LinInt_20210106.root";
//  G4String fileName = "../output_RunID48100_GaussSrc_standardE_100k_MuFe_NoSleeves_WithPMTSh_GaussSrc_20210106.root";
//  G4String fileName = "../output_RunID58100_GaussSrc_standardE_100k_MuFe_WithSleeves_NoPMTSh_GaussSrc_20210106_alternative.root";
//   G4String fileName = "../output_RunID68100_GaussSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_GaussSrc_20210106_alternative.root";


  // Fixed elevation
//  G4String fileName = "../output_RunID540030_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_NoPMTSh_pa20210111.root"  ;
//   G4String fileName = "../output_RunID640030_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
//   G4String fileName = "../output_RunID650030_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
//   G4String fileName = "../output_RunID660030_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";

//  G4String fileName = "../output_RunID54002413_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_NoPMTSh_pa20210223.root"  ;
//  G4String fileName = "../output_RunID64002413_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210223.root"  ;
//  G4String fileName = "../output_RunID65002413_fixTheta30_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210223.root"  ;

//   G4String fileName = "../output_RunID64502413_fixTheta2413_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210310.root"  ;
//   G4String fileName = "../output_RunID64502413_fixTheta3131_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210310.root"  ;
//   G4String fileName = "../output_RunID64502413_fixTheta4312_PointSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210310.root"  ;


  // Gaussian source
//   G4String fileName = "../output_RunID5400100_GaussSrc_standardE_100k_MuFe_WithSleeves_NoPMTSh_pa20210111.root"  ;
//   G4String fileName = "../output_RunID6400100_GaussSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
  G4String fileName = "../output_RunID6450100_GaussSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
//   G4String fileName = "../output_RunID6500100_GaussSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";
//   G4String fileName = "../output_RunID6600100_GaussSrc_standardE_100k_MuFe_WithSleeves_WithPMTSh_pa20210111.root";





  # else
  G4String fileName = "../dummy_output.root";
  # endif
  
  G4cout << "Output saved to file " << fileName << G4endl;
  analysisManager->OpenFile(fileName);
# endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TERunAction::EndOfRunAction(const G4Run* )
{
# if defined(ISO_CURVES_GENERATOR) || defined(ENERGIES_RANGE_ISO_CURVES_GENERATOR) || defined(FIXED_PARTICLES)
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
# endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
