#include "MyTruthAnalysis/TruthAna_NtupleManager.h"

TruthAna_NtupleManager::TruthAna_NtupleManager()
{
}

TruthAna_NtupleManager::TruthAna_NtupleManager(TruthAna_OptionParser *opt, TTree *tree)
{
}


TruthAna_NtupleManager::~TruthAna_NtupleManager()
{
}


void TruthAna_NtupleManager::SetNtupleVars(TruthAna_OptionParser *opt, TTree *tree)
{

  tree->Branch("n_mc_event_weight",&n_mc_event_weight);

  //TRUTH1 collections
  tree->Branch("n_neutrinos_pt" , &n_neutrinos_pt      );
  tree->Branch("n_neutrinos_eta", &n_neutrinos_eta     );
  tree->Branch("n_neutrinos_phi", &n_neutrinos_phi     );
  tree->Branch("n_neutrinos_m"  , &n_neutrinos_m       );
  tree->Branch("n_neutrinos_n"  , &n_neutrinos_n       );
  tree->Branch("n_muons_pt"     , &n_muons_pt          );
  tree->Branch("n_muons_eta"    , &n_muons_eta         );
  tree->Branch("n_muons_phi"    , &n_muons_phi         );
  tree->Branch("n_muons_m"      , &n_muons_m           );
  tree->Branch("n_muons_n"      , &n_muons_n           );
  tree->Branch("n_electrons_pt" , &n_electrons_pt      );
  tree->Branch("n_electrons_eta", &n_electrons_eta     );
  tree->Branch("n_electrons_phi", &n_electrons_phi     );
  tree->Branch("n_electrons_m"  , &n_electrons_m       );
  tree->Branch("n_electrons_n"  , &n_electrons_n       );
  tree->Branch("n_taus_pt"      , &n_taus_pt           );                                                       
  tree->Branch("n_taus_eta"     , &n_taus_eta          );                                                       
  tree->Branch("n_taus_phi"     , &n_taus_phi          );                                                       
  tree->Branch("n_taus_m"       , &n_taus_m            );                                                       
  tree->Branch("n_taus_n"       , &n_taus_n            );                                                       
  tree->Branch("n_wzjets_pt"    , &n_wzjets_pt         );
  tree->Branch("n_wzjets_eta"   , &n_wzjets_eta        );
  tree->Branch("n_wzjets_phi"   , &n_wzjets_phi        );
  tree->Branch("n_wzjets_m"     , &n_wzjets_m          );
  tree->Branch("n_wzjets_n"     , &n_wzjets_n          );

  tree->Branch("n_particle_wboson_pt"    , &n_particle_wboson_pt           );
  tree->Branch("n_particle_wboson_eta"   , &n_particle_wboson_eta          );
  tree->Branch("n_particle_wboson_phi"   , &n_particle_wboson_phi          );
  tree->Branch("n_particle_wboson_m"     , &n_particle_wboson_m            );
  tree->Branch("n_particle_wboson_n"     , &n_particle_wboson_n            );

  tree->Branch("n_parton_decay"     , &n_parton_decay              );
  tree->Branch("n_parton_higgs_pt"  , &n_parton_higgs_pt           );
  tree->Branch("n_parton_higgs_eta" , &n_parton_higgs_eta          );
  tree->Branch("n_parton_higgs_phi" , &n_parton_higgs_phi          );
  tree->Branch("n_parton_higgs_m"   , &n_parton_higgs_m            );
  tree->Branch("n_parton_higgs_n"   , &n_parton_higgs_n            );
  ////  
  tree->Branch("n_hadron_J_B_pt"  , &n_hadron_J_B_pt           );
  tree->Branch("n_hadron_J_B_eta" , &n_hadron_J_B_eta          );
  tree->Branch("n_hadron_J_B_phi" , &n_hadron_J_B_phi          );
  tree->Branch("n_hadron_J_B_m"   , &n_hadron_J_B_m            );
  
  tree->Branch("n_hadron_J_aB_pt"  , &n_hadron_J_aB_pt           );
  tree->Branch("n_hadron_J_aB_eta" , &n_hadron_J_aB_eta          );
  tree->Branch("n_hadron_J_aB_phi" , &n_hadron_J_aB_phi          );
  tree->Branch("n_hadron_J_aB_m"   , &n_hadron_J_aB_m            );
    
  tree->Branch("n_hadron_B_plus_pt"  , &n_hadron_B_plus_pt           );
  tree->Branch("n_hadron_B_plus_eta" , &n_hadron_B_plus_eta          );
  tree->Branch("n_hadron_B_plus_phi" , &n_hadron_B_plus_phi          );
  tree->Branch("n_hadron_B_plus_m"   , &n_hadron_B_plus_m            );
  
  tree->Branch("n_hadron_B_minus_pt"  , &n_hadron_B_minus_pt           );
  tree->Branch("n_hadron_B_minus_eta" , &n_hadron_B_minus_eta          );
  tree->Branch("n_hadron_B_minus_phi" , &n_hadron_B_minus_phi          );
  tree->Branch("n_hadron_B_minus_m"   , &n_hadron_B_minus_m            );    

  tree->Branch("n_hadron_K_B_pt"  , &n_hadron_K_B_pt           );
  tree->Branch("n_hadron_K_B_eta" , &n_hadron_K_B_eta          );
  tree->Branch("n_hadron_K_B_phi" , &n_hadron_K_B_phi          );
  tree->Branch("n_hadron_K_B_m"   , &n_hadron_K_B_m            );
  
  tree->Branch("n_hadron_K_aB_pt"  , &n_hadron_K_aB_pt           );
  tree->Branch("n_hadron_K_aB_eta" , &n_hadron_K_aB_eta          );
  tree->Branch("n_hadron_K_aB_phi" , &n_hadron_K_aB_phi          );
  tree->Branch("n_hadron_K_aB_m"   , &n_hadron_K_aB_m            );
    
  tree->Branch("n_hadron_J_B_Boosted_pt"  , &n_hadron_J_B_Boosted_pt           );
  tree->Branch("n_hadron_J_B_Boosted_eta" , &n_hadron_J_B_Boosted_eta          );
  tree->Branch("n_hadron_J_B_Boosted_phi" , &n_hadron_J_B_Boosted_phi          );
  tree->Branch("n_hadron_J_B_Boosted_m"   , &n_hadron_J_B_Boosted_m            );
  
  tree->Branch("n_hadron_J_aB_Boosted_pt"  , &n_hadron_J_aB_Boosted_pt           );
  tree->Branch("n_hadron_J_aB_Boosted_eta" , &n_hadron_J_aB_Boosted_eta          );
  tree->Branch("n_hadron_J_aB_Boosted_phi" , &n_hadron_J_aB_Boosted_phi          );
  tree->Branch("n_hadron_J_aB_Boosted_m"   , &n_hadron_J_aB_Boosted_m            );  
    
  tree->Branch("n_hadron_J_B_boost_pt"  , &n_hadron_J_B_boost_pt           );
  tree->Branch("n_hadron_J_B_boost_eta" , &n_hadron_J_B_boost_eta          );
  tree->Branch("n_hadron_J_B_boost_phi" , &n_hadron_J_B_boost_phi          );
  tree->Branch("n_hadron_J_B_boost_m"   , &n_hadron_J_B_boost_m            );
  
  tree->Branch("n_hadron_J_aB_boost_pt"  , &n_hadron_J_aB_boost_pt           );
  tree->Branch("n_hadron_J_aB_boost_eta" , &n_hadron_J_aB_boost_eta          );
  tree->Branch("n_hadron_J_aB_boost_phi" , &n_hadron_J_aB_boost_phi          );
  tree->Branch("n_hadron_J_aB_boost_m"   , &n_hadron_J_aB_boost_m            );
    


  tree->Branch("n_parton_u_pt"      , &n_parton_u_pt               );
  tree->Branch("n_parton_u_eta"     , &n_parton_u_eta              );
  tree->Branch("n_parton_u_phi"     , &n_parton_u_phi              );
  tree->Branch("n_parton_u_m"       , &n_parton_u_m                );
  tree->Branch("n_parton_u_n"       , &n_parton_u_n                );
  tree->Branch("n_parton_au_pt"     , &n_parton_au_pt              );
  tree->Branch("n_parton_au_eta"    , &n_parton_au_eta             );
  tree->Branch("n_parton_au_phi"    , &n_parton_au_phi             );
  tree->Branch("n_parton_au_m"      , &n_parton_au_m               );
  tree->Branch("n_parton_au_n"      , &n_parton_au_n               );

  tree->Branch("n_parton_d_pt"      , &n_parton_d_pt               );
  tree->Branch("n_parton_d_eta"     , &n_parton_d_eta              );
  tree->Branch("n_parton_d_phi"     , &n_parton_d_phi              );
  tree->Branch("n_parton_d_m"       , &n_parton_d_m                );
  tree->Branch("n_parton_d_n"       , &n_parton_d_n                );
  tree->Branch("n_parton_ad_pt"     , &n_parton_ad_pt              );
  tree->Branch("n_parton_ad_eta"    , &n_parton_ad_eta             );
  tree->Branch("n_parton_ad_phi"    , &n_parton_ad_phi             );
  tree->Branch("n_parton_ad_m"      , &n_parton_ad_m               );
  tree->Branch("n_parton_ad_n"      , &n_parton_ad_n               );

  tree->Branch("n_parton_c_pt"      , &n_parton_c_pt               );
  tree->Branch("n_parton_c_eta"     , &n_parton_c_eta              );
  tree->Branch("n_parton_c_phi"     , &n_parton_c_phi              );
  tree->Branch("n_parton_c_m"       , &n_parton_c_m                );
  tree->Branch("n_parton_c_n"       , &n_parton_c_n                );
  tree->Branch("n_parton_ac_pt"     , &n_parton_ac_pt              );
  tree->Branch("n_parton_ac_eta"    , &n_parton_ac_eta             );
  tree->Branch("n_parton_ac_phi"    , &n_parton_ac_phi             );
  tree->Branch("n_parton_ac_m"      , &n_parton_ac_m               );
  tree->Branch("n_parton_ac_n"      , &n_parton_ac_n               );

  tree->Branch("n_parton_s_pt"      , &n_parton_s_pt               );
  tree->Branch("n_parton_s_eta"     , &n_parton_s_eta              );
  tree->Branch("n_parton_s_phi"     , &n_parton_s_phi              );
  tree->Branch("n_parton_s_m"       , &n_parton_s_m                );
  tree->Branch("n_parton_s_n"       , &n_parton_s_n                );
  tree->Branch("n_parton_as_pt"     , &n_parton_as_pt              );
  tree->Branch("n_parton_as_eta"    , &n_parton_as_eta             );
  tree->Branch("n_parton_as_phi"    , &n_parton_as_phi             );
  tree->Branch("n_parton_as_m"      , &n_parton_as_m               );
  tree->Branch("n_parton_as_n"      , &n_parton_as_n               );

  tree->Branch("n_parton_top_pt"    , &n_parton_top_pt             );
  tree->Branch("n_parton_top_eta"   , &n_parton_top_eta            );
  tree->Branch("n_parton_top_phi"   , &n_parton_top_phi            );
  tree->Branch("n_parton_top_m"     , &n_parton_top_m              );
  tree->Branch("n_parton_top_n"     , &n_parton_top_n              );
  tree->Branch("n_parton_atop_pt"   , &n_parton_atop_pt            );
  tree->Branch("n_parton_atop_eta"  , &n_parton_atop_eta           );
  tree->Branch("n_parton_atop_phi"  , &n_parton_atop_phi           );
  tree->Branch("n_parton_atop_m"    , &n_parton_atop_m             );
  tree->Branch("n_parton_atop_n"    , &n_parton_atop_n             );

  tree->Branch("n_parton_b_pt"      , &n_parton_b_pt               );
  tree->Branch("n_parton_b_eta"     , &n_parton_b_eta              );
  tree->Branch("n_parton_b_phi"     , &n_parton_b_phi              );
  tree->Branch("n_parton_b_m"       , &n_parton_b_m                );
  tree->Branch("n_parton_b_n"       , &n_parton_b_n                );
  tree->Branch("n_parton_ab_pt"     , &n_parton_ab_pt              );
  tree->Branch("n_parton_ab_eta"    , &n_parton_ab_eta             );
  tree->Branch("n_parton_ab_phi"    , &n_parton_ab_phi             );
  tree->Branch("n_parton_ab_m"      , &n_parton_ab_m               );
  tree->Branch("n_parton_ab_n"      , &n_parton_ab_n               );

  tree->Branch("n_parton_wplus_pt"  , &n_parton_wplus_pt           );
  tree->Branch("n_parton_wplus_eta" , &n_parton_wplus_eta          );
  tree->Branch("n_parton_wplus_phi" , &n_parton_wplus_phi          );
  tree->Branch("n_parton_wplus_m"   , &n_parton_wplus_m            );
  tree->Branch("n_parton_wplus_n"   , &n_parton_wplus_n            );
  tree->Branch("n_parton_wminus_pt" , &n_parton_wminus_pt          );
  tree->Branch("n_parton_wminus_eta", &n_parton_wminus_eta         );
  tree->Branch("n_parton_wminus_phi", &n_parton_wminus_phi         );
  tree->Branch("n_parton_wminus_m"  , &n_parton_wminus_m           );
  tree->Branch("n_parton_wminus_n"  , &n_parton_wminus_n           );

  //bbHtt and ttHtt collections
  tree->Branch("n_parton_ttfromH_pt"       , &n_parton_ttfromH_pt 	              );
  tree->Branch("n_parton_ttfromH_eta"      , &n_parton_ttfromH_eta	              );
  tree->Branch("n_parton_ttfromH_phi"      , &n_parton_ttfromH_phi	              );
  tree->Branch("n_parton_ttfromH_m"        , &n_parton_ttfromH_m  	              );
  tree->Branch("n_parton_ttfromH_n"        , &n_parton_ttfromH_n  	              );
  tree->Branch("n_parton_ttfromH_poscosphi", &n_parton_ttfromH_poscosphi              );
  tree->Branch("n_parton_ttfromH_deltaeta" , &n_parton_ttfromH_deltaeta          );
  tree->Branch("n_parton_ttfromH_deltaphi" , &n_parton_ttfromH_deltaphi          );
  tree->Branch("n_parton_ttfromH_deltar"   , &n_parton_ttfromH_deltar            );

  //ttHtt specific
  tree->Branch("n_parton_ttnotfromH_pt"      , &n_parton_ttnotfromH_pt                );
  tree->Branch("n_parton_ttnotfromH_eta"     , &n_parton_ttnotfromH_eta               );
  tree->Branch("n_parton_ttnotfromH_phi"     , &n_parton_ttnotfromH_phi               );
  tree->Branch("n_parton_ttnotfromH_m"       , &n_parton_ttnotfromH_m                 );
  tree->Branch("n_parton_ttnotfromH_n"       , &n_parton_ttnotfromH_n                 );
  tree->Branch("n_parton_ttnotfromH_deltaeta", &n_parton_ttnotfromH_deltaeta          );
  tree->Branch("n_parton_ttnotfromH_deltaphi", &n_parton_ttnotfromH_deltaphi          );
  tree->Branch("n_parton_ttnotfromH_deltar"  , &n_parton_ttnotfromH_deltar            );

  //tcHbb specific
  tree->Branch("n_parton_tcnotfromH_pt"      , &n_parton_tcnotfromH_pt                );
  tree->Branch("n_parton_tcnotfromH_eta"     , &n_parton_tcnotfromH_eta               );
  tree->Branch("n_parton_tcnotfromH_phi"     , &n_parton_tcnotfromH_phi               );
  tree->Branch("n_parton_tcnotfromH_m"       , &n_parton_tcnotfromH_m                 );
  tree->Branch("n_parton_tcnotfromH_n"       , &n_parton_tcnotfromH_n                 );
  tree->Branch("n_parton_tcnotfromH_deltaeta", &n_parton_tcnotfromH_deltaeta          );
  tree->Branch("n_parton_tcnotfromH_deltaphi", &n_parton_tcnotfromH_deltaphi          );
  tree->Branch("n_parton_tcnotfromH_deltar"  , &n_parton_tcnotfromH_deltar            );

  //txQ specific
  tree->Branch("n_parton_sisterH_eta"        , &n_parton_sisterH_eta );
  tree->Branch("n_parton_sisterH_phi"        , &n_parton_sisterH_phi );
  tree->Branch("n_parton_sisterH_m"          , &n_parton_sisterH_m );
  tree->Branch("n_parton_sisterH_pt"         , &n_parton_sisterH_pt );
  tree->Branch("n_parton_sisterH_pdgid"      , &n_parton_sisterH_pdgid );
  tree->Branch("n_parton_fromothertop_pt"    , &n_parton_fromothertop_pt );
  tree->Branch("n_parton_fromothertop_eta"   , &n_parton_fromothertop_eta );
  tree->Branch("n_parton_fromothertop_phi"   , &n_parton_fromothertop_phi );
  tree->Branch("n_parton_fromothertop_m"     , &n_parton_fromothertop_m );
  tree->Branch("n_parton_fromothertop_pdgid" , &n_parton_fromothertop_pdgid );
  tree->Branch("tXq_nmatched"                , &tXq_nmatched );
  tree->Branch("tXq_MultipleMatching"        , &tXq_MultipleMatching );
  tree->Branch("tXq_b1Index"                 , &tXq_b1Index );
  tree->Branch("tXq_b2Index"                 , &tXq_b2Index );
  tree->Branch("tXq_b3Index"                 , &tXq_b3Index );
  tree->Branch("tXq_qIndex"                  , &tXq_qIndex );
  tree->Branch("b1b2dR"                      , &b1b2dR );
    
  tree->Branch("b1b2_boosted_dR"                      , &b1b2_boosted_dR );
    
  tree->Branch("JJ_Boosted_dR"                      , &JJ_Boosted_dR ); 
  tree->Branch("JJ_boost_dR"                      , &JJ_boost_dR ); 
  tree->Branch("cos_JJ_Boosted_dR"                      , &cos_JJ_Boosted_dR);   
  tree->Branch("cos_JJ_boost_dR"                      , &cos_JJ_boost_dR); 
  
  tree->Branch("e_ae_Boosted_dR"                      , &e_ae_Boosted_dR ); 
  tree->Branch("e_ae_boost_dR"                      , &e_ae_boost_dR ); 
  tree->Branch("cos_e_ae_Boosted_dR"                      , &cos_e_ae_Boosted_dR);   
  tree->Branch("cos_e_ae_boost_dR"                      , &cos_e_ae_boost_dR); 
    
  tree->Branch("ae_e_Boosted_dR"                      , &ae_e_Boosted_dR ); 
  tree->Branch("ae_e_boost_dR"                      , &ae_e_boost_dR ); 
  tree->Branch("cos_ae_e_Boosted_dR"                      , &cos_ae_e_Boosted_dR);   
  tree->Branch("cos_ae_e_boost_dR"                      , &cos_ae_e_boost_dR); 
    
  tree->Branch("e_ae_Boosted_dPhi"                      , &e_ae_Boosted_dPhi ); 
  tree->Branch("e_ae_boost_dPhi"                      , &e_ae_boost_dPhi ); 
  tree->Branch("cos_e_ae_Boosted_dPhi"                      , &cos_e_ae_Boosted_dPhi);   
  tree->Branch("cos_e_ae_boost_dPhi"                      , &cos_e_ae_boost_dPhi); 
    
  tree->Branch("ae_e_Boosted_dPhi"                      , &ae_e_Boosted_dPhi ); 
  tree->Branch("ae_e_boost_dPhi"                      , &ae_e_boost_dPhi ); 
  tree->Branch("cos_ae_e_Boosted_dPhi"                      , &cos_ae_e_Boosted_dPhi);   
  tree->Branch("cos_ae_e_boost_dPhi"                      , &cos_ae_e_boost_dPhi);  
  
  tree->Branch("e_ae_Boosted_dEta"                      , &e_ae_Boosted_dEta ); 
  tree->Branch("e_ae_boost_dEta"                      , &e_ae_boost_dEta ); 
  tree->Branch("cos_e_ae_Boosted_dEta"                      , &cos_e_ae_Boosted_dEta);   
  tree->Branch("cos_e_ae_boost_dEta"                      , &cos_e_ae_boost_dEta); 
    
  tree->Branch("ae_e_Boosted_dEta"                      , &ae_e_Boosted_dEta ); 
  tree->Branch("ae_e_boost_dEta"                      , &ae_e_boost_dEta ); 
  tree->Branch("cos_ae_e_Boosted_dEta"                      , &cos_ae_e_Boosted_dEta);   
  tree->Branch("cos_ae_e_boost_dEta"                      , &cos_ae_e_boost_dEta);    
    
  tree->Branch("e_ae_dR"                      , &e_ae_dR ); 
  tree->Branch("ae_e_dR"                      , &ae_e_dR ); 
  tree->Branch("e_ae_dPhi"                      , &e_ae_dPhi ); 
  tree->Branch("ae_e_dPhi"                      , &ae_e_dPhi ); 
  tree->Branch("e_ae_dEta"                      , &e_ae_dEta ); 
  tree->Branch("ae_e_dEta"                      , &ae_e_dEta ); 
    
  tree->Branch("Phi"                      , &Phi);
  tree->Branch("Cos_Phi"                      , &Cos_Phi);
    
  tree->Branch("JJ_Boosted_dPhi"                      , &JJ_Boosted_dPhi ); 
  tree->Branch("JJ_boost_dPhi"                      , &JJ_boost_dPhi ); 
  tree->Branch("cos_JJ_Boosted_dPhi"                      , &cos_JJ_Boosted_dPhi);   
  tree->Branch("cos_JJ_boost_dPhi"                      , &cos_JJ_boost_dPhi); 
    
  tree->Branch("JJ_Boosted_dEta"                      , &JJ_Boosted_dEta ); 
  tree->Branch("JJ_boost_dEta"                      , &JJ_boost_dEta ); 
  tree->Branch("cos_JJ_Boosted_dEta"                      , &cos_JJ_Boosted_dEta);   
  tree->Branch("cos_JJ_boost_dEta"                      , &cos_JJ_boost_dEta); 
    
  tree->Branch("Phi_plane_angle"                      , &Phi_plane_angle);   
    
    
  tree->Branch("b1b3dR"                      , &b1b3dR );
  tree->Branch("b2b3dR"                      , &b2b3dR );
  tree->Branch("b1qSdR"                      , &b1qSdR );
  tree->Branch("b2qSdR"                      , &b1qSdR );
  tree->Branch("minDRbb"                     , &minDRbb );
  tree->Branch("dRdR_Jpsi_K"                     , &dR_Jpsi_K );  

  //bbHtt specific
  tree->Branch("n_parton_bbnotfromH_pt"      , &n_parton_bbnotfromH_pt                );
  tree->Branch("n_parton_bbnotfromH_eta"     , &n_parton_bbnotfromH_eta               );
  tree->Branch("n_parton_bbnotfromH_phi"     , &n_parton_bbnotfromH_phi               );
  tree->Branch("n_parton_bbnotfromH_m"       , &n_parton_bbnotfromH_m                 );
  tree->Branch("n_parton_bbnotfromH_n"       , &n_parton_bbnotfromH_n                 );
  tree->Branch("n_parton_bbnotfromH_deltaeta", &n_parton_bbnotfromH_deltaeta          );
  tree->Branch("n_parton_bbnotfromH_deltaphi", &n_parton_bbnotfromH_deltaphi          );
  tree->Branch("n_parton_bbnotfromH_deltar"  , &n_parton_bbnotfromH_deltar            );

  //Generic collection
  /*  
  tree->Branch("n_parton_fromH_pt"       , &n_parton_fromH_pt 	              );
  tree->Branch("n_parton_fromH_eta"      , &n_parton_fromH_eta	              );
  tree->Branch("n_parton_fromH_phi"      , &n_parton_fromH_phi	              );
  tree->Branch("n_parton_fromH_m"        , &n_parton_fromH_m  	              );
  */
  //Alvaro
  tree->Branch("n_parton_fromH_pt_p"       , &n_parton_fromH_pt_p 	              );
  tree->Branch("n_parton_fromH_eta_p"      , &n_parton_fromH_eta_p	              );
  tree->Branch("n_parton_fromH_phi_p"      , &n_parton_fromH_phi_p	              );
  tree->Branch("n_parton_fromH_m_p"        , &n_parton_fromH_m_p  	              );
   
  tree->Branch("n_parton_fromH_pt_a"       , &n_parton_fromH_pt_a 	              );
  tree->Branch("n_parton_fromH_eta_a"      , &n_parton_fromH_eta_a	              );
  tree->Branch("n_parton_fromH_phi_a"      , &n_parton_fromH_phi_a	              );
  tree->Branch("n_parton_fromH_m_a"        , &n_parton_fromH_m_a  	              );
  //
    
  tree->Branch("n_parton_fromH_n"        , &n_parton_fromH_n  	              );
  tree->Branch("n_parton_fromH_deltaeta" , &n_parton_fromH_deltaeta           );
  tree->Branch("n_parton_fromH_deltaphi" , &n_parton_fromH_deltaphi           );
  tree->Branch("n_parton_fromH_deltar"   , &n_parton_fromH_deltar             );
  tree->Branch("n_parton_H_decays"       , &n_parton_H_decays                 );

  //Zb collection 
  tree->Branch("n_Zb_Z_pt"       , &n_Zb_Z_pt 	              );
  tree->Branch("n_Zb_Z_eta"      , &n_Zb_Z_eta	              );
  tree->Branch("n_Zb_Z_phi"      , &n_Zb_Z_phi	              );
  tree->Branch("n_Zb_Z_m"        , &n_Zb_Z_m  	              );
  tree->Branch("n_Zb_wzjets_pt"  , &n_Zb_wzjets_pt            );
  tree->Branch("n_Zb_wzjets_eta" , &n_Zb_wzjets_eta           );
  tree->Branch("n_Zb_wzjets_phi" , &n_Zb_wzjets_phi           );
  tree->Branch("n_Zb_wzjets_m"   , &n_Zb_wzjets_m             );
  tree->Branch("n_Zb_wzjets_n"   , &n_Zb_wzjets_n             );

  //H+ angular variables
  tree->Branch("n_parton_topdecaychain_lepton_pt"    , &n_parton_topdecaychain_lepton_pt      );
  tree->Branch("n_parton_topdecaychain_lepton_eta"   , &n_parton_topdecaychain_lepton_eta     );
  tree->Branch("n_parton_topdecaychain_lepton_phi"   , &n_parton_topdecaychain_lepton_phi     );
  tree->Branch("n_parton_topdecaychain_lepton_m"     , &n_parton_topdecaychain_lepton_m       );
  tree->Branch("n_parton_topdecaychain_neutrino_pt"  , &n_parton_topdecaychain_neutrino_pt    );
  tree->Branch("n_parton_topdecaychain_neutrino_eta" , &n_parton_topdecaychain_neutrino_eta   );
  tree->Branch("n_parton_topdecaychain_neutrino_phi" , &n_parton_topdecaychain_neutrino_phi   );
  tree->Branch("n_parton_topdecaychain_neutrino_m"   , &n_parton_topdecaychain_neutrino_m     );

  tree->Branch("n_HT", &n_HT );

  tree->Branch("n_angana_lep_top_dphi"           , &n_angana_lep_top_dphi	     );
  tree->Branch("n_angana_lep_top_dr"             , &n_angana_lep_top_dr		     );
  tree->Branch("n_angana_nu_top_dphi"            , &n_angana_nu_top_dphi      	     );
  tree->Branch("n_angana_nu_top_dr"              , &n_angana_nu_top_dr		     );
  tree->Branch("n_angana_lep_top_dphi_topcm"     , &n_angana_lep_top_dphi_topcm	     );
  tree->Branch("n_angana_nu_top_dphi_topcm"      , &n_angana_nu_top_dphi_topcm	     );
  
  //B+ decay variables
  tree->Branch("n_B_pt"    , &n_B_pt      );
  tree->Branch("n_B_eta"   , &n_B_eta     );
  tree->Branch("n_B_phi"   , &n_B_phi     );
  tree->Branch("n_B_m"     , &n_B_m       );
    
  tree->Branch("n_J_pt"  , &n_J_pt    );
  tree->Branch("n_J_eta" , &n_J_eta   );
  tree->Branch("n_J_phi" , &n_J_phi   );
  tree->Branch("n_J_m"   , &n_J_m     );  
    
  tree->Branch("n_K_pt"  , &n_K_pt    );
  tree->Branch("n_K_eta" , &n_K_eta   );
  tree->Branch("n_K_phi" , &n_K_phi   );
  tree->Branch("n_K_m"   , &n_K_m     );
    
  tree->Branch("n_K_J_dphi", &n_K_J_dphi);
  tree->Branch("n_K_J_deta", &n_K_J_deta);  
    
  return;
}


void TruthAna_NtupleManager::FillParticleCollection(std::vector<TLorentzVector> vtlv, std::string coll_name)
{
  int size=vtlv.size();
  
  for(int i=0; i<size; i++)
    {
      float pt =vtlv[i].Pt();
      float eta=vtlv[i].Rapidity();
      float phi=vtlv[i].Phi();
      float m  =vtlv[i].M();


      if(coll_name=="Topdecaychain_lepton")
	{
	  n_parton_topdecaychain_lepton_pt .push_back( pt  );
	  n_parton_topdecaychain_lepton_eta.push_back( eta );
	  n_parton_topdecaychain_lepton_phi.push_back( phi );
	  n_parton_topdecaychain_lepton_m  .push_back( m   );
	}
      else if(coll_name=="Topdecaychain_neutrino")
	{
	  n_parton_topdecaychain_neutrino_pt .push_back( pt  );
	  n_parton_topdecaychain_neutrino_eta.push_back( eta );
	  n_parton_topdecaychain_neutrino_phi.push_back( phi );
	  n_parton_topdecaychain_neutrino_m  .push_back( m   );
	}
      else if(coll_name=="Neutrinos")
	{
	  n_neutrinos_pt  .push_back( pt  );
	  n_neutrinos_eta .push_back( eta );
	  n_neutrinos_phi .push_back( phi );
	  n_neutrinos_m   .push_back( m   ); 
	  n_neutrinos_n   =size;
	}

      else if(coll_name=="Electrons")
	{
	  n_electrons_pt  .push_back( pt  );
	  n_electrons_eta .push_back( eta );
	  n_electrons_phi .push_back( phi );
	  n_electrons_m   .push_back( m   ); 
	  n_electrons_n   =size;
	}

      else if(coll_name=="Muons")
	{
	  n_muons_pt  .push_back( pt  );
	  n_muons_eta .push_back( eta );
	  n_muons_phi .push_back( phi );
	  n_muons_m   .push_back( m   ); 
	  n_muons_n   =size;
	}

      else if(coll_name=="Taus")
	{
	  n_taus_pt  .push_back( pt  );
	  n_taus_eta .push_back( eta );
	  n_taus_phi .push_back( phi );
	  n_taus_m   .push_back( m   ); 
	  n_taus_n   =size;
	}

      else if(coll_name=="AntiKt4TruthWZJets")
	{
	  n_wzjets_pt  .push_back( pt  );
	  n_wzjets_eta .push_back( eta );
	  n_wzjets_phi .push_back( phi );
	  n_wzjets_m   .push_back( m   ); 
	  n_wzjets_n   =size;
	}

      else if(coll_name=="ParticleWBoson")
	{
	  n_particle_wboson_pt  .push_back( pt  );
	  n_particle_wboson_eta .push_back( eta );
	  n_particle_wboson_phi .push_back( phi );
	  n_particle_wboson_m   .push_back( m   ); 
	  n_particle_wboson_n   =size;
	}

      else if(coll_name=="UpQuarks")
      	{
	  n_parton_u_pt  .push_back( pt  );
	  n_parton_u_eta .push_back( eta );
	  n_parton_u_phi .push_back( phi );
	  n_parton_u_m   .push_back( m   ); 
	  n_parton_u_n   =size;
	}
      
      else if(coll_name=="Anti-UpQuarks")
      	{
	  n_parton_au_pt  .push_back( pt  );
	  n_parton_au_eta .push_back( eta );
	  n_parton_au_phi .push_back( phi );
	  n_parton_au_m   .push_back( m   ); 
	  n_parton_au_n   =size;
	}
      
      else if(coll_name=="DownQuarks")
      	{
	  n_parton_d_pt  .push_back( pt  );
	  n_parton_d_eta .push_back( eta );
	  n_parton_d_phi .push_back( phi );
	  n_parton_d_m   .push_back( m   ); 
	  n_parton_d_n   =size;
	}
      
      else if(coll_name=="Anti-DownQuarks")
      	{
	  n_parton_ad_pt  .push_back( pt  );
	  n_parton_ad_eta .push_back( eta );
	  n_parton_ad_phi .push_back( phi );
	  n_parton_ad_m   .push_back( m   ); 
	  n_parton_ad_n   =size;
	}

      else if(coll_name=="CharmQuarks")
      	{
	  n_parton_c_pt  .push_back( pt  );
	  n_parton_c_eta .push_back( eta );
	  n_parton_c_phi .push_back( phi );
	  n_parton_c_m   .push_back( m   ); 
	  n_parton_c_n   =size;
	}
      
      else if(coll_name=="Anti-CharmQuarks")
      	{
	  n_parton_ac_pt  .push_back( pt  );
	  n_parton_ac_eta .push_back( eta );
	  n_parton_ac_phi .push_back( phi );
	  n_parton_ac_m   .push_back( m   ); 
	  n_parton_ac_n   =size;
	}

      else if(coll_name=="StrangeQuarks")
      	{
	  n_parton_s_pt  .push_back( pt  );
	  n_parton_s_eta .push_back( eta );
	  n_parton_s_phi .push_back( phi );
	  n_parton_s_m   .push_back( m   ); 
	  n_parton_s_n   =size;
	}
      
      else if(coll_name=="Anti-StrangeQuarks")
      	{
	  n_parton_as_pt  .push_back( pt  );
	  n_parton_as_eta .push_back( eta );
	  n_parton_as_phi .push_back( phi );
	  n_parton_as_m   .push_back( m   ); 
	  n_parton_as_n   =size;
	}

      else if(coll_name=="TopQuarks")
      	{
	  n_parton_top_pt  .push_back( pt  );
	  n_parton_top_eta .push_back( eta );
	  n_parton_top_phi .push_back( phi );
	  n_parton_top_m   .push_back( m   ); 
	  n_parton_top_n   =size;
	}
      
      else if(coll_name=="Anti-TopQuarks")
      	{
	  n_parton_atop_pt  .push_back( pt  );
	  n_parton_atop_eta .push_back( eta );
	  n_parton_atop_phi .push_back( phi );
	  n_parton_atop_m   .push_back( m   ); 
	  n_parton_atop_n   =size;
	}

      else if(coll_name=="BottomQuarks")
      	{
	  n_parton_b_pt  .push_back( pt  );
	  n_parton_b_eta .push_back( eta );
	  n_parton_b_phi .push_back( phi );
	  n_parton_b_m   .push_back( m   ); 
	  n_parton_b_n   =size;
	}
      
      else if(coll_name=="Anti-BottomQuarks")
      	{
	  n_parton_ab_pt  .push_back( pt  );
	  n_parton_ab_eta .push_back( eta );
	  n_parton_ab_phi .push_back( phi );
	  n_parton_ab_m   .push_back( m   ); 
	  n_parton_ab_n   =size;
	}

      else if(coll_name=="WPlus")
      	{
	  n_parton_wplus_pt .push_back( pt  );
	  n_parton_wplus_eta.push_back( eta );
	  n_parton_wplus_phi.push_back( phi );
	  n_parton_wplus_m  .push_back( m   ); 
	  n_parton_wplus_n  =size;
	}
      
      else if(coll_name=="WMinus")
      	{
	  n_parton_wminus_pt .push_back( pt  );
	  n_parton_wminus_eta.push_back( eta );
	  n_parton_wminus_phi.push_back( phi );
	  n_parton_wminus_m  .push_back( m   ); 
	  n_parton_wminus_n  =size;
	}
      
      else if(coll_name=="Higgs")
	{
	  n_parton_higgs_pt  .push_back( pt  );
	  n_parton_higgs_eta .push_back( eta );
	  n_parton_higgs_phi .push_back( phi );
	  n_parton_higgs_m   .push_back( m   ); 
	  n_parton_higgs_n   =size;               
	}
      ////
      else if(coll_name=="B_plus")
	{
	  n_hadron_B_plus_pt  .push_back( pt  );
	  n_hadron_B_plus_eta .push_back( eta );
	  n_hadron_B_plus_phi .push_back( phi );
	  n_hadron_B_plus_m   .push_back( m   );          
	}
      else if(coll_name=="B_minus")
	{
	  n_hadron_B_minus_pt  .push_back( pt  );
	  n_hadron_B_minus_eta .push_back( eta );
	  n_hadron_B_minus_phi .push_back( phi );
	  n_hadron_B_minus_m   .push_back( m   );     
	}    
      
      else if(coll_name=="J_B")
	{
	  n_hadron_J_B_pt  .push_back( pt  );
	  n_hadron_J_B_eta .push_back( eta );
	  n_hadron_J_B_phi .push_back( phi );
	  n_hadron_J_B_m   .push_back( m   );          
	}

      else if(coll_name=="J_aB")
	{
	  n_hadron_J_aB_pt  .push_back( pt  );
	  n_hadron_J_aB_eta .push_back( eta );
	  n_hadron_J_aB_phi .push_back( phi );
	  n_hadron_J_aB_m   .push_back( m   );          
	}
      
      else if(coll_name=="K_B")
	{
	  n_hadron_K_B_pt  .push_back( pt  );
	  n_hadron_K_B_eta .push_back( eta );
	  n_hadron_K_B_phi .push_back( phi );
	  n_hadron_K_B_m   .push_back( m   );          
	}

      else if(coll_name=="K_aB")
	{
	  n_hadron_K_aB_pt  .push_back( pt  );
	  n_hadron_K_aB_eta .push_back( eta );
	  n_hadron_K_aB_phi .push_back( phi );
	  n_hadron_K_aB_m   .push_back( m   );          
	}
      
      else if(coll_name=="J_B_Boosted")
	{
	  n_hadron_J_B_Boosted_pt  .push_back( pt  );
	  n_hadron_J_B_Boosted_eta .push_back( eta );
	  n_hadron_J_B_Boosted_phi .push_back( phi );
	  n_hadron_J_B_Boosted_m   .push_back( m   );          
	}

      else if(coll_name=="J_aB_Boosted")
	{
	  n_hadron_J_aB_Boosted_pt  .push_back( pt  );
	  n_hadron_J_aB_Boosted_eta .push_back( eta );
	  n_hadron_J_aB_Boosted_phi .push_back( phi );
	  n_hadron_J_aB_Boosted_m   .push_back( m   );          
	}
      
      else if(coll_name=="J_B_boost")
	{
	  n_hadron_J_B_boost_pt  .push_back( pt  );
	  n_hadron_J_B_boost_eta .push_back( eta );
	  n_hadron_J_B_boost_phi .push_back( phi );
	  n_hadron_J_B_boost_m   .push_back( m   );          
	}

      else if(coll_name=="J_aB_boost")
	{
	  n_hadron_J_aB_boost_pt  .push_back( pt  );
	  n_hadron_J_aB_boost_eta .push_back( eta );
	  n_hadron_J_aB_boost_phi .push_back( phi );
	  n_hadron_J_aB_boost_m   .push_back( m   );          
	}
      

      else if(coll_name=="Zb_Z")
	{
	  n_Zb_Z_pt  .push_back( pt  );
	  n_Zb_Z_eta .push_back( eta );
	  n_Zb_Z_phi .push_back( phi );
	  n_Zb_Z_m   .push_back( m   );
	}
      else if(coll_name=="Zb_wzjets")
	{
	  n_Zb_wzjets_pt  .push_back( pt  );
	  n_Zb_wzjets_eta .push_back( eta );
	  n_Zb_wzjets_phi .push_back( phi );
	  n_Zb_wzjets_m   .push_back( m   );
	  n_Zb_wzjets_n   =size;               
	}
      else
	{
	  std::cout<<"Trying to fill the non-existing particle collection "<< coll_name<< " exit now!"<<std::endl;
	  exit(0);
	}
    }
  return;
}

void TruthAna_NtupleManager::FillHCollection(std::vector<TLorentzVector> vtlv, std::string coll_name)
{
  if(coll_name=="TTFromH")
    {
      n_parton_ttfromH_pt .push_back(vtlv[2].Pt());
      n_parton_ttfromH_eta.push_back(vtlv[2].Rapidity());
      n_parton_ttfromH_phi.push_back(vtlv[2].Phi());
      n_parton_ttfromH_m  .push_back(vtlv[2].M());
      n_parton_ttfromH_n  =vtlv.size();
      n_parton_ttfromH_deltaeta.push_back(vtlv[0].Rapidity()-vtlv[1].Rapidity());  
      n_parton_ttfromH_deltaphi.push_back(vtlv[0].DeltaPhi(vtlv[1]));  
      n_parton_ttfromH_deltar  .push_back(vtlv[0].DeltaR(vtlv[1]));  
    }

  else if(coll_name=="TTNotFromH")
    {
      n_parton_ttnotfromH_pt .push_back(vtlv[2].Pt());
      n_parton_ttnotfromH_eta.push_back(vtlv[2].Rapidity());
      n_parton_ttnotfromH_phi.push_back(vtlv[2].Phi());
      n_parton_ttnotfromH_m  .push_back(vtlv[2].M());
      n_parton_ttnotfromH_n  =vtlv.size();
      n_parton_ttnotfromH_deltaeta.push_back(vtlv[0].Rapidity()-vtlv[1].Rapidity());  
      n_parton_ttnotfromH_deltaphi.push_back(vtlv[0].DeltaPhi(vtlv[1]));  
      n_parton_ttnotfromH_deltar  .push_back(vtlv[0].DeltaR(vtlv[1]));  
    }

  else if(coll_name=="BBNotFromH")
    {
      
      n_parton_bbnotfromH_pt .push_back(vtlv[2].Pt());
      n_parton_bbnotfromH_eta.push_back(vtlv[2].Rapidity());
      n_parton_bbnotfromH_phi.push_back(vtlv[2].Phi());
      n_parton_bbnotfromH_m  .push_back(vtlv[2].M());
      
      
      n_parton_bbnotfromH_n  =vtlv.size();
      
      n_parton_bbnotfromH_deltaeta.push_back(vtlv[0].Rapidity()-vtlv[1].Rapidity());  
      n_parton_bbnotfromH_deltaphi.push_back(vtlv[0].DeltaPhi(vtlv[1]));  
      n_parton_bbnotfromH_deltar  .push_back(vtlv[0].DeltaR(vtlv[1]));  
    }

  else if(coll_name=="TCNotFromH")
    {
      n_parton_tcnotfromH_pt .push_back(vtlv[2].Pt());
      n_parton_tcnotfromH_eta.push_back(vtlv[2].Rapidity());
      n_parton_tcnotfromH_phi.push_back(vtlv[2].Phi());
      n_parton_tcnotfromH_m  .push_back(vtlv[2].M());
      n_parton_tcnotfromH_n  =vtlv.size();
      n_parton_tcnotfromH_deltaeta.push_back(vtlv[0].Rapidity()-vtlv[1].Rapidity());  
      n_parton_tcnotfromH_deltaphi.push_back(vtlv[0].DeltaPhi(vtlv[1]));  
      n_parton_tcnotfromH_deltar  .push_back(vtlv[0].DeltaR(vtlv[1]));  
    }

  else if(coll_name=="FromH")
    {
      /*
      n_parton_fromH_pt_ .push_back(vtlv[2].Pt());
      n_parton_fromH_eta.push_back(vtlv[2].Rapidity());
      n_parton_fromH_phi.push_back(vtlv[2].Phi());
      n_parton_fromH_m  .push_back(vtlv[2].M());
      */
      //ALvaro
      n_parton_fromH_pt_p .push_back(vtlv[0].Pt());
      n_parton_fromH_eta_p.push_back(vtlv[0].Rapidity());
      n_parton_fromH_phi_p.push_back(vtlv[0].Phi());
      n_parton_fromH_m_p  .push_back(vtlv[0].M());
      
      n_parton_fromH_pt_a .push_back(vtlv[1].Pt());
      n_parton_fromH_eta_a.push_back(vtlv[1].Rapidity());
      n_parton_fromH_phi_a.push_back(vtlv[1].Phi());
      n_parton_fromH_m_a  .push_back(vtlv[1].M());
      //
      
      n_parton_fromH_n  =vtlv.size();
      //n_parton_fromH_n  =vtlv[0].pdgid;
      n_parton_fromH_deltaeta.push_back(vtlv[0].Rapidity()-vtlv[1].Rapidity());  
      n_parton_fromH_deltaphi.push_back(vtlv[0].DeltaPhi(vtlv[1]));  
      n_parton_fromH_deltar  .push_back(vtlv[0].DeltaR(vtlv[1]));  
    }
  else if(coll_name=="FromB")
   {
      n_B_pt.push_back(vtlv[0].Pt());
      n_B_eta.push_back(vtlv[0].Rapidity());
      n_B_phi.push_back(vtlv[0].Phi());
      n_B_m.push_back(vtlv[0].M());
      
      n_J_pt.push_back(vtlv[1].Pt());
      n_J_eta.push_back(vtlv[1].Rapidity());
      n_J_phi.push_back(vtlv[1].Phi());
      n_J_m.push_back(vtlv[1].M());
      
      n_K_pt.push_back(vtlv[2].Pt());
      n_K_eta.push_back(vtlv[2].Rapidity());
      n_K_phi.push_back(vtlv[2].Phi());
      n_K_m.push_back(vtlv[2].M());
      
      n_K_J_dphi.push_back(vtlv[1].DeltaPhi(vtlv[2]));
      n_K_J_deta.push_back(vtlv[1].Rapidity()-vtlv[2].Rapidity());      
      
      
      
  }
        
  else 	
    {
      std::cout<<"Trying to fill the non-existing particle collection "<< coll_name<< " exit now!"<<std::endl;
      exit(0);
    }
  
  return;
}

void TruthAna_NtupleManager::FillHCollection(std::vector<myParticle> vtlv, std::string coll_name)
{
  if(coll_name=="TCNotFromH")
    {
      n_parton_tcnotfromH_pt .push_back(vtlv[2].FourVector.Pt());
      n_parton_tcnotfromH_eta.push_back(vtlv[2].FourVector.Rapidity());
      n_parton_tcnotfromH_phi.push_back(vtlv[2].FourVector.Phi());
      n_parton_tcnotfromH_m  .push_back(vtlv[2].FourVector.M());
      n_parton_tcnotfromH_n  =vtlv.size();
      n_parton_tcnotfromH_deltaeta.push_back(vtlv[0].FourVector.Rapidity()-vtlv[1].FourVector.Rapidity());  
      n_parton_tcnotfromH_deltaphi.push_back(vtlv[0].FourVector.DeltaPhi(vtlv[1].FourVector));  
      n_parton_tcnotfromH_deltar  .push_back(vtlv[0].FourVector.DeltaR(vtlv[1].FourVector));  
    }
  else 	
    {
      std::cout<<"Trying to fill the non-existing particle collection "<< coll_name<< " exit now!"<<std::endl;
      exit(0);
    }
  
  return;
}

void TruthAna_NtupleManager::FillHCollection(myParticle vtlv, std::string coll_name)
{
if(coll_name=="SisH")
    {
      n_parton_sisterH_pt    = vtlv.FourVector.Pt();
      n_parton_sisterH_eta   = vtlv.FourVector.Rapidity();
      n_parton_sisterH_phi   = vtlv.FourVector.Phi();
      n_parton_sisterH_m     = vtlv.FourVector.M();
      n_parton_sisterH_pdgid = vtlv.pdgid;
    }
  else if(coll_name=="FromOtherTop")
    {
      n_parton_fromothertop_pt    = vtlv.FourVector.Pt();
      n_parton_fromothertop_eta   = vtlv.FourVector.Rapidity();
      n_parton_fromothertop_phi   = vtlv.FourVector.Phi();
      n_parton_fromothertop_m     = vtlv.FourVector.M();
      n_parton_fromothertop_pdgid = vtlv.pdgid;
    }
  else 	
    {
      std::cout<<"Trying to fill the non-existing particle collection "<< coll_name<< " exit now!"<<std::endl;
      exit(0);
    }
  return;
}


//Temporary implementation
void TruthAna_NtupleManager::FillEventCollection(float mc_event_weight)
{
  n_mc_event_weight=mc_event_weight;
  return;
}

void TruthAna_NtupleManager::FillScalarCollection(float variable, std::string coll_name)
{
  if     (coll_name=="Angana_lep_top_dphi"      ){ n_angana_lep_top_dphi      =variable;}
  else if(coll_name=="Angana_lep_top_dr"        ){ n_angana_lep_top_dr        =variable;}
  else if(coll_name=="Angana_nu_top_dphi"       ){ n_angana_nu_top_dphi       =variable;}
  else if(coll_name=="Angana_nu_top_dr"         ){ n_angana_nu_top_dr         =variable;}
  else if(coll_name=="Angana_lep_top_dphi_topcm"){ n_angana_lep_top_dphi_topcm=variable;}
  else if(coll_name=="Angana_nu_top_dphi_topcm" ){ n_angana_nu_top_dphi_topcm =variable;}
  else if(coll_name=="HT" ){ n_HT =variable;}
  else if(coll_name=="b1b2dR"                   ){ b1b2dR                     =variable;}
  else if(coll_name=="b1b2_boosted_dR"                   ){ b1b2_boosted_dR                     =variable;}  
    
  else if(coll_name=="JJ_Boosted_dR"                   ){ JJ_Boosted_dR                     =variable;}  
  else if(coll_name=="JJ_boost_dR"                   ){ JJ_boost_dR                     =variable;}  
  else if(coll_name=="cos_JJ_Boosted_dR"                   ){ cos_JJ_Boosted_dR                     =variable;}  
  else if(coll_name=="cos_JJ_boost_dR"                   ){ cos_JJ_boost_dR                     =variable;}   
  
  else if(coll_name=="e_ae_Boosted_dR"                   ){ e_ae_Boosted_dR                     =variable;}  
  else if(coll_name=="e_ae_boost_dR"                   ){ e_ae_boost_dR                     =variable;}  
  else if(coll_name=="cos_e_ae_Boosted_dR"                   ){ cos_e_ae_Boosted_dR                     =variable;}  
  else if(coll_name=="cos_e_ae_boost_dR"                   ){ cos_e_ae_boost_dR                     =variable;}     

  else if(coll_name=="ae_e_Boosted_dR"                   ){ ae_e_Boosted_dR                     =variable;}  
  else if(coll_name=="ae_e_boost_dR"                   ){ ae_e_boost_dR                     =variable;}  
  else if(coll_name=="cos_ae_e_Boosted_dR"                   ){ cos_ae_e_Boosted_dR                     =variable;}  
  else if(coll_name=="cos_ae_e_boost_dR"                   ){ cos_ae_e_boost_dR                     =variable;}  
    
  else if(coll_name=="e_ae_Boosted_dPhi"                   ){ e_ae_Boosted_dPhi                     =variable;}  
  else if(coll_name=="e_ae_boost_dPhi"                   ){ e_ae_boost_dPhi                     =variable;}  
  else if(coll_name=="cos_e_ae_Boosted_dPhi"                   ){ cos_e_ae_Boosted_dPhi                     =variable;}  
  else if(coll_name=="cos_e_ae_boost_dPhi"                   ){ cos_e_ae_boost_dPhi                     =variable;}     

  else if(coll_name=="ae_e_Boosted_dPhi"                   ){ ae_e_Boosted_dPhi                     =variable;}  
  else if(coll_name=="ae_e_boost_dPhi"                   ){ ae_e_boost_dPhi                     =variable;}  
  else if(coll_name=="cos_ae_e_Boosted_dPhi"                   ){ cos_ae_e_Boosted_dPhi                     =variable;}  
  else if(coll_name=="cos_ae_e_boost_dPhi"                   ){ cos_ae_e_boost_dPhi                     =variable;} 
    
  else if(coll_name=="e_ae_Boosted_dEta"                   ){ e_ae_Boosted_dEta                     =variable;}  
  else if(coll_name=="e_ae_boost_dEta"                   ){ e_ae_boost_dEta                     =variable;}  
  else if(coll_name=="cos_e_ae_Boosted_dEta"                   ){ cos_e_ae_Boosted_dEta                     =variable;}  
  else if(coll_name=="cos_e_ae_boost_dEta"                   ){ cos_e_ae_boost_dEta                     =variable;}     

  else if(coll_name=="ae_e_Boosted_dEta"                   ){ ae_e_Boosted_dEta                     =variable;}  
  else if(coll_name=="ae_e_boost_dEta"                   ){ ae_e_boost_dEta                     =variable;}  
  else if(coll_name=="cos_ae_e_Boosted_dEta"                   ){ cos_ae_e_Boosted_dEta                     =variable;}  
  else if(coll_name=="cos_ae_e_boost_dEta"                   ){ cos_ae_e_boost_dEta                     =variable;}   
    
  else if(coll_name=="e_ae_dR"                       ){ e_ae_dR                         =variable;}  
  else if(coll_name=="ae_e_dR"                       ){ ae_e_dR                         =variable;}  
  else if(coll_name=="e_ae_dPhi"                     ){ e_ae_dPhi                       =variable;}  
  else if(coll_name=="ae_e_dPhi"                     ){ ae_e_dPhi                       =variable;}  
  else if(coll_name=="e_ae_dEta"                     ){ e_ae_dEta                       =variable;}  
  else if(coll_name=="ae_e_dEta"                     ){ ae_e_dEta                       =variable;}  
    
  else if(coll_name=="Phi"                   ){ Phi                     =variable;}
  else if(coll_name=="Cos_Phi"                   ){ Cos_Phi                     =variable;}  
    
  else if(coll_name=="JJ_Boosted_dPhi"                   ){ JJ_Boosted_dPhi                     =variable;}  
  else if(coll_name=="JJ_boost_dPhi"                   ){ JJ_boost_dPhi                     =variable;}  
  else if(coll_name=="cos_JJ_Boosted_dPhi"                   ){ cos_JJ_Boosted_dPhi                     =variable;}  
  else if(coll_name=="cos_JJ_boost_dPhi"                   ){ cos_JJ_boost_dPhi                     =variable;} 

  else if(coll_name=="JJ_Boosted_dEta"                   ){ JJ_Boosted_dEta                     =variable;}  
  else if(coll_name=="JJ_boost_dEta"                   ){ JJ_boost_dEta                     =variable;}  
  else if(coll_name=="cos_JJ_Boosted_dEta"                   ){ cos_JJ_Boosted_dEta                     =variable;}  
  else if(coll_name=="cos_JJ_boost_dEta"                   ){ cos_JJ_boost_dEta                     =variable;}     
  else if(coll_name=="Phi_plane_angle"                   ){ Phi_plane_angle                     =variable;}  
    
  else if(coll_name=="b1b3dR"                   ){ b1b3dR                     =variable;}
  else if(coll_name=="b2b3dR"                   ){ b2b3dR                     =variable;}
  else if(coll_name=="b1qSdR"                   ){ b1qSdR                     =variable;}
  else if(coll_name=="b2qSdR"                   ){ b2qSdR                     =variable;}  
  else if(coll_name=="minDRbb"                  ){ minDRbb                    =variable;}
  else if(coll_name=="dR_Jpsi_K"                  ){ dR_Jpsi_K                    =variable;}  
  else
    {
      std::cout<<"Trying to fill the non-existing particle collection "<< coll_name<< " exit now!"<<std::endl;
      exit(0);
    }
  return;
}

void TruthAna_NtupleManager::FillScalarCollection(int variable, std::string coll_name)
{
  if     (coll_name=="tXq_nmatched"           ) {tXq_nmatched         = variable;}
  else if(coll_name=="tXq_MultipleMatching"   ) {tXq_MultipleMatching = variable;}
  else if(coll_name=="tXq_b1Index"            ) {tXq_b1Index          = variable;}
  else if(coll_name=="tXq_b2Index"            ) {tXq_b2Index          = variable;}
  else if(coll_name=="tXq_b3Index"            ) {tXq_b3Index          = variable;}
  else if(coll_name=="tXq_qIndex"             ) {tXq_qIndex           = variable;}
  else
    {
      std::cout<<"Int Collection: Trying to fill the non-existing particle collection "<< coll_name<< " exit now!"<<std::endl;
      exit(0);
    }
  return;
}

void TruthAna_NtupleManager::FillEventCollectionHack(std::vector<int> parton_H_decays)
{
  int size=parton_H_decays.size();
  
  for(int i=0; i<size; i++)
    {
      n_parton_H_decays.push_back(parton_H_decays[i]);
    }
}

void TruthAna_NtupleManager::ClearCollections()
{
  n_mc_event_weight=0.;
  
  n_neutrinos_pt .clear();                                                                           
  n_neutrinos_eta.clear();                                                                           
  n_neutrinos_phi.clear();                                                                           
  n_neutrinos_m  .clear();                                                                           
  n_neutrinos_n  =0;

  n_muons_pt .clear();                                                                               
  n_muons_eta.clear();                                                                               
  n_muons_phi.clear();                                                                               
  n_muons_m  .clear();                                                                               
  n_muons_n  =0;                                                                               

  n_electrons_pt .clear();                                                                           
  n_electrons_eta.clear();                                                                           
  n_electrons_phi.clear();                                                                           
  n_electrons_m  .clear();                                                                           
  n_electrons_n  =0;

  n_taus_pt .clear();                                                                                
  n_taus_eta.clear();                                                                                
  n_taus_phi.clear();                                                                                
  n_taus_m  .clear();           
  n_taus_n  =0;
  
  n_wzjets_pt .clear();                                                                              
  n_wzjets_eta.clear();                                                                              
  n_wzjets_phi.clear();                                                                              
  n_wzjets_m  .clear();         
  n_wzjets_n  =0;      

  n_particle_wboson_pt .clear(); 
  n_particle_wboson_eta.clear(); 
  n_particle_wboson_phi.clear(); 
  n_particle_wboson_m  .clear(); 
  n_particle_wboson_n  =0;       
  
  //Other collections
  
  //Generic parton collections

  n_parton_decay    =0;

  n_parton_higgs_pt .clear();                                                                               
  n_parton_higgs_eta.clear();                                                                               
  n_parton_higgs_phi.clear();                                                                               
  n_parton_higgs_m  .clear();                    
  n_parton_higgs_n  =0;

///
  n_hadron_B_plus_pt.clear();
  n_hadron_B_plus_eta.clear();
  n_hadron_B_plus_phi.clear();  
  n_hadron_B_plus_m.clear();  

  n_hadron_B_minus_pt.clear();
  n_hadron_B_minus_eta.clear();
  n_hadron_B_minus_phi.clear();  
  n_hadron_B_minus_m.clear();      

  n_hadron_J_B_pt.clear();
  n_hadron_J_B_eta.clear();
  n_hadron_J_B_phi.clear();  
  n_hadron_J_B_m.clear();    
    
  n_hadron_J_aB_pt.clear();
  n_hadron_J_aB_eta.clear();
  n_hadron_J_aB_phi.clear();  
  n_hadron_J_aB_m.clear();    
    
  n_hadron_K_B_pt.clear();
  n_hadron_K_B_eta.clear();
  n_hadron_K_B_phi.clear();  
  n_hadron_K_B_m.clear();    
    
  n_hadron_K_aB_pt.clear();
  n_hadron_K_aB_eta.clear();
  n_hadron_K_aB_phi.clear();  
  n_hadron_K_aB_m.clear(); 
    
  n_hadron_J_B_Boosted_pt.clear();
  n_hadron_J_B_Boosted_eta.clear();
  n_hadron_J_B_Boosted_phi.clear();  
  n_hadron_J_B_Boosted_m.clear();    
    
  n_hadron_J_aB_Boosted_pt.clear();
  n_hadron_J_aB_Boosted_eta.clear();
  n_hadron_J_aB_Boosted_phi.clear();  
  n_hadron_J_aB_Boosted_m.clear();
    
  n_hadron_J_B_boost_pt.clear();
  n_hadron_J_B_boost_eta.clear();
  n_hadron_J_B_boost_phi.clear();  
  n_hadron_J_B_boost_m.clear();    
    
  n_hadron_J_aB_boost_pt.clear();
  n_hadron_J_aB_boost_eta.clear();
  n_hadron_J_aB_boost_phi.clear();  
  n_hadron_J_aB_boost_m.clear();    
    
     

  n_parton_u_pt .clear();                                                                               
  n_parton_u_eta.clear();                                                                               
  n_parton_u_phi.clear();                                                                               
  n_parton_u_m  .clear();                    
  n_parton_u_n  =0;

  n_parton_au_pt .clear();                                                                               
  n_parton_au_eta.clear();                                                                               
  n_parton_au_phi.clear();                                                                               
  n_parton_au_m  .clear();                    
  n_parton_au_n =0;

  n_parton_d_pt .clear();                                                                               
  n_parton_d_eta.clear();                                                                               
  n_parton_d_phi.clear();                                                                               
  n_parton_d_m  .clear();                    
  n_parton_d_n  =0;

  n_parton_ad_pt .clear();                                                                               
  n_parton_ad_eta.clear();                                                                               
  n_parton_ad_phi.clear();                                                                               
  n_parton_ad_m  .clear();                    
  n_parton_ad_n =0;

  n_parton_c_pt .clear();                                                                               
  n_parton_c_eta.clear();                                                                               
  n_parton_c_phi.clear();                                                                               
  n_parton_c_m  .clear();                    
  n_parton_c_n  =0;

  n_parton_ac_pt .clear();                                                                               
  n_parton_ac_eta.clear();                                                                               
  n_parton_ac_phi.clear();                                                                               
  n_parton_ac_m  .clear();                    
  n_parton_ac_n =0;

  n_parton_s_pt .clear();                                                                               
  n_parton_s_eta.clear();                                                                               
  n_parton_s_phi.clear();                                                                               
  n_parton_s_m  .clear();                    
  n_parton_s_n  =0;

  n_parton_as_pt .clear();                                                                               
  n_parton_as_eta.clear();                                                                               
  n_parton_as_phi.clear();                                                                               
  n_parton_as_m  .clear();                    
  n_parton_as_n =0;

  n_parton_top_pt .clear();                                                                               
  n_parton_top_eta.clear();                                                                               
  n_parton_top_phi.clear();                                                                               
  n_parton_top_m  .clear();                    
  n_parton_top_n  =0;

  n_parton_atop_pt .clear();                                                                               
  n_parton_atop_eta.clear();                                                                               
  n_parton_atop_phi.clear();                                                                               
  n_parton_atop_m  .clear();                    
  n_parton_atop_n   =0;      

  n_parton_b_pt .clear();                                                                               
  n_parton_b_eta.clear();                                                                               
  n_parton_b_phi.clear();                                                                               
  n_parton_b_m  .clear();                    
  n_parton_b_n  =0;

  n_parton_ab_pt .clear();                                                                               
  n_parton_ab_eta.clear();                                                                               
  n_parton_ab_phi.clear();                                                                               
  n_parton_ab_m  .clear();                    
  n_parton_ab_n =0;

  n_parton_wplus_pt  .clear(); 
  n_parton_wplus_eta .clear();
  n_parton_wplus_phi .clear();
  n_parton_wplus_m   .clear();
  n_parton_wplus_n   =0;      

  n_parton_wminus_pt .clear();
  n_parton_wminus_eta.clear();
  n_parton_wminus_phi.clear();
  n_parton_wminus_m  .clear();
  n_parton_wminus_n  =0;      
  
  //bbHtt and ttHtt collections
  n_parton_ttfromH_pt .clear();                                                                           
  n_parton_ttfromH_eta.clear();                                                                             
  n_parton_ttfromH_phi.clear();                                                                             
  n_parton_ttfromH_m  .clear();                    
  n_parton_ttfromH_n  =0;
  n_parton_ttfromH_poscosphi .clear();                    
  n_parton_ttfromH_deltaeta  .clear();                    
  n_parton_ttfromH_deltaphi  .clear();                    
  n_parton_ttfromH_deltar    .clear();                    

  //ttHtt specific
  n_parton_ttnotfromH_pt .clear();                                                                           
  n_parton_ttnotfromH_eta.clear();                                                                         
  n_parton_ttnotfromH_phi.clear();                                                                         
  n_parton_ttnotfromH_m  .clear();                    
  n_parton_ttnotfromH_n  =0;
  n_parton_ttnotfromH_deltaeta  .clear();                    
  n_parton_ttnotfromH_deltaphi  .clear();                    
  n_parton_ttnotfromH_deltar    .clear();                    

  //bbHtt specific
  n_parton_bbnotfromH_pt .clear();                                                                           
  n_parton_bbnotfromH_eta.clear();                                                                         
  n_parton_bbnotfromH_phi.clear();                                                                         
  n_parton_bbnotfromH_m  .clear();                    
  n_parton_bbnotfromH_n  =0;       
  n_parton_bbnotfromH_deltaeta.clear();                    
  n_parton_bbnotfromH_deltaphi.clear();                    
  n_parton_bbnotfromH_deltar  .clear();                    

  //tcHbb specific
  n_parton_tcnotfromH_pt .clear();                                                                          
  n_parton_tcnotfromH_eta.clear();                                                                        
  n_parton_tcnotfromH_phi.clear();                                                                        
  n_parton_tcnotfromH_m  .clear();                   
  n_parton_tcnotfromH_n  =0;          
  n_parton_tcnotfromH_deltaeta.clear();
  n_parton_tcnotfromH_deltaphi.clear();
  n_parton_tcnotfromH_deltar  .clear();
  n_parton_H_decays.clear();

  //tXq specific
  n_parton_sisterH_pt = -99;
  n_parton_sisterH_eta = -99;
  n_parton_sisterH_phi = -99;
  n_parton_sisterH_m = -99;
  n_parton_sisterH_pdgid = -99;
  n_parton_fromothertop_pt = -99;
  n_parton_fromothertop_eta = -99;
  n_parton_fromothertop_phi = -99;
  n_parton_fromothertop_m = -99;
  n_parton_fromothertop_pdgid = -99;
  tXq_nmatched = -1;
  tXq_MultipleMatching = -1;
  tXq_b1Index = -2;
  tXq_b2Index = -2;
  tXq_b3Index = -2;
  tXq_qIndex = -2;
  b1b2dR = 100;
  b1b2_boosted_dR = 100;
    
  JJ_Boosted_dR = 100; 
  JJ_boost_dR = 100;  
  cos_JJ_Boosted_dR = 100;
  cos_JJ_boost_dR = 100; 
    
  e_ae_Boosted_dR = 100; 
  e_ae_boost_dR = 100;  
  cos_e_ae_Boosted_dR = 100;
  cos_e_ae_boost_dR = 100;  
    
  ae_e_Boosted_dR = 100; 
  ae_e_boost_dR = 100;  
  cos_ae_e_Boosted_dR = 100;
  cos_ae_e_boost_dR = 100;
    
  e_ae_Boosted_dPhi = 100; 
  e_ae_boost_dPhi = 100;  
  cos_e_ae_Boosted_dPhi = 100;
  cos_e_ae_boost_dPhi = 100;  
    
  ae_e_Boosted_dPhi = 100; 
  ae_e_boost_dPhi = 100;  
  cos_ae_e_Boosted_dPhi = 100;
  cos_ae_e_boost_dPhi = 100;
  
  e_ae_Boosted_dEta = 100; 
  e_ae_boost_dEta = 100;  
  cos_e_ae_Boosted_dEta = 100;
  cos_e_ae_boost_dEta = 100;  
    
  ae_e_Boosted_dEta = 100; 
  ae_e_boost_dEta = 100;  
  cos_ae_e_Boosted_dEta = 100;
  cos_ae_e_boost_dEta = 100; 
    
  e_ae_dR  = 100;   
  ae_e_dR  = 100;   
  e_ae_dPhi  = 100; 
  ae_e_dPhi  = 100; 
  e_ae_dEta  = 100; 
  ae_e_dEta  = 100; 
    
  Phi = 100;
  Cos_Phi = 100;  
    
  JJ_Boosted_dPhi = 100; 
  JJ_boost_dPhi = 100;  
  cos_JJ_Boosted_dPhi = 100;
  cos_JJ_boost_dPhi = 100;  
    
  JJ_Boosted_dEta = 100; 
  JJ_boost_dEta = 100;  
  cos_JJ_Boosted_dEta = 100;
  cos_JJ_boost_dEta = 100;  
  
  Phi_plane_angle = 100;    
    
  b1b3dR = 100;
  b2b3dR = 100;
  b1qSdR = 100;
  b2qSdR = 100;  
  minDRbb = 100;
  dR_Jpsi_K = 100; 
  /*  
  n_parton_fromH_pt .clear();
  n_parton_fromH_eta.clear();
  n_parton_fromH_phi.clear();
  n_parton_fromH_m  .clear();
  */
  //Alvaro  
  n_parton_fromH_pt_p .clear();
  n_parton_fromH_eta_p.clear();
  n_parton_fromH_phi_p.clear();
  n_parton_fromH_m_p  .clear();
   
  n_parton_fromH_pt_a .clear();
  n_parton_fromH_eta_a.clear();
  n_parton_fromH_phi_a.clear();
  n_parton_fromH_m_a  .clear();
  //  
    
  n_parton_fromH_n  =0;      
  n_parton_fromH_deltaeta.clear(); 
  n_parton_fromH_deltaphi.clear(); 
  n_parton_fromH_deltar  .clear();

  //Zb collections
  n_Zb_Z_pt  .clear();
  n_Zb_Z_eta .clear();
  n_Zb_Z_phi .clear();
  n_Zb_Z_m   .clear();
  n_Zb_wzjets_pt .clear();
  n_Zb_wzjets_eta.clear();
  n_Zb_wzjets_phi.clear();
  n_Zb_wzjets_m  .clear();
  n_Zb_wzjets_n  =0;

  n_parton_topdecaychain_lepton_pt .clear();
  n_parton_topdecaychain_lepton_eta.clear();
  n_parton_topdecaychain_lepton_phi.clear();
  n_parton_topdecaychain_lepton_m  .clear();

  n_parton_topdecaychain_neutrino_pt .clear();
  n_parton_topdecaychain_neutrino_eta.clear();
  n_parton_topdecaychain_neutrino_phi.clear();
  n_parton_topdecaychain_neutrino_m  .clear();

  n_HT =0.;

  n_angana_lep_top_dphi =0.;
  n_angana_lep_top_dr   =0.;
  n_angana_nu_top_dphi  =0.;
  n_angana_nu_top_dr    =0.;

  n_angana_lep_top_dphi_topcm =0.;
  n_angana_nu_top_dphi_topcm  =0.;
    
  //B+
  n_B_pt.clear();
  n_B_eta.clear();
  n_B_phi.clear();
  n_B_m.clear();
    
  n_J_pt.clear();
  n_J_eta.clear();
  n_J_phi.clear();
  n_J_m.clear();
    
  n_K_pt.clear();
  n_K_eta.clear();
  n_K_phi.clear();
  n_K_m.clear(); 
    
  n_K_J_dphi.clear();
  n_K_J_deta.clear();
    
    
  return;
}
