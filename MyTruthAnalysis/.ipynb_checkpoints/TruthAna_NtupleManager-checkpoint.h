#ifndef MYTRUTHANALYSIS_TRUTHANA_NTUPLEMANAGER_H
#define MYTRUTHANALYSIS_TRUTHANA_NTUPLEMANAGER_H

#include <vector>
#include <TTree.h>

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "MyTruthAnalysis/MyClasses.h"


class TruthAna_OptionParser;

class TruthAna_NtupleManager
{
 public: 
  TruthAna_NtupleManager();
  TruthAna_NtupleManager(TruthAna_OptionParser *opt, TTree *tree);
  ~TruthAna_NtupleManager();

  float n_mc_event_weight ; //!                                                                          

  //TRUTH1 specific collections
  std::vector<float> n_neutrinos_pt ; //!                                                                          
  std::vector<float> n_neutrinos_eta; //!                                                                          
  std::vector<float> n_neutrinos_phi; //!                                                                          
  std::vector<float> n_neutrinos_m  ; //!                                                                          
  int n_neutrinos_n  ; //!                                                                          
  std::vector<float> n_muons_pt ; //!                                                                              
  std::vector<float> n_muons_eta; //!                                                                              
  std::vector<float> n_muons_phi; //!                                                                              
  std::vector<float> n_muons_m  ; //!                                                                              
  int n_muons_n  ; //!                                                                              
  std::vector<float> n_electrons_pt ; //!                                                                          
  std::vector<float> n_electrons_eta; //!                                                                          
  std::vector<float> n_electrons_phi; //!                                                                          
  std::vector<float> n_electrons_m  ; //!                                                                          
  int n_electrons_n  ; //!                                                                          
  std::vector<float> n_taus_pt ; //!                                                                               
  std::vector<float> n_taus_eta; //!                                                                               
  std::vector<float> n_taus_phi; //!                                                                               
  std::vector<float> n_taus_m  ; //!          
  int n_taus_n  ; //!          
  std::vector<float> n_wzjets_pt ; //!                                                                             
  std::vector<float> n_wzjets_eta; //!                                                                             
  std::vector<float> n_wzjets_phi; //!                                                                             
  std::vector<float> n_wzjets_m  ; //!        
  int n_wzjets_n  ; //!        
  std::vector<float> n_particle_wboson_pt ; //!                                                                          
  std::vector<float> n_particle_wboson_eta; //!                                                                            
  std::vector<float> n_particle_wboson_phi; //!                                                                           
  std::vector<float> n_particle_wboson_m  ; //!        
  int n_particle_wboson_n  ; //!        


  //Other collections
  
  std::vector<float> n_hadron_B_plus_pt ; //!                                                                         
  std::vector<float> n_hadron_B_plus_eta; //!                                                                         
  std::vector<float> n_hadron_B_plus_phi; //!                                                                         
  std::vector<float> n_hadron_B_plus_m  ; //!  
  
  std::vector<float> n_hadron_B_minus_pt ; //!                                                                         
  std::vector<float> n_hadron_B_minus_eta; //!                                                                         
  std::vector<float> n_hadron_B_minus_phi; //!                                                                         
  std::vector<float> n_hadron_B_minus_m  ; //!  
  
  std::vector<float> n_hadron_J_B_pt ; //!                                                                       
  std::vector<float> n_hadron_J_B_eta; //!                                                                         
  std::vector<float> n_hadron_J_B_phi; //!                                                                         
  std::vector<float> n_hadron_J_B_m  ; //!  
  
  std::vector<float> n_hadron_K_aB_pt ; //!                                                                       
  std::vector<float> n_hadron_K_aB_eta; //!                                                                         
  std::vector<float> n_hadron_K_aB_phi; //!                                                                         
  std::vector<float> n_hadron_K_aB_m  ; //!  
  
  std::vector<float> n_hadron_K_B_pt ; //!                                                                       
  std::vector<float> n_hadron_K_B_eta; //!                                                                         
  std::vector<float> n_hadron_K_B_phi; //!                                                                         
  std::vector<float> n_hadron_K_B_m  ; //!  
  
  std::vector<float> n_hadron_J_aB_pt ; //!                                                                       
  std::vector<float> n_hadron_J_aB_eta; //!                                                                         
  std::vector<float> n_hadron_J_aB_phi; //!                                                                         
  std::vector<float> n_hadron_J_aB_m  ; //!  
  
  std::vector<float> n_hadron_J_B_Boosted_pt ; //!                                                                     
  std::vector<float> n_hadron_J_B_Boosted_eta; //!                                                                     
  std::vector<float> n_hadron_J_B_Boosted_phi; //!                                                                     
  std::vector<float> n_hadron_J_B_Boosted_m  ; //!  
  
  std::vector<float> n_hadron_J_aB_Boosted_pt ; //!                                                                   
  std::vector<float> n_hadron_J_aB_Boosted_eta; //!                                                                   
  std::vector<float> n_hadron_J_aB_Boosted_phi; //!                                                                   
  std::vector<float> n_hadron_J_aB_Boosted_m  ; //!  
  
  std::vector<float> n_hadron_J_B_boost_pt ; //!                                                                     
  std::vector<float> n_hadron_J_B_boost_eta; //!                                                                     
  std::vector<float> n_hadron_J_B_boost_phi; //!                                                                     
  std::vector<float> n_hadron_J_B_boost_m  ; //!  
  
  std::vector<float> n_hadron_J_aB_boost_pt ; //!                                                                   
  std::vector<float> n_hadron_J_aB_boost_eta; //!                                                                   
  std::vector<float> n_hadron_J_aB_boost_phi; //!                                                                   
  std::vector<float> n_hadron_J_aB_boost_m  ; //!  
  
  //Generic parton collections
  int n_parton_decay; //!                                                                              
  std::vector<float> n_parton_higgs_pt ; //!                                                                              
  std::vector<float> n_parton_higgs_eta; //!                                                                              
  std::vector<float> n_parton_higgs_phi; //!                                                                              
  std::vector<float> n_parton_higgs_m  ; //!                   
  int n_parton_higgs_n  ; //!                   

  std::vector<float> n_parton_u_pt ; //!                                                                              
  std::vector<float> n_parton_u_eta; //!                                                                              
  std::vector<float> n_parton_u_phi; //!                                                                              
  std::vector<float> n_parton_u_m  ; //!                   
  int n_parton_u_n  ; //!                  
  std::vector<float> n_parton_au_pt ; //!                                                                              
  std::vector<float> n_parton_au_eta; //!                                                                              
  std::vector<float> n_parton_au_phi; //!                                                                              
  std::vector<float> n_parton_au_m  ; //!                   
  int n_parton_au_n  ; //!                 
  std::vector<float> n_parton_d_pt ; //!                                                                              
  std::vector<float> n_parton_d_eta; //!                                                                              
  std::vector<float> n_parton_d_phi; //!                                                                              
  std::vector<float> n_parton_d_m  ; //!                   
  int n_parton_d_n  ; //!                   
  std::vector<float> n_parton_ad_pt ; //!                                                                              
  std::vector<float> n_parton_ad_eta; //!                                                                              
  std::vector<float> n_parton_ad_phi; //!                                                                              
  std::vector<float> n_parton_ad_m  ; //!                   
  int n_parton_ad_n  ; //!                   
  std::vector<float> n_parton_c_pt ; //!                                                                              
  std::vector<float> n_parton_c_eta; //!                                                                              
  std::vector<float> n_parton_c_phi; //!                                                                              
  std::vector<float> n_parton_c_m  ; //!                   
  int n_parton_c_n  ; //!                    
  std::vector<float> n_parton_ac_pt ; //!                                                                              
  std::vector<float> n_parton_ac_eta; //!                                                                              
  std::vector<float> n_parton_ac_phi; //!                                                                              
  std::vector<float> n_parton_ac_m  ; //!                   
  int n_parton_ac_n  ; //!                   
  std::vector<float> n_parton_s_pt ; //!                                                                              
  std::vector<float> n_parton_s_eta; //!                                                                              
  std::vector<float> n_parton_s_phi; //!                                                                              
  std::vector<float> n_parton_s_m  ; //!                   
  int n_parton_s_n  ; //!                    
  std::vector<float> n_parton_as_pt ; //!                                                                              
  std::vector<float> n_parton_as_eta; //!                                                                              
  std::vector<float> n_parton_as_phi; //!                                                                              
  std::vector<float> n_parton_as_m  ; //!                   
  int n_parton_as_n  ; //!                   
  std::vector<float> n_parton_top_pt ; //!                                                                              
  std::vector<float> n_parton_top_eta; //!                                                                              
  std::vector<float> n_parton_top_phi; //!                                                                              
  std::vector<float> n_parton_top_m  ; //!                   
  int n_parton_top_n  ; //!                   
  std::vector<float> n_parton_atop_pt ; //!                                                                              
  std::vector<float> n_parton_atop_eta; //!                                                                              
  std::vector<float> n_parton_atop_phi; //!                                                                              
  std::vector<float> n_parton_atop_m  ; //!                   
  int n_parton_atop_n  ; //!    
  std::vector<float> n_parton_b_pt ; //!                                                                              
  std::vector<float> n_parton_b_eta; //!                                                                              
  std::vector<float> n_parton_b_phi; //!                                                                              
  std::vector<float> n_parton_b_m  ; //!                   
  int n_parton_b_n  ; //!                   
  std::vector<float> n_parton_ab_pt ; //!                                                                              
  std::vector<float> n_parton_ab_eta; //!                                                                              
  std::vector<float> n_parton_ab_phi; //!                                                                              
  std::vector<float> n_parton_ab_m  ; //!                   
  int n_parton_ab_n  ; //!                   
    
  std::vector<float> n_parton_wplus_pt ; //!                                                                              
  std::vector<float> n_parton_wplus_eta; //!                                                                              
  std::vector<float> n_parton_wplus_phi; //!                                                                              
  std::vector<float> n_parton_wplus_m  ; //!                   
  int n_parton_wplus_n  ; //!    
  std::vector<float> n_parton_wminus_pt ; //!                                                                              
  std::vector<float> n_parton_wminus_eta; //!                                                                              
  std::vector<float> n_parton_wminus_phi; //!                                                                              
  std::vector<float> n_parton_wminus_m  ; //!                   
  int n_parton_wminus_n  ; //!    
  
  //bbHtt and ttHtt collections
  std::vector<float> n_parton_ttfromH_pt ; //!                                                                          
  std::vector<float> n_parton_ttfromH_eta; //!                                                                            
  std::vector<float> n_parton_ttfromH_phi; //!                                                                            
  std::vector<float> n_parton_ttfromH_m  ; //!                   
  int n_parton_ttfromH_n  ; //!    
  std::vector<float> n_parton_ttfromH_poscosphi  ; //!                   
  std::vector<float> n_parton_ttfromH_deltaeta  ; //!                   
  std::vector<float> n_parton_ttfromH_deltaphi  ; //!                   
  std::vector<float> n_parton_ttfromH_deltar  ; //!                   

  //ttHtt specific
  std::vector<float> n_parton_ttnotfromH_pt ; //!                                                                          
  std::vector<float> n_parton_ttnotfromH_eta; //!                                                                        
  std::vector<float> n_parton_ttnotfromH_phi; //!                                                                        
  std::vector<float> n_parton_ttnotfromH_m  ; //!                   
  int n_parton_ttnotfromH_n  ; //!    
  std::vector<float> n_parton_ttnotfromH_deltaeta  ; //!                   
  std::vector<float> n_parton_ttnotfromH_deltaphi  ; //!                   
  std::vector<float> n_parton_ttnotfromH_deltar  ; //!                   

  //bbHtt specific
  std::vector<float> n_parton_bbnotfromH_pt ; //!                                                                          
  std::vector<float> n_parton_bbnotfromH_eta; //!                                                                        
  std::vector<float> n_parton_bbnotfromH_phi; //!                                                                        
  std::vector<float> n_parton_bbnotfromH_m  ; //!                   
  int n_parton_bbnotfromH_n  ; //!    
  std::vector<float> n_parton_bbnotfromH_deltaeta  ; //!                   
  std::vector<float> n_parton_bbnotfromH_deltaphi  ; //!                   
  std::vector<float> n_parton_bbnotfromH_deltar  ; //!                   

  //tcHbb specific
  std::vector<float> n_parton_tcnotfromH_pt ; //!                                                                          
  std::vector<float> n_parton_tcnotfromH_eta; //!                                                                        
  std::vector<float> n_parton_tcnotfromH_phi; //!                                                                        
  std::vector<float> n_parton_tcnotfromH_m  ; //!                   
  int n_parton_tcnotfromH_n  ; //!    
  std::vector<float> n_parton_tcnotfromH_deltaeta  ; //!                   
  std::vector<float> n_parton_tcnotfromH_deltaphi  ; //!                   
  std::vector<float> n_parton_tcnotfromH_deltar  ; //!                   
  /*
  std::vector<float> n_parton_fromH_pt ; //!                                                                          
  std::vector<float> n_parton_fromH_eta; //!                                                                            
  std::vector<float> n_parton_fromH_phi; //!                                                                            
  std::vector<float> n_parton_fromH_m  ; //!  
  */
  //Alvaro
  
  std::vector<float> n_parton_fromH_pt_p ; //!                                                                          
  std::vector<float> n_parton_fromH_eta_p; //!                                                                            
  std::vector<float> n_parton_fromH_phi_p; //!                                                                            
  std::vector<float> n_parton_fromH_m_p  ; //!  
  
  
  
  std::vector<float> n_parton_fromH_pt_a ; //!                                                                          
  std::vector<float> n_parton_fromH_eta_a; //!                                                                            
  std::vector<float> n_parton_fromH_phi_a; //!                                                                            
  std::vector<float> n_parton_fromH_m_a  ; //!  
  //
  int n_parton_fromH_n  ; //!    
  std::vector<float> n_parton_fromH_deltaeta  ; //!                   
  std::vector<float> n_parton_fromH_deltaphi  ; //!                   
  std::vector<float> n_parton_fromH_deltar    ; //!                   

  //tXq specific
  float n_parton_sisterH_pt; //!
  float n_parton_sisterH_eta; //!
  float n_parton_sisterH_phi; //!
  float n_parton_sisterH_m; //!
  int n_parton_sisterH_pdgid; //!
  float n_parton_fromothertop_pt; //!
  float n_parton_fromothertop_eta; //!
  float n_parton_fromothertop_phi; //!
  float n_parton_fromothertop_m; //!
  int n_parton_fromothertop_pdgid; //!
  
  int tXq_nmatched; //!
  int tXq_MultipleMatching; //!
  int tXq_b1Index; //!
  int tXq_b2Index; //!
  int tXq_b3Index; //!
  int tXq_qIndex; //!

  float b1b2dR;//!
  float b1b2_boosted_dR;//!
  
  float JJ_Boosted_dR;//!
  float JJ_boost_dR;//!
  float cos_JJ_Boosted_dR;//!
  float cos_JJ_boost_dR;//!
  
  float e_ae_Boosted_dR;//!
  float e_ae_boost_dR;//!
  float cos_e_ae_Boosted_dR;//!
  float cos_e_ae_boost_dR;//!
  
  float ae_e_Boosted_dR;//!
  float ae_e_boost_dR;//!
  float cos_ae_e_Boosted_dR;//!
  float cos_ae_e_boost_dR;//!
  
  float e_ae_Boosted_dPhi;//!
  float e_ae_boost_dPhi;//!
  float cos_e_ae_Boosted_dPhi;//!
  float cos_e_ae_boost_dPhi;//!
  
  float ae_e_Boosted_dPhi;//!
  float ae_e_boost_dPhi;//!
  float cos_ae_e_Boosted_dPhi;//!
  float cos_ae_e_boost_dPhi;//!
  
  float e_ae_Boosted_dEta;//!
  float e_ae_boost_dEta;//!
  float cos_e_ae_Boosted_dEta;//!
  float cos_e_ae_boost_dEta;//!
  
  float ae_e_Boosted_dEta;//!
  float ae_e_boost_dEta;//!
  float cos_ae_e_Boosted_dEta;//!
  float cos_ae_e_boost_dEta;//!
  
  float Phi;
  float Cos_Phi;
  
  float JJ_Boosted_dPhi;//!
  float JJ_boost_dPhi;//!
  float cos_JJ_Boosted_dPhi;//!
  float cos_JJ_boost_dPhi;//!
  
  float JJ_Boosted_dEta;//!
  float JJ_boost_dEta;//!
  float cos_JJ_Boosted_dEta;//!
  float cos_JJ_boost_dEta;//!
  
  float Phi_plane_angle;//!
  
  float b1b3dR;//!
  float b2b3dR;//!
  float b1qSdR;
  float b2qSdR;
  float minDRbb;//!
  float dR_Jpsi_K;//!

  //Low mass Hplus specific
  std::vector<float> n_parton__pt ; //!                                                                          
  std::vector<float> n_parton__eta; //!                                                                        
  std::vector<float> n_parton__phi; //!                                                                        
  std::vector<float> n_parton__m  ; //!                   
  int n_parton__n  ; //!    
  std::vector<float> n_parton__deltaeta  ; //!                   
  std::vector<float> n_parton__deltaphi  ; //!                   
  std::vector<float> n_parton__deltar  ; //!                   

  //Zb collections
  std::vector<float> n_Zb_Z_pt  ; //!                   
  std::vector<float> n_Zb_Z_eta ; //!                   
  std::vector<float> n_Zb_Z_phi ; //!                   
  std::vector<float> n_Zb_Z_m   ; //!                   
  std::vector<float> n_Zb_wzjets_pt ; //!                   
  std::vector<float> n_Zb_wzjets_eta; //!                   
  std::vector<float> n_Zb_wzjets_phi; //!                   
  std::vector<float> n_Zb_wzjets_m  ; //!                   
  int  n_Zb_wzjets_n  ; //!
  
  //B+ collections
  std::vector<float> n_B_pt;//!
  std::vector<float> n_B_eta;//!
  std::vector<float> n_B_phi;//!
  std::vector<float> n_B_m;//!
  
  std::vector<float> n_J_pt;//!
  std::vector<float> n_J_eta;//!
  std::vector<float> n_J_phi;//!
  std::vector<float> n_J_m;//!
  
  std::vector<float> n_K_pt;//!
  std::vector<float> n_K_eta;//!
  std::vector<float> n_K_phi;//!
  std::vector<float> n_K_m;//!
  
  std::vector<float> n_K_J_dphi;//!
  std::vector<float> n_K_J_deta;//!
  

  std::vector<int> n_parton_H_decays ; //!                   

  std::vector<float> n_parton_topdecaychain_lepton_pt; //!    
  std::vector<float> n_parton_topdecaychain_lepton_eta; //!    
  std::vector<float> n_parton_topdecaychain_lepton_phi; //!    
  std::vector<float> n_parton_topdecaychain_lepton_m; //!    

  std::vector<float> n_parton_topdecaychain_neutrino_pt; //!    
  std::vector<float> n_parton_topdecaychain_neutrino_eta; //!    
  std::vector<float> n_parton_topdecaychain_neutrino_phi; //!    
  std::vector<float> n_parton_topdecaychain_neutrino_m; //!    

  float n_HT; //!

  float n_angana_lep_top_dphi; //!
  float n_angana_lep_top_dr; //!
  float n_angana_nu_top_dphi; //!
  float n_angana_nu_top_dr; //!
  float n_angana_lep_top_dphi_topcm; //!
  float n_angana_nu_top_dphi_topcm; //!

  void SetNtupleVars(TruthAna_OptionParser *opt, TTree *tree);
  void FillParticleCollection(std::vector<TLorentzVector> vtlv, std::string coll_name);
  void FillScalarCollection(float variable, std::string coll_name);
  void FillScalarCollection(int variable, std::string coll_name);
  void FillHCollection(std::vector<TLorentzVector> vtlv, std::string coll_name);
  void FillHCollection(std::vector<myParticle> vtlv, std::string coll_name);
  void FillHCollection(myParticle vtlv, std::string coll_name);
  //Temporary implementation
  void FillEventCollection(float mc_event_weight);
  void FillEventCollectionHack(std::vector<int> parton_H_decays);

  void ClearCollections();

 private:
  TruthAna_OptionParser *opt;

};

#endif 
