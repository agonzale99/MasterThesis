#ifndef MyTruthAnalysis_MyxAODTruthAnalysis_H
#define MyTruthAnalysis_MyxAODTruthAnalysis_H

#include <EventLoop/Algorithm.h>

#include "xAODTruth/TruthParticleContainer.h"
#include "xAODJet/JetContainer.h"
#include "xAODEventInfo/EventInfo.h"

#include <TTree.h>
#include <TMath.h>

#include "MyTruthAnalysis/MyClasses.h"
#include "MyTruthAnalysis/TruthAna_OptionParser.h"
#include "MyTruthAnalysis/TruthAna_NtupleManager.h"
#include "MyTruthAnalysis/TruthAna_GenericSelector.h"
#include "MyTruthAnalysis/TruthAna_ttHttSelector.h"
#include "MyTruthAnalysis/TruthAna_bbHttSelector.h"
#include "MyTruthAnalysis/TruthAna_tcHbbSelector.h"
#include "MyTruthAnalysis/TruthAna_LowMassHplusSelector.h"
#include "MyTruthAnalysis/TruthAna_ZbSelector.h"
#include "MyTruthAnalysis/TruthAna_CutFlowTools.h"

class TruthAna_OptionParser;
class TruthAna_NtupleManager;
class TruthAna_GenericSelector;
class TruthAna_ttHttSelector;
class TruthAna_bbHttSelector;
class TruthAna_tcHbbSelector;
class TruthAna_LowMassHplusSelector;
class TruthAna_ZbSelector;
class TruthAna_CutFlowTools;

class MyxAODTruthAnalysis : public EL::Algorithm
{

public:

  MyxAODTruthAnalysis ();

  enum class PROCESS{TT_H_TT,BB_H_TT,TC_H_BB,ZB,TT,LOWMASSHPLUS,GENERIC};
  enum class DERIVATION{TRUTH0,TRUTH1};

  PROCESS proc;
  DERIVATION deriv;

  // these are the functions inherited from Algorithm
  virtual EL::StatusCode setupJob (EL::Job& job);
  virtual EL::StatusCode fileExecute ();
  virtual EL::StatusCode histInitialize ();
  virtual EL::StatusCode changeInput (bool firstFile);
  virtual EL::StatusCode initialize ();
  virtual EL::StatusCode execute ();
  virtual EL::StatusCode postExecute ();
  virtual EL::StatusCode finalize ();
  virtual EL::StatusCode histFinalize ();

  std::string outputName;

  EL::StatusCode LoadContainers();
  EL::StatusCode ParticleFiller();
  EL::StatusCode NtupleFiller();
  EL::StatusCode ClearParticles();
  bool CutFlow();
  
  TTree *tree; //!
  const xAOD::EventInfo              *m_eventInfoCont; //!                                                         
  const xAOD::TruthParticleContainer *m_TruthParticles; //!                                                      
  const xAOD::TruthParticleContainer *m_TruthNeutrinos; //!                                                      
  const xAOD::TruthParticleContainer *m_TruthMuons; //!                                                          
  const xAOD::TruthParticleContainer *m_TruthElectrons; //!                                                      
  const xAOD::TruthParticleContainer *m_TruthTaus; //!                                                           
  const xAOD::JetContainer           *m_AntiKt4TruthWZJets; //!    
  
  //
  
  const xAOD::TruthParticleContainer *m_LabelBHadrons_in; //! 
  
  const xAOD::TruthParticleContainer_v1 *m_LabelBHadrons_fin; //! 
  
  
  
  

  TruthAna_OptionParser* m_OptionParser; //! 
  TruthAna_NtupleManager* m_NtupleManager; //!
  TruthAna_GenericSelector* m_GenericSelector; //!
  TruthAna_ttHttSelector* m_ttHttSelector; //!
  TruthAna_bbHttSelector* m_bbHttSelector; //!
  TruthAna_tcHbbSelector* m_tcHbbSelector; //!
  TruthAna_LowMassHplusSelector* m_LowMassHplusSelector; //!
  TruthAna_ZbSelector* m_ZbSelector; //!
  TruthAna_CutFlowTools* m_CutFlowTools; //!

  std::vector<int> v_leptons_pdgid; //!
  std::vector<int> v_neutrinos_pdgid; //!

  std::vector<TLorentzVector> v_neutrinos; //! 
  std::vector<TLorentzVector> v_muons; //! 
  std::vector<TLorentzVector> v_electrons; //!
  std::vector<TLorentzVector> v_taus; //! 
  std::vector<TLorentzVector> v_AntiKt4TruthWZJets; //!  
  std::vector<TLorentzVector> v_particle_wboson; //! 
  std::vector<TLorentzVector> v_parton_u; //!  
  std::vector<TLorentzVector> v_parton_au; //!  
  std::vector<TLorentzVector> v_parton_d; //!  
  std::vector<TLorentzVector> v_parton_ad; //!  
  std::vector<TLorentzVector> v_parton_c; //!  
  std::vector<TLorentzVector> v_parton_ac; //!  
  std::vector<TLorentzVector> v_parton_s; //!  
  std::vector<TLorentzVector> v_parton_as; //!  
  std::vector<TLorentzVector> v_parton_top; //!  
  std::vector<TLorentzVector> v_parton_atop; //!  
  std::vector<TLorentzVector> v_parton_b; //!  
  std::vector<TLorentzVector> v_parton_ab; //!  
  std::vector<TLorentzVector> v_parton_wplus; //!  
  std::vector<TLorentzVector> v_parton_wminus; //!  
  std::vector<TLorentzVector> v_parton_higgs; //!  
  std::vector<TLorentzVector> v_parton_ttfromH; //!  
  std::vector<TLorentzVector> v_parton_ttnotfromH; //!  
  std::vector<TLorentzVector> v_parton_bbfromH; //!  
  std::vector<myParticle> v_parton_tcnotfromH; //!  
  myParticle v_parton_sisterH; //!  
  myParticle v_parton_fromothertop; //!
  std::vector<TLorentzVector> v_parton_Hplus_decays; //!    
  std::vector<TLorentzVector> v_Zb_Z; //!  
  std::vector<TLorentzVector> v_Zb_wzjets; //!  
  std::vector<TLorentzVector> v_parton_fromH; //!  
  std::vector<int> v_H_decay_pdgid; //!  
  std::vector<TLorentzVector> v_parton_topdecaychain_lepton; //!    
  std::vector<TLorentzVector> v_parton_topdecaychain_neutrino; //!    

  std::vector< std::pair< TLorentzVector,int > > v_TopDecayChain; //!    

//////Alvaro
  std::vector<TLorentzVector> v_B0_hadrons; //! 
  std::vector<TLorentzVector> v_Bplus_hadrons_decay; //! 
  std::vector<TLorentzVector> v_Jpsi_fromB;//!
  std::vector<TLorentzVector> v_B_decay;//!
  std::vector<TLorentzVector> v_K_fromB;//!
  
  std::vector<TLorentzVector> B_plus;//!
  std::vector<TLorentzVector> B_minus;//!
  
  std::vector<TLorentzVector> J_aB;//!
  std::vector<TLorentzVector> J_B;//!
  
  std::vector<TLorentzVector> e_aB;//! electron
  std::vector<TLorentzVector> e_B;//! 
  
  std::vector<TLorentzVector> ae_aB;//! positron
  std::vector<TLorentzVector> ae_B;//! 
  
  
  
  
  std::vector<TLorentzVector> J_aB_boost;//!al B frame
  std::vector<TLorentzVector> J_B_boost;//!
  
  std::vector<TLorentzVector> e_aB_boost;//!
  std::vector<TLorentzVector> e_B_boost;//!al B frame
  
  std::vector<TLorentzVector> ae_aB_boost;//!
  std::vector<TLorentzVector> ae_B_boost;//!al B frame
  
  std::vector<TLorentzVector> J_aB_Boosted;//!
  std::vector<TLorentzVector> J_B_Boosted;//!al X frame
  
  std::vector<TLorentzVector> e_aB_Boosted;//!
  std::vector<TLorentzVector> e_B_Boosted;//!al X frame
  
  std::vector<TLorentzVector> ae_aB_Boosted;//!
  std::vector<TLorentzVector> ae_B_Boosted;//!al X frame
  
  std::vector<TLorentzVector> K_aB;//!
  std::vector<TLorentzVector> K_B;//!
  
  float mc_event_weight; //!
  float v_HT; //!
  //H+ angular analysis 
  float v_angana_lep_top_dphi; //!
  float v_angana_lep_top_dr; //!
  float v_angana_nu_top_dphi; //!
  float v_angana_nu_top_dr; //!
  float v_angana_lep_top_dphi_topcm; //!
  float v_angana_nu_top_dphi_topcm; //!

  int nJets; //!
  int nmatched; //!
  int MultipleMatching; //!
  int b1_JetIndex; //!
  int b2_JetIndex; //!
  int b3_JetIndex; //!
  int qS_JetIndex; //!

  float b1b2dR;//!
  float b1b2_boosted_dR;//!
  
  float JJ_boosted_dR;//!
  float cos_JJ_boosted_dR;//!
  float JJ_Boosted_dR;//!
  float cos_JJ_Boosted_dR;//!
  float JJ_boost_dR;//!
  float cos_JJ_boost_dR;//!
  
  float JJ_boosted_dPhi;//!
  float cos_JJ_boosted_dPhi;//!
  float JJ_Boosted_dPhi;//!
  float cos_JJ_Boosted_dPhi;//!
  float JJ_boost_dPhi;//!
  float cos_JJ_boost_dPhi;//!
  
  float JJ_boosted_dEta;//!
  float cos_JJ_boosted_dEta;//!
  float JJ_Boosted_dEta;//!
  float cos_JJ_Boosted_dEta;//!
  float JJ_boost_dEta;//!
  float cos_JJ_boost_dEta;//!
  
  //first lepton coming from B secondo from aB
  //X frame
  float e_ae_Boosted_dR;
  float ae_e_Boosted_dR;
  float cos_e_ae_Boosted_dR;
  float cos_ae_e_Boosted_dR;
  
  float e_ae_Boosted_dPhi;
  float ae_e_Boosted_dPhi;
  float cos_e_ae_Boosted_dPhi;
  float cos_ae_e_Boosted_dPhi;
  
  float e_ae_Boosted_dEta;
  float ae_e_Boosted_dEta;
  float cos_e_ae_Boosted_dEta;
  float cos_ae_e_Boosted_dEta;
  
  
  //Bframe
  float e_ae_boost_dR;
  float ae_e_boost_dR;
  float cos_e_ae_boost_dR;
  float cos_ae_e_boost_dR;
  
  float e_ae_boost_dPhi;
  float ae_e_boost_dPhi;
  float cos_e_ae_boost_dPhi;
  float cos_ae_e_boost_dPhi;
  
  float e_ae_boost_dEta;
  float ae_e_boost_dEta;
  float cos_e_ae_boost_dEta;
  float cos_ae_e_boost_dEta;
  
  //non boosted
  float e_ae_dR;
  float ae_e_dR;
  
  float e_ae_dPhi;
  float ae_e_dPhi;
  
  float e_ae_dEta;
  float ae_e_dEta;
  
  
  
  float Phi;
  float Cos_Phi;
  
  
  
  
  float Phi_plane_angle;//!  
  float Theta_plane_angle;//!
  
  
  float b1b3dR;//!
  float b2b3dR;//!
  float b1qSdR;//!
  float b2qSdR;//!
  float minDRbb;//!
  
  float dR_B_Jpsi;//!
  float dR_B_K;//!
  float dR_Jpsi_K;//!

  // this is needed to distribute the algorithm to the workers
  ClassDef(MyxAODTruthAnalysis, 1);

};

#endif
