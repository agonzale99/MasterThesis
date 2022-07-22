#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <MyTruthAnalysis/MyxAODTruthAnalysis.h>

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "xAODRootAccess/tools/Message.h"

// EDM includes:
#include "xAODEventInfo/EventInfo.h"

#include <TFile.h>
#include <TMath.h>

// this is needed to distribute the algorithm to the workers
ClassImp(MyxAODTruthAnalysis)

/// Helper macro for checking xAOD::TReturnCode return values
#define EL_RETURN_CHECK( CONTEXT, EXP )                     \
  do {							    \
    if( ! EXP.isSuccess() ) {				    \
      Error( CONTEXT,					    \
	     XAOD_MESSAGE( "Failed to execute: %s" ),	    \
	     #EXP );					    \
      return EL::StatusCode::FAILURE;			    \
    }							    \
  } while( false )



MyxAODTruthAnalysis :: MyxAODTruthAnalysis()
{
  proc =PROCESS::TC_H_BB; //LOWMASSHPLUS;
  deriv=DERIVATION::TRUTH1;

  m_OptionParser         = new TruthAna_OptionParser;  
  m_NtupleManager        = new TruthAna_NtupleManager; 
  m_GenericSelector      = new TruthAna_GenericSelector; 
  m_ttHttSelector        = new TruthAna_ttHttSelector; 
  m_bbHttSelector        = new TruthAna_bbHttSelector; 
  m_tcHbbSelector        = new TruthAna_tcHbbSelector; 
  m_LowMassHplusSelector = new TruthAna_LowMassHplusSelector; 
  m_ZbSelector           = new TruthAna_ZbSelector; 
  m_CutFlowTools         = new TruthAna_CutFlowTools; 

  v_leptons_pdgid.push_back(11); 
  v_leptons_pdgid.push_back(13); 
  v_leptons_pdgid.push_back(15); 
  v_neutrinos_pdgid.push_back(-12);  
  v_neutrinos_pdgid.push_back(-14);  
  v_neutrinos_pdgid.push_back(-16);  

  m_CutFlowTools->addCutFlow("CutFlow","All|OneMu|OneEl|OneTau|OneNu|Lep+Nu|FourJets|FiveJets|>=SixJets");
  m_CutFlowTools->addCutFlow("CutFlow weighted","All|OneMu|OneEl|OneTau|OneNu|Lep+Nu|FourJets|FiveJets|>=SixJets");
  
}



EL::StatusCode MyxAODTruthAnalysis :: setupJob (EL::Job& job)
{
  job.useXAOD ();
  EL_RETURN_CHECK( "setupJob()", xAOD::Init() ); // call before opening first file

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: histInitialize ()
{
  TFile *outputFile = wk()->getOutputFile (outputName);
  tree = new TTree ("tree", "tree");
  tree->SetDirectory (outputFile);
  
  m_NtupleManager->SetNtupleVars(m_OptionParser,tree);

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: fileExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: changeInput (bool firstFile)
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: initialize ()
{
  //Just print infos
  xAOD::TEvent* event = wk()->xaodEvent(); 
  Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: execute ()
{
  
  //Info("initialize()", "Number of events = %lli", event->getEntries() ); // print long long int
  
  //Clear vectors to be saved in the tree
  m_NtupleManager->ClearCollections();
  //Load xAOD containers
  LoadContainers();
  //Fill TLorentzVectors of all particles
  ParticleFiller();
  //Fill ntuple
  NtupleFiller();
  //Checking some basic cuts, may use it to enable the tree filling
  CutFlow();
  //Clear vectors to be saved in the tree
  ClearParticles();
  
  tree->Fill();

  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: postExecute ()
{
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: finalize ()
{
  m_CutFlowTools->printAllCutFlow();
  return EL::StatusCode::SUCCESS;
}



EL::StatusCode MyxAODTruthAnalysis :: histFinalize ()
{
  return EL::StatusCode::SUCCESS;
}


EL::StatusCode MyxAODTruthAnalysis :: LoadContainers()
{
  xAOD::TEvent* event = wk()->xaodEvent(); // you should have already added this as described before
  m_eventInfoCont           =0; EL_RETURN_CHECK("execute()",event->retrieve(m_eventInfoCont     , "EventInfo"));
  m_TruthParticles          =0; EL_RETURN_CHECK("execute()",event->retrieve(m_TruthParticles    , "TruthParticles"));
 //   
  //m_LabelBHadrons_in          =0; EL_RETURN_CHECK("execute()",event->retrieve(m_LabelBHadrons_in    , "TruthLabelBHadronsInitial"));  
  //m_LabelBHadrons_fin          =0; EL_RETURN_CHECK("execute()",event->retrieve(m_LabelBHadrons_fin    , "TruthLabelBHadronsInitial"));    

  mc_event_weight=m_eventInfoCont->mcEventWeight();

  if(deriv==DERIVATION::TRUTH1)
    {
      m_TruthNeutrinos      =0; EL_RETURN_CHECK("execute()",event->retrieve(m_TruthNeutrinos    ,"TruthNeutrinos"));
      m_TruthMuons          =0; EL_RETURN_CHECK("execute()",event->retrieve(m_TruthMuons        ,"TruthMuons"));       
      m_TruthElectrons      =0; EL_RETURN_CHECK("execute()",event->retrieve(m_TruthElectrons    ,"TruthElectrons"));   
      m_TruthTaus           =0; EL_RETURN_CHECK("execute()",event->retrieve(m_TruthTaus         ,"TruthTaus"));         
      m_AntiKt4TruthWZJets  =0; EL_RETURN_CHECK("execute()",event->retrieve(m_AntiKt4TruthWZJets,"AntiKt4TruthDressedWZJets")); //AntiKt4TruthWZJets ADRI
    }

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode MyxAODTruthAnalysis :: ParticleFiller ()
{
  
  //Fills with the GenericSelector
  //TRUTH1 ONLY
  if(deriv==DERIVATION::TRUTH1)
    {
      v_neutrinos = m_GenericSelector->GenericParticle(m_TruthNeutrinos,(float)25000.,(float)10.,
						       m_GenericSelector->CUT_PDG::NO,-1,m_GenericSelector->CUT_STATUS::NO,-1);
      v_muons     = m_GenericSelector->GenericParticle(m_TruthMuons,(float)25000.,(float)2.5,
						       m_GenericSelector->CUT_PDG::NO,-1,m_GenericSelector->CUT_STATUS::NO,-1);
      v_electrons = m_GenericSelector->GenericParticle(m_TruthElectrons,(float)25000.,(float)2.5,
						       m_GenericSelector->CUT_PDG::NO,-1,m_GenericSelector->CUT_STATUS::NO,-1);
      v_taus      = m_GenericSelector->GenericParticle(m_TruthTaus,(float)25000.,(float)2.5,
						       m_GenericSelector->CUT_PDG::NO,-1,m_GenericSelector->CUT_STATUS::NO,-1);
      
      v_AntiKt4TruthWZJets = m_GenericSelector->GenericJet(m_AntiKt4TruthWZJets,(float)25000.,(float)2.5);
      sort(v_AntiKt4TruthWZJets.begin(), v_AntiKt4TruthWZJets.end(), ComparePt());

      v_HT =  m_GenericSelector->CalculateHT(v_AntiKt4TruthWZJets, v_electrons, v_muons, v_taus, v_neutrinos);

      //If run on Truth1 you won't need to go any further
      //return EL::StatusCode::SUCCESS; //ADRI
    }

  //If run on Truth0 derivations keep going as you will need to access more information and store additional variables

  //Final state matrix-element particles

   /* int counter = -1;
    for(auto vcont : *m_TruthParticles){
      if (vcont->auxdata<int>("status") == 22 || vcont->auxdata<int>("status") == 23){
      counter++;
      int par_pdgid=vcont->auxdata<int>("pdgId");
      float loc_px    =vcont->px(); float loc_py    =vcont->py();
      float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
      TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e);
      std::cout << "PARTICLE " << counter << " " << par_pdgid << " " << vcont->auxdata<int>("status") << " " << 0.001*tlv.Pt() << " " << 0.001*tlv.M() << std::endl;    
      }
    }*/


  std::vector< std::pair< TLorentzVector,int > > AllStablePartons = m_GenericSelector->GetAllPartons(m_TruthParticles, m_GenericSelector->CUT_STATUS::YES, 23);
  v_parton_u    = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES,  2);
  v_parton_au   = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES, -2);
  v_parton_d    = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES,  1);
  v_parton_ad   = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES, -1);
  v_parton_c    = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES,  4);
  v_parton_ac   = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES, -4);
  v_parton_s    = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES,  3);
  v_parton_as   = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES, -3);
  v_parton_b    = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES,  5);
  v_parton_ab   = m_GenericSelector->RetrieveParton(AllStablePartons, m_GenericSelector->CUT_PDG::YES, -5);

  //Intermediate particles
  std::vector< std::pair< TLorentzVector,int > > AllIntermediatePartons = m_GenericSelector->GetAllPartons(m_TruthParticles, m_GenericSelector->CUT_STATUS::YES, 22);
  v_parton_top   = m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::YES,  6 );
  v_parton_atop  = m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::YES, -6 );
  v_parton_wplus = m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::YES,  24);
  v_parton_wminus= m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::YES, -24);
  //35 id ALP
   
    
   //Alvaro
   //Hadrons
   //std::vector< std::pair< TLorentzVector,int > > v_initial_B_hadrons= m_GenericSelector->GetAllPartons(m_LabelBHadrons_in, m_GenericSelector->CUT_STATUS::NO, 23); 
    //std::vector< std::pair< TLorentzVector,int > > v_final_B_hadrons= m_GenericSelector->GetAllPartons(m_LabelBHadrons_fin, m_GenericSelector->CUT_STATUS::NO, 23); 
  //std::cout << "hola" << std::endl;
    //v_final_B_hadrons = m_GenericSelector->GenericParticle(m_LabelBHadrons_fin,(float)25000.,(float)10.,
//						       m_GenericSelector->CUT_PDG::NO,-1,m_GenericSelector->CUT_STATUS::NO,-1);
    
    //std::vector< std::pair< TLorentzVector,int > > AllStablePartons = m_GenericSelector->GetAllPartons(m_TruthParticles, m_GenericSelector->CUT_STATUS::YES, 23);
    
    //std::vector< std::pair< TLorentzVector,int > > AllPartons = m_GenericSelector->GetAllPartons(m_TruthParticles, m_GenericSelector->CUT_STATUS::NO, 23);
    //v_B0_hadrons    = m_GenericSelector->RetrieveParton(AllPartons, m_GenericSelector->CUT_PDG::YES,  511);
    //v_Bplus_hadrons    = m_GenericSelector->RetrieveParton(AllPartons, m_GenericSelector->CUT_PDG::YES,  521);
    
    
  if(proc==PROCESS::TT_H_TT ||  proc==PROCESS::BB_H_TT)v_parton_higgs = m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::ABS,  25);
  else if(proc==PROCESS::TC_H_BB)                      v_parton_higgs = m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::ABS,  25);
  else if(proc==PROCESS::LOWMASSHPLUS)                 v_parton_higgs = m_GenericSelector->RetrieveParton(AllIntermediatePartons, m_GenericSelector->CUT_PDG::ABS,  37);

  //Process specific fills 
  if(proc==PROCESS::TT_H_TT)
    {
      v_parton_ttfromH    = m_ttHttSelector->DecayH(m_TruthParticles,(float)0.,(float)10.,m_ttHttSelector->DECAY_H::ALL);
      v_parton_ttnotfromH = m_ttHttSelector->AssociatedToH(m_TruthParticles,(float)0.,(float)10.,m_ttHttSelector->ASSOCIATED_H::TT);
    }

  else if(proc==PROCESS::BB_H_TT)
    {
      v_parton_ttfromH    = m_bbHttSelector->DecayH(m_TruthParticles,(float)0.,(float)10.,m_bbHttSelector->DECAY_H::ALL);
      v_H_decay_pdgid     = m_bbHttSelector->GetHiggsDecays(m_TruthParticles,25);
      //to be added
      //Implemetation has to make sure that there are no b-quarks from the chains H->b->gb->ggb..
      //v_parton_bbnotfromH = m_bbHttSelector->AssociatedToH(m_TruthParticles,(float)0.,(float)10.,m_bbHttSelector->ASSOCIATED_H::BB);
    }

  else if(proc==PROCESS::TC_H_BB)
    {
      //W boson particle level 
      if(deriv==DERIVATION::TRUTH1 && v_electrons.size()!=0 && v_neutrinos.size()!=0 )
	v_particle_wboson= m_tcHbbSelector->GetParticleW(v_electrons,v_neutrinos);
      myParticle testing = myParticle();
      TLorentzVector test; test.SetPxPyPzE(1.,1.,1.,1.);
      testing.FourVector = test;
      v_parton_bbfromH    = m_tcHbbSelector->DecayH(m_TruthParticles,(float)0.,(float)10.,m_tcHbbSelector->DECAY_H::ALL);
      //sort(v_parton_bbfromH.begin(), v_parton_bbfromH.end(), ComparePt());
      //std::cout << "DEBUG TC_H_BB: " << v_parton_bbfromH[0].Pt()*0.001 << " " << v_parton_bbfromH[0].Eta() << " " <<v_parton_bbfromH[0].Phi() << " " << v_parton_bbfromH[1].Pt()*0.001 << " " << v_parton_bbfromH[1].Eta() << " " <<v_parton_bbfromH[1].Phi() << std::endl;
      v_parton_sisterH    = m_tcHbbSelector->SisterH(m_TruthParticles);
      //std::cout << "->" << v_parton_sisterH.pdgid << " " << v_parton_sisterH.FourVector.Pt()*0.001 << " " << v_parton_sisterH.FourVector.Eta() << " " << v_parton_sisterH.FourVector.Phi()  << std::endl;
      v_parton_fromothertop = m_tcHbbSelector->QuarkOtherTop(m_TruthParticles);   
      //std::cout << "->" << v_parton_fromothertop.pdgid << " " << v_parton_fromothertop.FourVector.Pt()*0.001 << " " << v_parton_fromothertop.FourVector.Eta() << " " << v_parton_fromothertop.FourVector.Phi() << std::endl;
      v_parton_tcnotfromH = m_tcHbbSelector->AssociatedToH_P(m_TruthParticles,(float)0.,(float)10.,m_tcHbbSelector->ASSOCIATED_H::TC);
      v_H_decay_pdgid     = m_tcHbbSelector->GetHiggsDecays(m_TruthParticles,25);

      nJets = v_AntiKt4TruthWZJets.size();

      if (v_parton_sisterH.FourVector.Pt() == v_parton_fromothertop.FourVector.Pt() || v_parton_bbfromH[0].Pt()==v_parton_sisterH.FourVector.Pt() || v_parton_bbfromH[0].Pt()==v_parton_fromothertop.FourVector.Pt() || v_parton_bbfromH[1].Pt()==v_parton_sisterH.FourVector.Pt() || v_parton_bbfromH[1].Pt()==v_parton_fromothertop.FourVector.Pt()){
        std::cout << "ERROR TWO PARTICLES EQUAL" << std::endl;
        exit(0);
      }

      //std::cout << "Jets " << nJets << std::endl;
      MultipleMatching = 0;
      b1_JetIndex =  -1;
      b2_JetIndex =  -1;
      qS_JetIndex =  -1;
      b3_JetIndex =  -1;

      b1b2dR = 10.0;
      b1b2_boosted_dR = 10.0;
      
      JJ_Boosted_dR = 10.0;
      JJ_boost_dR = 10.0;
      cos_JJ_Boosted_dR = 10.0;
      cos_JJ_boost_dR = 10.0;
      
      e_ae_Boosted_dR = 10.0;
      e_ae_boost_dR = 10.0;
      cos_e_ae_Boosted_dR = 10.0;
      cos_e_ae_boost_dR = 10.0;
      
      ae_e_Boosted_dR = 10.0;
      ae_e_boost_dR = 10.0;
      cos_ae_e_Boosted_dR = 10.0;
      cos_ae_e_boost_dR = 10.0;
      
      ae_e_dR = 10.0;
      e_ae_dR = 10.0;
      ae_e_dPhi = 10.0;
      e_ae_dPhi = 10.0;
      ae_e_dEta = 10.0;
      e_ae_dEta = 10.0;
      
      e_ae_Boosted_dPhi = 10.0;
      e_ae_boost_dPhi = 10.0;
      cos_e_ae_Boosted_dPhi = 10.0;
      cos_e_ae_boost_dPhi = 10.0;
      
      ae_e_Boosted_dPhi = 10.0;
      ae_e_boost_dPhi = 10.0;
      cos_ae_e_Boosted_dPhi = 10.0;
      cos_ae_e_boost_dPhi = 10.0;
      
      e_ae_Boosted_dEta = 10.0;
      e_ae_boost_dEta = 10.0;
      cos_e_ae_Boosted_dEta = 10.0;
      cos_e_ae_boost_dEta = 10.0;
      
      ae_e_Boosted_dEta = 10.0;
      ae_e_boost_dEta = 10.0;
      cos_ae_e_Boosted_dEta = 10.0;
      cos_ae_e_boost_dEta = 10.0;
      
      Phi = 10.0;
      Cos_Phi = 10.0;
      
      JJ_Boosted_dPhi = 10.0;
      JJ_boost_dPhi = 10.0;
      cos_JJ_Boosted_dPhi = 10.0;
      cos_JJ_boost_dPhi = 10.0;
      
      JJ_Boosted_dEta = 10.0;
      JJ_boost_dEta = 10.0;
      cos_JJ_Boosted_dEta = 10.0;
      cos_JJ_boost_dEta = 10.0;
      
      Phi_plane_angle = 10.0;
      b1b3dR = 10.0;
      b2b3dR = 10.0;
      b1qSdR = 10.0;
      b2qSdR = 10.0;
      minDRbb = 100.0;
      nmatched = 0;
      double DeltaR = 0.3;
      dR_Jpsi_K = 10.0;
      

      double DRb1 =101.0;
      double DRb1min = 100.0;      
      double DRb1Min = 99.0;      
      double DRb2 =101.0;
      double DRb2min = 100.0;
      double DRb2Min = 99.0;
      double DRb3 =101.0;
      double DRb3min = 100.0;
      double DRb3Min = 99.0;
      double DRqS =101.0;
      double DRqSmin = 100.0;
      double DRqSMin = 99.0;
      int index = -1;
      /*for(auto vjet: v_AntiKt4TruthWZJets){
        index++;
        DRb1 = vjet.DeltaR(v_parton_bbfromH[0]);
        DRb2 = vjet.DeltaR(v_parton_bbfromH[1]);
        DRb3 = vjet.DeltaR(v_parton_fromothertop.FourVector);
        DRqS = vjet.DeltaR(v_parton_sisterH.FourVector);

        if (DRb1 < DRb1min) {DRb1min=DRb1;}
        if (DRb2 < DRb2min) {DRb2min=DRb2;}
        if (DRb3 < DRb3min) {DRb3min=DRb3;}
        if (DRqS < DRqSmin) {DRqSmin=DRqS;}

        if (DRb1min < DeltaR && DRb1min < DRb1Min && DRb1min < DRb2min && DRb1min < DRb3min && DRb1min < DRqSmin ){
          b1_JetIndex=index;
        }
      }

      if (DRb1min < 10) { nmatched++;}
      if (DRb2min < 10) { nmatched++;}
      if (DRb2min < 10) { nmatched++;}
      if (DRqSmin < 10) { nmatched++;}
      */
      index = -1;
      for(auto vjet: v_AntiKt4TruthWZJets){
        index++;
      //std::cout << index << " " << vjet.Pt()*0.001 << " " << vjet.Eta() << " " << vjet.Phi() << std::endl;
        //b1->particle  
        if ( vjet.DeltaR(v_parton_bbfromH[0]) < DeltaR)                { b1_JetIndex=index; nmatched++; }
        //b2->antiparticle  
        if ( vjet.DeltaR(v_parton_bbfromH[1]) < DeltaR)                { b2_JetIndex=index; nmatched++; }
        if ( vjet.DeltaR(v_parton_sisterH.FourVector) < DeltaR)        { qS_JetIndex=index; nmatched++; }
        if ( vjet.DeltaR(v_parton_fromothertop.FourVector) < DeltaR)   { b3_JetIndex=index; nmatched++; }
      }

      if (b1_JetIndex >= 0 && b2_JetIndex >= 0 && b1_JetIndex == b2_JetIndex ) MultipleMatching++;
      if (b1_JetIndex >= 0 && qS_JetIndex >= 0 && b1_JetIndex == qS_JetIndex ) MultipleMatching++;
      if (b1_JetIndex >= 0 && b3_JetIndex >= 0 && b1_JetIndex == b3_JetIndex ) MultipleMatching++;
      if (b2_JetIndex >= 0 && qS_JetIndex >= 0 && b2_JetIndex == qS_JetIndex ) MultipleMatching++;
      if (b2_JetIndex >= 0 && b3_JetIndex >= 0 && b2_JetIndex == b3_JetIndex ) MultipleMatching++;
      if (qS_JetIndex >= 0 && b3_JetIndex >= 0 && qS_JetIndex == b3_JetIndex ) MultipleMatching++;

      if (b1_JetIndex >= 0 && b2_JetIndex >= 0 && b1_JetIndex != b2_JetIndex) b1b2dR = v_AntiKt4TruthWZJets[b1_JetIndex].DeltaR(v_AntiKt4TruthWZJets[b2_JetIndex]);
      if (b1_JetIndex >= 0 && b3_JetIndex >= 0 && b1_JetIndex != b3_JetIndex) b1b3dR = v_AntiKt4TruthWZJets[b1_JetIndex].DeltaR(v_AntiKt4TruthWZJets[b3_JetIndex]);
      if (b2_JetIndex >= 0 && b3_JetIndex >= 0 && b2_JetIndex != b3_JetIndex) b2b3dR = v_AntiKt4TruthWZJets[b2_JetIndex].DeltaR(v_AntiKt4TruthWZJets[b3_JetIndex]);
      if (b1_JetIndex >= 0 && qS_JetIndex >= 0 && b1_JetIndex != qS_JetIndex) b1qSdR = v_AntiKt4TruthWZJets[b1_JetIndex].DeltaR(v_AntiKt4TruthWZJets[qS_JetIndex]);
      if (b2_JetIndex >= 0 && qS_JetIndex >= 0 && b2_JetIndex != qS_JetIndex) b2qSdR = v_AntiKt4TruthWZJets[b2_JetIndex].DeltaR(v_AntiKt4TruthWZJets[qS_JetIndex]);

      if (b1b2dR < 10 && b1b2dR < minDRbb) minDRbb = b1b2dR ;
      if (b1b3dR < 10 && b1b3dR < minDRbb) minDRbb = b1b3dR ;
      if (b2b3dR < 10 && b2b3dR < minDRbb) minDRbb = b2b3dR ;

      //std::cout << "matches:" << nmatched << " MM:" << MultipleMatching << " b1:"<< b1_JetIndex << " b2:"<< b2_JetIndex << " qS:" << qS_JetIndex << " b3:" << b3_JetIndex << " minDRbb:" << minDRbb <<  " b1b2dR:" << b1b2dR << " b1b3dR:"<< b1b3dR << " b2b3dR:" << b2b3dR << std::endl;
      
      /*int count = -1;
      const xAOD::TruthParticle* parent;
      std::cout << " " << std::endl;
      for(auto vcont : *m_TruthParticles){
        count++;
        int apdgid = fabs(vcont->auxdata<int>("pdgId"));
        if (!(vcont->auxdata<int>("status")==23 || ( fabs(vcont->auxdata<int>("pdgId"))<=18 && fabs(vcont->auxdata<int>("pdgId"))>=11 ))) continue;
        std::cout << " A: " << count << " " << vcont->pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;
        for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++){
          parent = vcont->parent(par_index);
          if(!parent) continue;
          int parent_pdgId=parent->pdgId();
          std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
        } 
	}*/
      //
      
      //Angular difference in H rest frame of b ab
      TLorentzVector b_p_boosted = v_parton_bbfromH[0]; //b
      TLorentzVector b_a_boosted = v_parton_bbfromH[1]; // anti-b
      TVector3 H_boost_vector = v_parton_bbfromH[2].BoostVector();
      
      b_p_boosted.Boost(-H_boost_vector); //booost b 
      b_a_boosted .Boost(-H_boost_vector); //boost ab
      
      b1b2_boosted_dR = b_p_boosted.DeltaR(b_a_boosted); //computation of angular differnce
      
      
      
      
      
      
      
      
      
      
      //J/psi J/psi angular difference in h/ALP rest frame
      
      
      int Bindex=0;
      int contt=0; //cont to know number of times we have both decays in the same event 
      int conttt=0;
      //std::cout<<std::endl<<"Hola"<<std::endl;
      //std::cout<<std::endl<< "New event: " << std::endl << std::endl;
      for(auto vcont : *m_TruthParticles){
        //std::cout<<std::endl<<"New Event"<<std::endl;  
        float locp_px    =vcont->px(); float locp_py    =vcont->py();
        float locp_pz    =vcont->pz(); float locp_e     =vcont->e() ;
        TLorentzVector tlv_parent; tlv_parent.SetPxPyPzE(locp_px,locp_py,locp_pz,locp_e);
        
        int apdgid = vcont->auxdata<int>("pdgId");  
        
          
        //we will take only B+ and B- coming from H
        
          
        //521==B+(u ab) match with ab from H 
        //511==B0  
        if((apdgid == 521) && (tlv_parent.DeltaR(v_parton_bbfromH[1]) < 0.5 )){
  
            
            for(unsigned int child_index=0;child_index < vcont->nChildren(); child_index++){
                const xAOD::TruthParticle* child = vcont->child(child_index);
                if(child->pdgId()==443)
                {
                    contt++;
                    float loc_px    =child->px(); float loc_py    =child->py();
                    float loc_pz    =child->pz(); float loc_e     =child->e() ;
                    TLorentzVector tlv_B; tlv_B.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
                    
                    J_B.push_back(tlv_B);
                    B_plus.push_back(tlv_parent);
                    
                    //find e+ e- from J/psi decays
                    
                    //std::cout<<std::endl<<"J/psi(pdg id "<<child->pdgId()<<")decays from B0 (pdgid "<<vcont->pdgId()<<") :"<<std::endl;
                    
                    for(unsigned int gchild_index=0;gchild_index < child->nChildren(); gchild_index++){
                        
                        const xAOD::TruthParticle* gchild = child->child(gchild_index);
                        /*
                        std::cout<<"id: "<<gchild->auxdata<int>("pdgId")<<" ";  
                        std::cout<<"charge: "<<gchild->charge()<<" ";
                        std::cout<<"status: "<<gchild->status()<<" ";
                        std::cout<<std::endl;
                        */
                        float loc_px    =gchild->px(); float loc_py    =gchild->py();
                        float loc_pz    =gchild->pz(); float loc_e     =gchild->e() ;
                        TLorentzVector tlv_e; tlv_e.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
                        
                        if(gchild->charge() < 0) e_B.push_back(tlv_e);
                        if(gchild->charge() > 0) ae_B.push_back(tlv_e);
                          
                            
                    }
                    
                    
                }
                else if(child->pdgId()==321)
                {
                    
                    float loc_px    =child->px(); float loc_py    =child->py();
                    float loc_pz    =child->pz(); float loc_e     =child->e() ;
                    TLorentzVector tlv_B; tlv_B.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
                    
                    K_B.push_back(tlv_B);
                }
            }
        //Just to print if the decays are correctly forced   
        //Don't print need contt==1  
        if(contt==100){
            std::cout<<std::endl<<"B+ decays: ";
            for(unsigned int child_index=0;child_index < vcont->nChildren(); child_index++){
                const xAOD::TruthParticle* child = vcont->child(child_index);
                //std::cout<<child->pdgId()<<" ";
                
                }
                
        }
        
        Bindex++;           
        }
          
        //-521==B-(au b) match with b from H    
        else if ((apdgid == -521) && (tlv_parent.DeltaR(v_parton_bbfromH[0]) < 0.5 )) {
            
            
        
            for(unsigned int child_index=0;child_index < vcont->nChildren(); child_index++){
                const xAOD::TruthParticle* child = vcont->child(child_index);
                if(child->pdgId()==443)
                {
                    conttt++;
                    float loc_px    =child->px(); float loc_py    =child->py();
                    float loc_pz    =child->pz(); float loc_e     =child->e() ;
                    TLorentzVector tlv_aB; tlv_aB.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e);
                    J_aB.push_back(tlv_aB);
                    
                    B_minus.push_back(tlv_parent);
                    
                    for(unsigned int gchild_index=0;gchild_index < child->nChildren(); gchild_index++){
                        
                        const xAOD::TruthParticle* gchild = child->child(gchild_index);
                        /*
                        std::cout<<"id: "<<gchild->auxdata<int>("pdgId")<<" ";  
                        std::cout<<"charge: "<<gchild->charge()<<" ";
                        std::cout<<"status: "<<gchild->status()<<" ";
                        std::cout<<std::endl;
                        */
                        float loc_px    =gchild->px(); float loc_py    =gchild->py();
                        float loc_pz    =gchild->pz(); float loc_e     =gchild->e() ;
                        TLorentzVector tlv_e; tlv_e.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
                        
                        if(gchild->charge() < 0) e_aB.push_back(tlv_e);
                        if(gchild->charge() > 0) ae_aB.push_back(tlv_e);
                    }
                }
                else if(child->pdgId()==-321)
                {
                    
                    float loc_px    =child->px(); float loc_py    =child->py();
                    float loc_pz    =child->pz(); float loc_e     =child->e() ;
                    TLorentzVector tlv_aB; tlv_aB.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e);
                    K_aB.push_back(tlv_aB);
                }
                
        }
        //Just to print if the decays are correctly forced   
        //Don't print need conttt == 1
        if(conttt==10){
            //std::cout<<std::endl<<"B- decays: ";
            for(unsigned int child_index=0;child_index < vcont->nChildren(); child_index++){
            const xAOD::TruthParticle* child = vcont->child(child_index);
            //std::cout<<child->pdgId()<<" ";
            }
            
            
            
        }
        
        Bindex++;      
        }
      
        
      }
      
      //std::cout<<Bindex<<std::endl;
      //now we have B_plus and B_minus we looking for we verify we have one on each event
      //Only take those events with B+ and B-
      //std::cout<<std::endl;
      
      
      //
      //Boosted to X frame
      //
      if(Bindex==2 && contt>0 && conttt>0) 
      {
          //std::cout<<std::endl<<"Valid event";
          
          
          
          
          J_B_Boosted.push_back(J_B[0]); // J/psi from B0          
          J_aB_Boosted.push_back(J_aB[0]); // J/psi from aB0
          
          //TVector3 H_boost_vector = v_parton_bbfromH[2].BoostVector(); //We take the Higgs from which b ab decay
      
          J_B_Boosted[0].Boost(-H_boost_vector); //booost b 
          J_aB_Boosted[0].Boost(-H_boost_vector); //boost ab
      
          JJ_Boosted_dR = J_B_Boosted[0].DeltaR(J_aB_Boosted[0]); //computation of angular differnce
          cos_JJ_Boosted_dR = TMath::Cos(JJ_Boosted_dR); 
          
          JJ_Boosted_dPhi = J_B_Boosted[0].DeltaPhi(J_aB_Boosted[0]);
          cos_JJ_Boosted_dPhi = TMath::Cos(JJ_Boosted_dPhi);
          
          JJ_Boosted_dEta = J_B_Boosted[0].PseudoRapidity()-J_aB_Boosted[0].PseudoRapidity();
          cos_JJ_Boosted_dEta = TMath::Cos(JJ_Boosted_dEta);
          
          
          //---------------electrons-------------//
          
          e_ae_dR = e_B[0].DeltaR(ae_aB[0]); //computation of angular differnce   
          ae_e_dR = ae_B[0].DeltaR(e_aB[0]); //computation of angular differnce
          
          e_ae_dPhi = e_B[0].DeltaPhi(ae_aB[0]);
          ae_e_dPhi = ae_B[0].DeltaPhi(e_aB[0]);
             
          e_ae_dEta = e_B[0].PseudoRapidity()-ae_aB[0].PseudoRapidity();
          ae_e_dEta = ae_B[0].PseudoRapidity()-e_aB[0].PseudoRapidity();
          /*
          std::cout<<e_B[0].PseudoRapidity()<<" "<<ae_aB[0].PseudoRapidity()<<std::endl;
          std::cout<<ae_B[0].PseudoRapidity()<<" "<<e_aB[0].PseudoRapidity()<<std::endl;
          std::cout<<"new event"<<std::endl;
          
          std::cout<<J_B[0].Pseudorapidity()<<" "<<J_aB[0].Pseudorapidity()<<std::endl;
          std::cout<<J_B[0].Pseudorapidity()<<" "<<J_aB[0].Pseudorapidity()<<std::endl;
          */
          
          
          
          //Non Boosting
          
          
          //
          e_B_Boosted.push_back(e_B[0]); // e from B          
          e_aB_Boosted.push_back(e_aB[0]); // e from aB
          
          ae_B_Boosted.push_back(ae_B[0]); // e+ from B          
          ae_aB_Boosted.push_back(ae_aB[0]); // e+ from aB
        
          //TVector3 H_boost_vector = v_parton_bbfromH[2].BoostVector(); //We take the Higgs from which b ab decay
      
          e_B_Boosted[0].Boost(-H_boost_vector); //booost e from B 
          e_aB_Boosted[0].Boost(-H_boost_vector); //boost e from aB
          
          ae_B_Boosted[0].Boost(-H_boost_vector); //booost e+ from B 
          ae_aB_Boosted[0].Boost(-H_boost_vector); //boost e+ from aB
      
          e_ae_Boosted_dR = e_B_Boosted[0].DeltaR(ae_aB_Boosted[0]); //computation of angular differnce
          cos_e_ae_Boosted_dR = TMath::Cos(e_ae_Boosted_dR); 
          
          ae_e_Boosted_dR = ae_B_Boosted[0].DeltaR(e_aB_Boosted[0]); //computation of angular differnce
          cos_ae_e_Boosted_dR = TMath::Cos(ae_e_Boosted_dR); 
          
          
          e_ae_Boosted_dPhi = e_B_Boosted[0].DeltaPhi(ae_aB_Boosted[0]);
          cos_e_ae_Boosted_dPhi = TMath::Cos(e_ae_Boosted_dPhi);
          
          ae_e_Boosted_dPhi = ae_B_Boosted[0].DeltaPhi(e_aB_Boosted[0]);
          cos_ae_e_Boosted_dPhi = TMath::Cos(ae_e_Boosted_dPhi);
          
          
          e_ae_Boosted_dEta = e_B_Boosted[0].PseudoRapidity()-ae_aB_Boosted[0].PseudoRapidity();
          cos_e_ae_Boosted_dEta = TMath::Cos(e_ae_Boosted_dEta);
          
          ae_e_Boosted_dEta = ae_B_Boosted[0].PseudoRapidity()-e_aB_Boosted[0].PseudoRapidity();
          cos_ae_e_Boosted_dEta = TMath::Cos(ae_e_Boosted_dEta);
          
          
          
          
          //std::cout<<cos_JJ_boosted_dR<<" ";
      
      
      // Boosted to B+ B- frame each
      
      
          J_B_boost.push_back(J_B[0]); // J/psi from B+
          J_aB_boost.push_back(J_aB[0]); // J/psi from B-
          
          TVector3 B_plus_boost_vector = B_plus[0].BoostVector(); 
          TVector3 B_minus_boost_vector = B_minus[0].BoostVector();
          
          J_B_boost[0].Boost(-B_plus_boost_vector);
          J_aB_boost[0].Boost(-B_minus_boost_vector);
          
          JJ_boost_dR = J_B_boost[0].DeltaR(J_aB_boost[0]);
          cos_JJ_boost_dR = TMath::Cos(JJ_boost_dR); 
          
          JJ_boost_dPhi = J_B_boost[0].DeltaPhi(J_aB_boost[0]);
          cos_JJ_boost_dPhi = TMath::Cos(JJ_boost_dPhi);
          
          JJ_boost_dEta = J_B_boost[0].PseudoRapidity()-J_aB_boost[0].PseudoRapidity();
          cos_JJ_boost_dEta = TMath::Cos(JJ_boost_dEta);
          //std::cout<<cos_JJ_boost_dR<<" ";
          
          //----------electrons--------------//
          
          e_B_boost.push_back(e_B[0]); // e from B+
          e_aB_boost.push_back(e_aB[0]); // e from B-
          
          ae_B_boost.push_back(e_B[0]); // e+ from B+
          ae_aB_boost.push_back(e_aB[0]); // e+ from B-
          
          //TVector3 B_plus_boost_vector = B_plus[0].BoostVector(); 
          //TVector3 B_minus_boost_vector = B_minus[0].BoostVector();
          
          e_B_boost[0].Boost(-B_plus_boost_vector);
          e_aB_boost[0].Boost(-B_minus_boost_vector);
          
          ae_B_boost[0].Boost(-B_plus_boost_vector);
          ae_aB_boost[0].Boost(-B_minus_boost_vector);
          
          
          e_ae_boost_dR = e_B_boost[0].DeltaR(ae_aB_boost[0]); //computation of angular differnce
          cos_e_ae_boost_dR = TMath::Cos(e_ae_boost_dR); 
          
          ae_e_boost_dR = ae_B_boost[0].DeltaR(e_aB_boost[0]); //computation of angular differnce
          cos_ae_e_boost_dR = TMath::Cos(ae_e_boost_dR); 
          
          e_ae_boost_dPhi = e_B_boost[0].DeltaPhi(ae_aB_boost[0]);
          cos_e_ae_boost_dPhi = TMath::Cos(e_ae_boost_dPhi);
          
          ae_e_boost_dPhi = ae_B_boost[0].DeltaPhi(e_aB_boost[0]);
          cos_ae_e_boost_dPhi = TMath::Cos(ae_e_boost_dPhi);
          
          
          e_ae_boost_dEta = e_B_boost[0].PseudoRapidity()-ae_aB_boost[0].PseudoRapidity();
          cos_e_ae_boost_dEta = TMath::Cos(e_ae_boost_dEta);
          
          ae_e_boost_dEta = ae_B_boost[0].PseudoRapidity()-e_aB_boost[0].PseudoRapidity();
          cos_ae_e_boost_dEta = TMath::Cos(ae_e_boost_dEta);
          
          //Angle between the vectors
          
          Phi = ae_B_boost[0].Angle(e_aB_boost[0].Vect());
          Cos_Phi = TMath::Cos(Phi);
          
          
      }
      
      
      //if(contt>0||conttt>0) std::cout<<std::endl<<std::endl<<"new event: ";
      
      //TAU TAU CP study https://arxiv.org/pdf/1308.2674.pdf
      
      //J_B -> J/psi coming from B+
      //J_ab -> J/psi coming from B-
      //K_B -> K+ coming B+
      //K_aB -> K- coming B-
      //We need to use the charged products
      if(Bindex==2 && contt>0 && conttt>0) 
      {
          //non boosted 3vectors
          TVector3 unit_Jplus = (J_B[0].BoostVector()).Unit();
          TVector3 unit_Jminus = (J_aB[0].BoostVector()).Unit();
          
          //Zero-Momentum Frame
          TVector3 ZMF = (J_B[0]+J_aB[0]).BoostVector();
          
          TLorentzVector J_B_boosted = J_B[0]; 
          TLorentzVector J_aB_boosted = J_aB[0]; 
          
          J_B_boosted.Boost(-ZMF);
          J_aB_boosted.Boost(-ZMF);
          
          //Boosted 3vectors
          TVector3 unit_Jplus_boosted = (J_B_boosted.BoostVector()).Unit();
          TVector3 unit_Jminus_boosted = (J_B_boosted.BoostVector()).Unit();
          
          //Now we need to compute the parallel and perpendicular components of the boosted in respect to the non boosted
          
          //Parallel Component
          float parallel_Kplus_boosted = (unit_Jplus_boosted * unit_Jplus);
          float parallel_Kminus_boosted = (unit_Jminus_boosted * unit_Jminus);
          
          //Perpendicular = (vector - parallel)
          TVector3 orthogonal_Kplus_boosted = (unit_Jplus_boosted - unit_Jplus).Unit();
          TVector3 orthogonal_Kminus_boosted = (unit_Jminus_boosted - unit_Jminus).Unit();
          
          //Compute Phi
          Phi_plane_angle = TMath::ACos(orthogonal_Kplus_boosted.Dot(orthogonal_Kminus_boosted));
          
          //Compute theta
          Theta_plane_angle = unit_Jminus_boosted.Dot(orthogonal_Kplus_boosted.Cross(orthogonal_Kminus_boosted));
          
      }
      
      
      
      
          
          
      
          
          
            
      
   
       } 
      

  else if(proc==PROCESS::LOWMASSHPLUS)
    {
      v_parton_Hplus_decays  = m_LowMassHplusSelector->DecayH(m_TruthParticles,(float)0.,(float)10., m_LowMassHplusSelector->DECAY_H::ALL);
      v_H_decay_pdgid        = m_LowMassHplusSelector->GetHiggsDecays(m_TruthParticles,37);
      
      //for spin correlation observables

      v_TopDecayChain                 = m_LowMassHplusSelector->GetTopDecayChain(m_TruthParticles);
      v_parton_topdecaychain_lepton   = m_LowMassHplusSelector->GetParticleInTopDecays(v_TopDecayChain, v_leptons_pdgid  );
      v_parton_topdecaychain_neutrino = m_LowMassHplusSelector->GetParticleInTopDecays(v_TopDecayChain, v_neutrinos_pdgid);

      v_angana_lep_top_dphi	  =m_LowMassHplusSelector->CalculateAngle(v_TopDecayChain, m_LowMassHplusSelector->ANGLE::LEP_TOP_DPHI);      
      v_angana_lep_top_dr	  =m_LowMassHplusSelector->CalculateAngle(v_TopDecayChain, m_LowMassHplusSelector->ANGLE::LEP_TOP_DR);        
      v_angana_nu_top_dphi	  =m_LowMassHplusSelector->CalculateAngle(v_TopDecayChain, m_LowMassHplusSelector->ANGLE::NU_TOP_DPHI);       
      v_angana_nu_top_dr	  =m_LowMassHplusSelector->CalculateAngle(v_TopDecayChain, m_LowMassHplusSelector->ANGLE::NU_TOP_DR);         
      v_angana_lep_top_dphi_topcm =m_LowMassHplusSelector->CalculateAngle(v_TopDecayChain, m_LowMassHplusSelector->ANGLE::LEP_TOP_DPHI_TOPCM);
      v_angana_nu_top_dphi_topcm  =m_LowMassHplusSelector->CalculateAngle(v_TopDecayChain, m_LowMassHplusSelector->ANGLE::NU_TOP_DPHI_TOPCM); 
    }

  else if(proc==PROCESS::ZB)
    {
      if     (v_muons.size()==2)      v_Zb_Z = m_ZbSelector->GetZboson(v_muons);
      else if(v_electrons.size()==2)  v_Zb_Z = m_ZbSelector->GetZboson(v_electrons);
      //cut on abs pdgid
      std::vector<TLorentzVector> v_b = m_GenericSelector->GenericParticle(m_TruthParticles,(float)0.,(float)10.,
									   m_GenericSelector->CUT_PDG::ABS,5,
									   m_GenericSelector->CUT_STATUS::YES,23);
      
      if(v_b.size()!=0)v_Zb_wzjets=m_ZbSelector->GetBJets(v_b,v_AntiKt4TruthWZJets);
    }

  return EL::StatusCode::SUCCESS;

}


EL::StatusCode MyxAODTruthAnalysis :: NtupleFiller ()
{
  if(v_neutrinos         .size()!=0)m_NtupleManager->FillParticleCollection(v_neutrinos,         "Neutrinos"         );
  if(v_electrons         .size()!=0)m_NtupleManager->FillParticleCollection(v_electrons,         "Electrons"         );
  if(v_muons             .size()!=0)m_NtupleManager->FillParticleCollection(v_muons,             "Muons"             );
  if(v_taus              .size()!=0)m_NtupleManager->FillParticleCollection(v_taus,              "Taus"              );
  if(v_AntiKt4TruthWZJets.size()!=0)m_NtupleManager->FillParticleCollection(v_AntiKt4TruthWZJets,"AntiKt4TruthWZJets");
  if(v_particle_wboson   .size()!=0)m_NtupleManager->FillParticleCollection(v_particle_wboson,   "ParticleWBoson"    );

  if(v_parton_u          .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_u,          "UpQuarks"          );
  if(v_parton_au         .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_au,         "Anti-UpQuarks"     );
  if(v_parton_d          .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_d,          "DownQuarks"        );
  if(v_parton_ad         .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_ad,         "Anti-DownQuarks"   );
  if(v_parton_c          .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_c,          "CharmQuarks"       );
  if(v_parton_ac         .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_ac,         "Anti-CharmQuarks"  );
  if(v_parton_s          .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_s,          "StrangeQuarks"     );
  if(v_parton_as         .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_as,         "Anti-StrangeQuarks");
  if(v_parton_top        .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_top,        "TopQuarks"         );
  if(v_parton_atop       .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_atop,       "Anti-TopQuarks"    );
  if(v_parton_b          .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_b,          "BottomQuarks"      );
  if(v_parton_ab         .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_ab,         "Anti-BottomQuarks" );
  if(v_parton_wplus      .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_wplus,      "WPlus"             );
  if(v_parton_wminus     .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_wminus,     "WMinus"            );
  if(v_parton_higgs      .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_higgs,      "Higgs"             );
    //cosas del decay
  if(B_plus      .size()!=0)m_NtupleManager->FillParticleCollection(B_plus,      "B_plus"             ); 
  if(B_minus      .size()!=0)m_NtupleManager->FillParticleCollection(B_minus,      "B_minus"             );   
   
  if(J_B      .size()!=0)m_NtupleManager->FillParticleCollection(J_B,      "J_B"             );   
  if(J_aB      .size()!=0)m_NtupleManager->FillParticleCollection(J_aB,      "J_aB"             ); 
  
  if(J_B_Boosted      .size()!=0)m_NtupleManager->FillParticleCollection(J_B_Boosted,      "J_B_Boosted"             );   
  if(J_aB_Boosted      .size()!=0)m_NtupleManager->FillParticleCollection(J_aB_Boosted,      "J_aB_Boosted"             ); 
  
  if(J_B_boost      .size()!=0)m_NtupleManager->FillParticleCollection(J_B_boost,      "J_B_boost"             );  
  if(J_aB_boost      .size()!=0)m_NtupleManager->FillParticleCollection(J_aB_boost,      "J_aB_boost"             );
    
  if(K_B      .size()!=0)m_NtupleManager->FillParticleCollection(K_B,      "K_B"             );   
  if(K_aB      .size()!=0)m_NtupleManager->FillParticleCollection(K_aB,      "K_aB"             );   

  //tttt/bbtt specific
  if(v_parton_ttfromH     .size()!=0)m_NtupleManager->FillHCollection(v_parton_ttfromH,           "FromH"            );
  if(v_parton_ttnotfromH  .size()!=0)m_NtupleManager->FillHCollection(v_parton_ttnotfromH,        "TTNotFromH"       );
  //tcbb specific
  if(v_parton_bbfromH     .size()!=0)m_NtupleManager->FillHCollection(v_parton_bbfromH,           "FromH"            );
  if(v_parton_tcnotfromH  .size()!=0)m_NtupleManager->FillHCollection(v_parton_tcnotfromH,        "TCNotFromH"       );
  //tXq specific
  if(v_parton_sisterH     .pdgid!=-99)m_NtupleManager->FillHCollection(v_parton_sisterH,           "SisH"             );
  if(v_parton_fromothertop.pdgid!=-99)m_NtupleManager->FillHCollection(v_parton_fromothertop,      "FromOtherTop"     );

  //Hplus
  if(v_parton_Hplus_decays.size()!=0)m_NtupleManager->FillHCollection(v_parton_Hplus_decays,      "FromH"            );
  //B+ Decays
  if(v_Bplus_hadrons_decay.size()!=0)m_NtupleManager->FillHCollection(v_Bplus_hadrons_decay,      "FromB"            );  

  if(v_parton_topdecaychain_lepton  .size()!=0)m_NtupleManager->FillParticleCollection(v_parton_topdecaychain_lepton  , "Topdecaychain_lepton"   );
  if(v_parton_topdecaychain_neutrino.size()!=0)m_NtupleManager->FillParticleCollection(v_parton_topdecaychain_neutrino, "Topdecaychain_neutrino" );
 
  if(v_HT !=0)m_NtupleManager->FillScalarCollection(v_HT, "HT" );

  if(v_angana_lep_top_dphi	  !=0)m_NtupleManager->FillScalarCollection(v_angana_lep_top_dphi,       "Angana_lep_top_dphi"         );
  if(v_angana_lep_top_dr	  !=0)m_NtupleManager->FillScalarCollection(v_angana_lep_top_dr,         "Angana_lep_top_dr"           );
  if(v_angana_nu_top_dphi	  !=0)m_NtupleManager->FillScalarCollection(v_angana_nu_top_dphi,        "Angana_nu_top_dphi"          );
  if(v_angana_nu_top_dr	          !=0)m_NtupleManager->FillScalarCollection(v_angana_nu_top_dr,          "Angana_nu_top_dr"            );
  if(v_angana_lep_top_dphi_topcm  !=0)m_NtupleManager->FillScalarCollection(v_angana_lep_top_dphi_topcm, "Angana_lep_top_dphi_topcm"   );
  if(v_angana_nu_top_dphi_topcm   !=0)m_NtupleManager->FillScalarCollection(v_angana_nu_top_dphi_topcm,  "Angana_nu_top_dphi_topcm"    );

  m_NtupleManager->FillScalarCollection(nmatched,          "tXq_nmatched");
  m_NtupleManager->FillScalarCollection(MultipleMatching, "tXq_MultipleMatching");
  m_NtupleManager->FillScalarCollection(b1_JetIndex,       "tXq_b1Index");
  m_NtupleManager->FillScalarCollection(b2_JetIndex,       "tXq_b2Index");
  m_NtupleManager->FillScalarCollection(b3_JetIndex,      "tXq_b3Index");
  m_NtupleManager->FillScalarCollection(qS_JetIndex,      "tXq_qIndex");
  m_NtupleManager->FillScalarCollection(b1b2dR,            "b1b2dR");
  m_NtupleManager->FillScalarCollection(b1b2_boosted_dR,            "b1b2_boosted_dR");  
    
  m_NtupleManager->FillScalarCollection(JJ_Boosted_dR,            "JJ_Boosted_dR"); 
  m_NtupleManager->FillScalarCollection(JJ_boost_dR,            "JJ_boost_dR"); 
  m_NtupleManager->FillScalarCollection(cos_JJ_Boosted_dR,            "cos_JJ_Boosted_dR"); 
  m_NtupleManager->FillScalarCollection(cos_JJ_boost_dR,            "cos_JJ_boost_dR"); 
    
  m_NtupleManager->FillScalarCollection(JJ_Boosted_dPhi,            "JJ_Boosted_dPhi"); 
  m_NtupleManager->FillScalarCollection(JJ_boost_dPhi,            "JJ_boost_dPhi"); 
  m_NtupleManager->FillScalarCollection(cos_JJ_Boosted_dPhi,            "cos_JJ_Boosted_dPhi"); 
  m_NtupleManager->FillScalarCollection(cos_JJ_boost_dPhi,            "cos_JJ_boost_dPhi"); 
  
  m_NtupleManager->FillScalarCollection(JJ_Boosted_dEta,            "JJ_Boosted_dEta"); 
  m_NtupleManager->FillScalarCollection(JJ_boost_dEta,            "JJ_boost_dEta"); 
  m_NtupleManager->FillScalarCollection(cos_JJ_Boosted_dEta,            "cos_JJ_Boosted_dEta"); 
  m_NtupleManager->FillScalarCollection(cos_JJ_boost_dEta,            "cos_JJ_boost_dEta");  
    
  m_NtupleManager->FillScalarCollection(Phi_plane_angle  ,            "Phi_plane_angle"); 
    
  m_NtupleManager->FillScalarCollection(e_ae_Boosted_dR,            "e_ae_Boosted_dR"); 
  m_NtupleManager->FillScalarCollection(e_ae_boost_dR,            "e_ae_boost_dR"); 
  m_NtupleManager->FillScalarCollection(ae_e_Boosted_dR,            "ae_e_Boosted_dR"); 
  m_NtupleManager->FillScalarCollection(ae_e_boost_dR,            "ae_e_boost_dR");   
    
  m_NtupleManager->FillScalarCollection(e_ae_dR,            "e_ae_dR"); 
  m_NtupleManager->FillScalarCollection(ae_e_dR,            "ae_e_dR");   
   
  m_NtupleManager->FillScalarCollection(cos_e_ae_Boosted_dR,            "cos_e_ae_Boosted_dR"); 
  m_NtupleManager->FillScalarCollection(cos_e_ae_boost_dR,            "cos_e_ae_boost_dR"); 
  m_NtupleManager->FillScalarCollection(cos_ae_e_Boosted_dR,            "cos_ae_e_Boosted_dR"); 
  m_NtupleManager->FillScalarCollection(cos_ae_e_boost_dR,            "cos_ae_e_boost_dR");
  
  m_NtupleManager->FillScalarCollection(e_ae_Boosted_dPhi,            "e_ae_Boosted_dPhi"); 
  m_NtupleManager->FillScalarCollection(e_ae_boost_dPhi,            "e_ae_boost_dPhi"); 
  m_NtupleManager->FillScalarCollection(ae_e_Boosted_dPhi,            "ae_e_Boosted_dPhi"); 
  m_NtupleManager->FillScalarCollection(ae_e_boost_dPhi,            "ae_e_boost_dPhi");   
    
  m_NtupleManager->FillScalarCollection(e_ae_dPhi,            "e_ae_dPhi"); 
  m_NtupleManager->FillScalarCollection(ae_e_dPhi,            "ae_e_dPhi");   
   
  m_NtupleManager->FillScalarCollection(cos_e_ae_Boosted_dPhi,            "cos_e_ae_Boosted_dPhi"); 
  m_NtupleManager->FillScalarCollection(cos_e_ae_boost_dPhi,            "cos_e_ae_boost_dPhi"); 
  m_NtupleManager->FillScalarCollection(cos_ae_e_Boosted_dPhi,            "cos_ae_e_Boosted_dPhi"); 
  m_NtupleManager->FillScalarCollection(cos_ae_e_boost_dPhi,            "cos_ae_e_boost_dPhi");
    
  m_NtupleManager->FillScalarCollection(e_ae_Boosted_dEta,            "e_ae_Boosted_dEta"); 
  m_NtupleManager->FillScalarCollection(e_ae_boost_dEta,            "e_ae_boost_dEta"); 
  m_NtupleManager->FillScalarCollection(ae_e_Boosted_dEta,            "ae_e_Boosted_dEta"); 
  m_NtupleManager->FillScalarCollection(ae_e_boost_dEta,            "ae_e_boost_dEta"); 
    
  m_NtupleManager->FillScalarCollection(e_ae_dEta,            "e_ae_dEta"); 
  m_NtupleManager->FillScalarCollection(ae_e_dEta,            "ae_e_dEta");  
   
  m_NtupleManager->FillScalarCollection(cos_e_ae_Boosted_dEta,            "cos_e_ae_Boosted_dEta"); 
  m_NtupleManager->FillScalarCollection(cos_e_ae_boost_dEta,            "cos_e_ae_boost_dEta"); 
  m_NtupleManager->FillScalarCollection(cos_ae_e_Boosted_dEta,            "cos_ae_e_Boosted_dEta"); 
  m_NtupleManager->FillScalarCollection(cos_ae_e_boost_dEta,            "cos_ae_e_boost_dEta");
    
  m_NtupleManager->FillScalarCollection(Phi,            "Phi");
  m_NtupleManager->FillScalarCollection(Cos_Phi,            "Cos_Phi");   
    
  m_NtupleManager->FillScalarCollection(b1b3dR,            "b1b3dR");
  m_NtupleManager->FillScalarCollection(b2b3dR,            "b2b3dR");
  m_NtupleManager->FillScalarCollection(b1qSdR,            "b1qSdR");
  m_NtupleManager->FillScalarCollection(b1qSdR,            "b1qSdR");
  m_NtupleManager->FillScalarCollection(minDRbb,           "minDRbb");
  m_NtupleManager->FillScalarCollection(dR_Jpsi_K,           "dR_Jpsi_K");  
   
    

  if(v_parton_bbfromH.size()!=0 && v_parton_Hplus_decays.size()!=0){
    std::cout<<"Got both a neutral and a charged Higgs, exit now ! "<<std::endl;
    exit(0);
  }

  //Zb specific
  if(v_Zb_Z              .size()!=0){
    std::cout<<"mhm 1"<<std::endl;
    m_NtupleManager->FillParticleCollection(v_Zb_Z,              "Zb_Z"              );
  }
  
  if(v_Zb_wzjets         .size()!=0){
    std::cout<<"mhm 2"<<std::endl;    
    m_NtupleManager->FillParticleCollection(v_Zb_wzjets,         "Zb_wzjets"         );
  }
  //to be added
  //m_NtupleManager->FillHCollection(v_parton_bbnotfromH,        "BBNotFromH"        );

  //temporary implementation
  m_NtupleManager->FillEventCollection(mc_event_weight);
  m_NtupleManager->FillEventCollectionHack(v_H_decay_pdgid);

  return EL::StatusCode::SUCCESS;
}


EL::StatusCode MyxAODTruthAnalysis :: ClearParticles ()
{
  //particle level 
  v_neutrinos.clear(); 
  v_muons.clear(); 
  v_electrons.clear();
  v_taus.clear(); 
  v_AntiKt4TruthWZJets.clear();
  v_particle_wboson.clear();
  //parton level
  v_parton_u.clear(); 
  v_parton_au.clear(); 
  v_parton_d.clear(); 
  v_parton_ad.clear(); 
  v_parton_c.clear(); 
  v_parton_ac.clear(); 
  v_parton_s.clear(); 
  v_parton_as.clear(); 
  v_parton_top.clear(); 
  v_parton_atop.clear(); 
  v_parton_b.clear(); 
  v_parton_ab.clear(); 
  v_parton_wplus.clear();  
  v_parton_wminus.clear();  
  v_parton_higgs.clear(); 
    
  B_plus.clear();
  B_minus.clear();
    
  J_B.clear();
  J_aB.clear();
    
  K_B.clear();
  K_aB.clear(); 
    
  J_B_boost.clear();
  J_aB_boost.clear(); 
    
  J_B_Boosted.clear();
  J_aB_Boosted.clear(); 
    
  e_B_boost.clear();
  e_aB_boost.clear(); 
    
  ae_B_Boosted.clear();
  ae_aB_Boosted.clear(); 
    
  ae_B_boost.clear();
  ae_aB_boost.clear(); 
    
  e_B_Boosted.clear();
  e_aB_Boosted.clear();
    
  e_aB.clear(); 
  ae_aB.clear();  
    
  e_B.clear(); 
  ae_B.clear();    
    
    
    
    
  //process specific
  v_parton_ttfromH.clear(); 
  v_parton_ttnotfromH.clear();
  v_parton_bbfromH.clear(); 
  v_parton_tcnotfromH.clear(); 
  v_parton_Hplus_decays.clear(); 
  v_parton_fromH.clear(); 
  v_Zb_Z.clear();
  v_Zb_wzjets.clear();
  v_H_decay_pdgid.clear();
  v_parton_topdecaychain_lepton.clear();
  v_parton_topdecaychain_neutrino.clear();
  
  v_TopDecayChain.clear();

  mc_event_weight=0;
  v_HT=0.;
  v_angana_lep_top_dphi=0.;
  v_angana_lep_top_dr=0.;
  v_angana_nu_top_dphi=0.;
  v_angana_nu_top_dr=0.;
  v_angana_lep_top_dphi_topcm=0.;
  v_angana_nu_top_dphi_topcm=0.;

  return EL::StatusCode::SUCCESS;

}


bool MyxAODTruthAnalysis :: CutFlow ()
{
  float w = mc_event_weight;

  m_CutFlowTools->cutFlow("CutFlow","All" ,1.);
  m_CutFlowTools->cutFlow("CutFlow weighted","All" ,w);

  if(deriv==DERIVATION::TRUTH1)
    {

      if(v_muons.size()>=1)
	{
	  m_CutFlowTools->cutFlow("CutFlow","OneMu" ,1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted","OneMu" ,w);
	}
  
      if(v_electrons.size()>=1)
	{
	  m_CutFlowTools->cutFlow("CutFlow","OneEl" ,1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted","OneEl" ,w );
	}
  
      if(v_taus.size()>=1)
	{
	  m_CutFlowTools->cutFlow("CutFlow","OneTau",1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted","OneTau",w );
	}

      if(v_neutrinos.size()>=1)
	{
	  m_CutFlowTools->cutFlow("CutFlow","OneNu" ,1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted","OneNu" ,w );
	}
  
      if(! (v_muons.size()>=1 || v_electrons.size()>=1 || v_taus.size()>=1) ) return false;
      if(! (v_neutrinos.size()>=1) ) return false;
  
      m_CutFlowTools->cutFlow("CutFlow","Lep+Nu",1.);
      m_CutFlowTools->cutFlow("CutFlow weighted","Lep+Nu",w );
 
      if(v_AntiKt4TruthWZJets.size()==4)
	{
	  m_CutFlowTools->cutFlow("CutFlow","FourJets" ,1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted","FourJets" ,w );
	}
  
      if(v_AntiKt4TruthWZJets.size()==5)
	{
	  m_CutFlowTools->cutFlow("CutFlow","FiveJets" ,1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted","FiveJets" ,w );
	}

      if(v_AntiKt4TruthWZJets.size()>=6)
	{
	  m_CutFlowTools->cutFlow("CutFlow",">=SixJets",1.);
	  m_CutFlowTools->cutFlow("CutFlow weighted",">=SixJets",w );
	}

      if(v_AntiKt4TruthWZJets.size()>=4) return true;
      else return false;

    }
 
  return true;
}
