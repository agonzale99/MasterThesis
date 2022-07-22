#include "MyTruthAnalysis/TruthAna_GenericSelector.h"
//For std::pair
#include <utility>

TruthAna_GenericSelector::TruthAna_GenericSelector()
{
}


TruthAna_GenericSelector::~TruthAna_GenericSelector()
{
}

//Meant to be used for getting all final state partons
std::vector< std::pair< TLorentzVector,int > > TruthAna_GenericSelector::GetAllPartons(const xAOD::TruthParticleContainer *cont, CUT_STATUS cut_status, int status)
{
  std::vector< std::pair< TLorentzVector , int > > generic_stable_parton;
    
    for(auto vcont : *cont){
      int par_pdgid=vcont->auxdata<int>("pdgId");
      float loc_px    =vcont->px(); float loc_py    =vcont->py();
      float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
      TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 

      if( cut_status==CUT_STATUS::YES && ( ! (vcont->auxdata<int>("status")== status ) ) ) continue;
    
      //generic_stable_parton.first.push_back(tlv);                                                                   
      //generic_stable_parton.second.push_back(par_pdgid);                                                                       
      //if((par_pdgid == 443 || par_pdgid == 521) && cut_status==CUT_STATUS::NO) std::cout << par_pdgid << " " << vcont->auxdata<int>("status") << std::endl;
        
      //
      //if(par_pdgid == 443) std::cout << vcont -> nChildren() << std::endl; 
        
      generic_stable_parton.push_back( std::make_pair(tlv,par_pdgid) );

    }

  return generic_stable_parton;
}

std::vector<TLorentzVector>  TruthAna_GenericSelector::RetrieveParton(std::vector< std::pair< TLorentzVector,int > > generic_stable_parton, CUT_PDG cut_pdg, int pdgid)
{
  std::vector<TLorentzVector> tlv_vec;
  
  for(auto auto_generic_stable_parton : generic_stable_parton ){
    int par_pdgid=auto_generic_stable_parton.second;
    TLorentzVector tlv=auto_generic_stable_parton.first;
    
    if     ( cut_pdg==CUT_PDG::YES && ( ! (par_pdgid== pdgid ) ) ) continue;
    else if( cut_pdg==CUT_PDG::ABS && ( ! (fabs(par_pdgid)== pdgid ) ) ) continue;
    //else { std::cout<<"No PDGID cut is applied, please check, aborting now"<<std::endl; exit(0); }

      
    tlv_vec.push_back(tlv);
  }
  
  return tlv_vec;
}

std::vector<TLorentzVector> TruthAna_GenericSelector::GenericParticle(const xAOD::TruthParticleContainer *cont, float pt, float eta, CUT_PDG cut_pdg, int pdgid, CUT_STATUS cut_status, int status)
{
  std::vector<TLorentzVector> tlv_vec;
    //std::cout << "hola" << std::endl;
  for(auto vcont : *cont){
      
    int par_pdgid=vcont->auxdata<int>("pdgId");  
    float loc_px    =vcont->pt(); float loc_py    =vcont->py();
    float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
    TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 

    
    
    if( tlv.Pt()<pt )continue;
    if( fabs( tlv.Rapidity() )>eta )continue;
    
    //special treatment for the H/A bosons
    if( par_pdgid==35 || par_pdgid== 36 || fabs(par_pdgid== 37) || par_pdgid==25 ) par_pdgid=25;
 
    if     ( cut_pdg==CUT_PDG::YES && ( ! (par_pdgid== pdgid ) ) ) continue;
    else if( cut_pdg==CUT_PDG::ABS && ( ! (fabs(par_pdgid)== pdgid ) ) ) continue;
    
    if( cut_status==CUT_STATUS::YES && ( ! (vcont->auxdata<int>("status")== status ) ) ) continue;
    //if(par_pdgid == 443 || par_pdgid == 521) std::cout << par_pdgid << " " << vcont->auxdata<int>("status") << std::endl;
      
    tlv_vec.push_back(tlv);                                                                   
  }
  
  return tlv_vec;
}

std::vector<TLorentzVector> TruthAna_GenericSelector::GenericJet(const xAOD::JetContainer *cont, float pt, float eta)
{
  std::vector<TLorentzVector> tlv_vec;
  
  for(auto vcont : *cont){
    float loc_px    =vcont->px(); float loc_py    =vcont->py();
    float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
    TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e);
    //std::cout << vcont->numConstituents() << " " <<  tlv.Pt() << " " << tlv.Eta() << " " << tlv.Phi() << std::endl;
    //for(auto stuff : vcont->getConstituents()){
    //  std::cout << stuff->pt() << " " << stuff->eta() << " " << stuff->phi() << std::endl;
    //}
    if( tlv.Pt()<pt )continue;
    if( fabs( tlv.Rapidity() )>eta )continue;
    tlv_vec.push_back(tlv);                                                                   
  }
  return tlv_vec;
}

float TruthAna_GenericSelector::CalculateHT(std::vector<TLorentzVector> cont_jet, std::vector<TLorentzVector> cont_electron, std::vector<TLorentzVector> cont_muon, std::vector<TLorentzVector> cont_tau, std::vector<TLorentzVector> cont_neutrino)
{
  float loc_HT=0.;
  
  for (auto jet : cont_jet)           { loc_HT += jet.Pt(); }
  for (auto electron : cont_electron) { loc_HT += electron.Pt(); }
  for (auto muon : cont_muon)         { loc_HT += muon.Pt(); }
  for (auto tau : cont_tau)           { loc_HT += tau.Pt(); }
  for (auto neutrino : cont_jet)      { loc_HT += neutrino.Pt(); }
  
  return loc_HT;
}



