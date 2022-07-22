#include "MyTruthAnalysis/TruthAna_LowMassHplusSelector.h"
#include <TLorentzVector.h>


//Does a loop looking for a top quark, then it looks for a W boson, then it goes for the neutrino and leptons

std::vector< std::pair< TLorentzVector,int > > TruthAna_LowMassHplusSelector::GetTopDecayChain(const xAOD::TruthParticleContainer *cont) 
{
  std::vector< std::pair< TLorentzVector,int > > particle_list;

  TLorentzVector tlv_atop;       
  TLorentzVector tlv_Wboson;   
  TLorentzVector tlv_bquark;   
  TLorentzVector tlv_lepton;   
  TLorentzVector tlv_aneutrino;

  int pdgid_atop=-100;     
  int pdgid_Wboson=-100;   
  int pdgid_lepton=-100;   
  int pdgid_aneutrino=-100;

  //Loop over all particles 
  for(auto vtops : *cont){
    
    //look for an anti-top because the process generated is H+ so the top will decay in anti-bottom and H+, while the atop will decay as SM like particle 
    bool is_atop = bool(vtops->auxdata<int>("pdgId") == -6);
    if(! is_atop ) continue;

    //Getting the top to be used for the loop
    const xAOD::TruthParticle* top_particle =  vtops;
    pdgid_atop=top_particle->auxdata<int>("pdgId");

  LOOP_TOP_CHILD: 
    for(unsigned int top_child_index=0; top_child_index < top_particle->nChildren(); top_child_index++){
      
      const xAOD::TruthParticle* top_child = top_particle->child(top_child_index);
      pdgid_Wboson=top_child->auxdata<int>("pdgId");
      //std::cout<<"Found a top quark child with pdgid"<<std::endl;
      //Excluding top which undergo t->tg      
      if( fabs(pdgid_Wboson)==6 ){
	top_particle = top_child;
	goto LOOP_TOP_CHILD;
      }
      
      if( fabs(pdgid_Wboson)!=24 )continue;
      
    LOOP_W_BOSON:
      for(unsigned int W_child_index=0; W_child_index < top_child->nChildren(); W_child_index++){

	//std::cout<<"Looping on W decay products"<<std::endl;
	const xAOD::TruthParticle* W_child = top_child->child(W_child_index);
	  
	int W_child_pdgId =W_child->pdgId();
	//If the child is a W boson it is not a good one
	if( fabs(W_child_pdgId)==24 ){
	  top_child=W_child;
	  goto LOOP_W_BOSON;
	}

	bool is_leptonic_W= bool( fabs(W_child_pdgId)>=11 && fabs(W_child_pdgId)<=16 );  
	if( ! is_leptonic_W )continue;
	  
	//Fill here all relevant information
	//Fill also the W and top when you find a neutrino, fill the charged lepton only otherwise
	if( int(fabs(W_child_pdgId))%2 == 0 ){
	  pdgid_Wboson    = top_child->auxdata<int>("pdgId");
	  pdgid_aneutrino = W_child->auxdata<int>("pdgId");	    
	  tlv_atop        = GetTLorentzVector(top_particle);
	  tlv_Wboson      = GetTLorentzVector(top_child);
	  tlv_aneutrino   = GetTLorentzVector(W_child);
	}
	else if( int(fabs(W_child_pdgId))%2 == 1 ){
	  pdgid_lepton = W_child->auxdata<int>("pdgId");	    
	  tlv_lepton   = GetTLorentzVector(W_child);
	}
      }
    }
  }
  
  particle_list.push_back( std::make_pair(tlv_atop     ,pdgid_atop     ) );
  particle_list.push_back( std::make_pair(tlv_Wboson   ,pdgid_Wboson   ) );
  particle_list.push_back( std::make_pair(tlv_aneutrino,pdgid_aneutrino) );
  particle_list.push_back( std::make_pair(tlv_lepton   ,pdgid_lepton   ) );

  if(particle_list.size()!=4){
    std::cout<<" Error, the top decay chain does not have the expected number of entries, exiting "<<std::endl;
    exit(0);
  }

  return particle_list;
  
}

float TruthAna_LowMassHplusSelector::CalculateAngle(std::vector< std::pair< TLorentzVector,int > > particle_list, ANGLE angle)
{
  float variable=0.;
  
  TLorentzVector tlv_top;
  TLorentzVector tlv_Wboson;
  TLorentzVector tlv_lepton;
  TLorentzVector tlv_neutrino;

  for(auto particle : particle_list){
    
    if( int(fabs(particle.second))==6 )      tlv_top    = particle.first;
    else if( int(fabs(particle.second))==24) tlv_Wboson = particle.first;
    else if( fabs(particle.second)<24 && int(fabs(particle.second))%2==1 ) tlv_lepton  = particle.first;
    else if( fabs(particle.second)<24 && int(fabs(particle.second))%2==0 ) tlv_neutrino= particle.first;
    else{
      std::cout<<"Something wrong in the top decay chain tracking, exiting "<<std::endl;
      exit(0);
    }
  }
  
  bool is_TOPCM = bool(angle==ANGLE::LEP_TOP_DPHI_TOPCM || angle==ANGLE::NU_TOP_DPHI_TOPCM);
  if(is_TOPCM){
    tlv_Wboson  .Boost(-tlv_top.BoostVector());
    tlv_lepton  .Boost(-tlv_top.BoostVector());
    tlv_neutrino.Boost(-tlv_top.BoostVector());
    //Boost the top last 
    tlv_top     .Boost(-tlv_top.BoostVector());
  }

  if(angle==ANGLE::LEP_TOP_DPHI){
    variable=tlv_top.DeltaPhi(tlv_lepton);
  }
  else if(angle==ANGLE::LEP_TOP_DR){
    variable=tlv_top.DeltaR(tlv_lepton);
      }
  else if(angle==ANGLE::NU_TOP_DPHI){
    variable=tlv_top.DeltaPhi(tlv_neutrino);
      }
  else if(angle==ANGLE::NU_TOP_DR){
    variable=tlv_top.DeltaR(tlv_neutrino);
  }
  else if(angle==ANGLE::LEP_TOP_DPHI_TOPCM){
    variable=tlv_top.DeltaPhi(tlv_lepton);
  }
  else if(angle==ANGLE::NU_TOP_DPHI_TOPCM){
    variable=tlv_top.DeltaPhi(tlv_neutrino);
  }
  else{
    std::cout<<"ANGLE option not recognized exit now"<<std::endl;
    exit(0);
  }

  return variable;
}

TLorentzVector TruthAna_LowMassHplusSelector::GetTLorentzVector(const xAOD::TruthParticle* particle){
  
  TLorentzVector tlv_particle; 
  
  float p_px  = particle->auxdata<float>("px"); 
  float p_py  = particle->auxdata<float>("py"); 
  float p_pz  = particle->auxdata<float>("pz"); 
  float p_e   = particle->auxdata<float>("e"); 

  tlv_particle.SetPxPyPzE(p_px,p_py,p_pz,p_e);

  return tlv_particle; 
}

std::vector<TLorentzVector> TruthAna_LowMassHplusSelector::GetParticleInTopDecays(std::vector< std::pair< TLorentzVector,int > > particle_list, std::vector<int> pdgid_list){

  std::vector<TLorentzVector> v_tlv_particle; 

  for(auto particle : particle_list){
    for(unsigned int i=0; i<pdgid_list.size(); i++){
      if(particle.second==pdgid_list[i]){
	v_tlv_particle.push_back(particle.first);
      }
    }
  }
  
  if( v_tlv_particle.size() != 1 ){
    std::cout<<"No particles found from top decays"<<std::endl;
    exit(0);
  }
  
  return v_tlv_particle; 
}


float TruthAna_LowMassHplusSelector::CalculateCosThetaStar(const xAOD::TruthParticleContainer *cont) 
{
  float cos_theta_star;
  std::vector<TLorentzVector> tops;
  TLorentzVector ttbar;
  TLorentzVector top;
  TLorentzVector atop;

  for(auto vcont : *cont){
    int par_pdgid =vcont->auxdata<int>("pdgId");
    int par_status=vcont->auxdata<int>("status");
    if( fabs(par_pdgid) != 6  ) continue;
    if(par_status != 22) continue;
    float loc_px    =vcont->px(); float loc_py    =vcont->py();
    float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
    TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e);
    tops.push_back(tlv);
    if     (par_pdgid==  6 ) top =tlv;
    else if(par_pdgid== -6 ) atop=tlv;
    else {std::cout<<"In CalculateCosThetaStar > Expected a top or an anti-top, nothing is found"<<std::endl; exit(0);}
  }

  if(tops.size()!=2){ std::cout<<"In CalculateCosThetaStar > Expected to tops, found"<< tops.size() <<std::endl; exit(0); }
  
  ttbar=tops[0]+tops[1];

  top.Boost(-ttbar.BoostVector());
  atop.Boost(-ttbar.BoostVector());

  cos_theta_star=TMath::Cos(atop.Vect().Angle(ttbar.Vect()));

  return cos_theta_star;
}
