#include "MyTruthAnalysis/TruthAna_BaseSelector.h"
#include <TLorentzVector.h>
#include "MyTruthAnalysis/MyClasses.h"

TruthAna_BaseSelector::TruthAna_BaseSelector(){}


TruthAna_BaseSelector::~TruthAna_BaseSelector(){}


std::vector<TLorentzVector>  TruthAna_BaseSelector::DecayH(const xAOD::TruthParticleContainer *cont, float pt, float eta, DECAY_H decay_h)
{
  std::vector<TLorentzVector> tlv_vec;
  TLorentzVector decdec;
    //ALvaro
  TLorentzVector particle;
  TLorentzVector antiparticle;
    //
  decdec.SetPtEtaPhiM(10,1,2,1);
  int count = -1;
  //std::cout << " HERE " << std::endl;
  for(auto vcont : *cont){
    count++;
    float loc_px    =vcont->px(); float loc_py    =vcont->py();
    float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
    TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
    
    //int apdgid = fabs(vcont->auxdata<int>("pdgId"));
    //Alvaro
    int apdgid = vcont->auxdata<int>("pdgId");
      
    if( tlv.Pt()<pt )continue;
    if( fabs( tlv.Rapidity() )>eta )continue;

    if     ( decay_h==DECAY_H::TT  && (!(apdgid==top_pid))) continue;
    else if( decay_h==DECAY_H::BB  && (!(fabs(apdgid==bottom_pid)))) continue;
    else if( decay_h==DECAY_H::TB  && (!(apdgid==top_pid||apdgid==bottom_pid))) continue;
    //take all particles as parent but the Higgs bosons
    else if( decay_h==DECAY_H::ALL && (fabs(apdgid)==35||fabs(apdgid)==36||fabs(apdgid)==37||fabs(apdgid)==25) ) continue;

    //Perhaps that's not needed.. status will change from generator to generator
    //Check that you always have two particles form the H decays

    /*
    //intermediated particle, eg. tt decays from H
    if( decay_h==DECAY_H::TT && ( ! ( vcont->auxdata<int>("status")== 22 ) ) ) continue;
    //final state particle, eg. bb decays from H
    if( decay_h==DECAY_H::BB && ( ! ( vcont->auxdata<int>("status")== 23 ) ) ) continue;
    */
    
    //if (vcont->auxdata<int>("status") != 23) continue;
    //std::cout << "ADRI " << count << " " << tlv.Pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;
    bool has_a_higgs_parent=false;
    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
        //std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
        if(fabs(parent_pdgId)==35||fabs(parent_pdgId)==36||fabs(parent_pdgId)==37||fabs(parent_pdgId)==25) has_a_higgs_parent=true;
      }
    if(!has_a_higgs_parent) continue; 
    //Alvaro  
    if(apdgid>0) particle=tlv;
    if(apdgid<0) antiparticle=tlv;
    //  
    //tlv_vec.push_back(tlv);
  }
    //
    tlv_vec.push_back(particle);
    tlv_vec.push_back(antiparticle);
  if(!(tlv_vec.size()==2))
    {
      std::cout<<"Expected two H decay products, found n="<<tlv_vec.size()<<"), aborting"<<std::endl;
      exit(0);
    }
  else   
    {
      decdec=tlv_vec[0]+tlv_vec[1];
      //std::cout << "Found with DR:" << tlv_vec[0].DeltaR(tlv_vec[1]) << std::endl;
    }

  
  tlv_vec.push_back(decdec);
  
  return tlv_vec;
}


std::vector<TLorentzVector>  TruthAna_BaseSelector::AssociatedToH(const xAOD::TruthParticleContainer *cont, float pt, float eta, ASSOCIATED_H associated_h)
{
  std::vector<TLorentzVector> tlv_vec;
  TLorentzVector assoc;

  for(auto vcont : *cont){
    float loc_px    =vcont->px(); float loc_py    =vcont->py();
    float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
    TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 

    int apdgid = fabs(vcont->auxdata<int>("pdgId"));

    if( tlv.Pt()<pt )continue;
    if( fabs( tlv.Rapidity() )>eta )continue;

    if     ( associated_h==ASSOCIATED_H::TT && (!(apdgid==top_pid))) continue;
    else if( associated_h==ASSOCIATED_H::BB && (!(apdgid==bottom_pid))) continue;
    else if( associated_h==ASSOCIATED_H::TB && (!(apdgid==top_pid||apdgid==bottom_pid))) continue;
    else if( associated_h==ASSOCIATED_H::TC && (!(apdgid==top_pid||apdgid==charm_pid))) continue;
    //Debug printouts?
    //else if( associated_h==ASSOCIATED_H::INCLUSIVE ) continue; 

    if( ! ( vcont->auxdata<int>("status")== 22 ) ) continue;

    bool has_a_higgsORtop_parent=false;
    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
        if( parent_pdgId==35||parent_pdgId==36||parent_pdgId==37||parent_pdgId==25||parent_pdgId==top_pid)has_a_higgsORtop_parent=true;
      }
    if(has_a_higgsORtop_parent) continue;
    tlv_vec.push_back(tlv);                                                                   
  }
  if(!(tlv_vec.size()==2))
    {
      std::cout<<"Expected two objects in association with H, found ntops="<<tlv_vec.size()<<"), aborting"<<std::endl;
      exit(0);
    }
  assoc=tlv_vec[0]+tlv_vec[1];
  tlv_vec.push_back(assoc);
  
  return tlv_vec;
}

std::vector<myParticle>  TruthAna_BaseSelector::AssociatedToH_P(const xAOD::TruthParticleContainer *cont, float pt, float eta, ASSOCIATED_H associated_h)
{
  std::vector<myParticle> tlv_vec;
  myParticle assoc;
  myParticle aux;
  int count = -1;
  for(auto vcont : *cont){
    float loc_px    =vcont->px(); float loc_py    =vcont->py();
    float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
    TLorentzVector tlv; tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
    int apdgid = fabs(vcont->auxdata<int>("pdgId"));

    /*
    if( ! ( vcont->auxdata<int>("status")== 22 || vcont->auxdata<int>("status")== 21 || vcont->auxdata<int>("status")== 23 ) ) continue;
    count++;
    std::cout << "ADRI " << count << " " << tlv.Pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;

    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
        std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
      }
    */
    if( tlv.Pt()<pt )continue;
    if( fabs( tlv.Rapidity() )>eta )continue;

    if     ( associated_h==ASSOCIATED_H::TT && (!(apdgid==top_pid))) continue;
    else if( associated_h==ASSOCIATED_H::BB && (!(apdgid==bottom_pid))) continue;
    else if( associated_h==ASSOCIATED_H::TB && (!(apdgid==top_pid||apdgid==bottom_pid))) continue;
    else if( associated_h==ASSOCIATED_H::TC && (!(apdgid==top_pid||apdgid==charm_pid))) continue;
    //Debug printouts?
    //else if( associated_h==ASSOCIATED_H::INCLUSIVE ) continue; 

    if( ! ( vcont->auxdata<int>("status")== 22 ) ) continue;

    bool has_a_higgsORtop_parent=false;
    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
        if( parent_pdgId==35||parent_pdgId==36||parent_pdgId==37||parent_pdgId==25||parent_pdgId==top_pid)has_a_higgsORtop_parent=true;
      }
    if(has_a_higgsORtop_parent) continue;
    aux.FourVector=tlv;
    aux.pdgid = vcont->auxdata<int>("pdgId");
    tlv_vec.push_back(aux);                                                                   
  }
  if(!(tlv_vec.size()==2))
    {
      std::cout<<"Expected two objects in association with H, found ntops="<<tlv_vec.size()<<"), aborting"<<std::endl;
      exit(0);
    }
  assoc.FourVector=tlv_vec[0].FourVector+tlv_vec[1].FourVector;
  assoc.pdgid=-99;
  tlv_vec.push_back(assoc);
  
  return tlv_vec;
}

myParticle TruthAna_BaseSelector::SisterH(const xAOD::TruthParticleContainer *cont)
{
  myParticle tlp;
  TLorentzVector tlv; tlv.SetPtEtaPhiM(0,0,0,0);
  tlp.FourVector = tlv;
  tlp.pdgid = -99; 
  const xAOD::TruthParticle* parent;
  bool found_top = false;
  //int count = -1;
  for(auto vcont : *cont){
    //count++;  
    int pdgid = vcont->auxdata<int>("pdgId");

    /*
    int apdgid = fabs(vcont->auxdata<int>("pdgId"));
    if( ! ( vcont->auxdata<int>("status")== 22 || vcont->auxdata<int>("status")== 21 || vcont->auxdata<int>("status")== 23 ) ) continue;
    count++;
    std::cout << "ADRI " << count << " " << tlv.Pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;

    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
        std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
      }
    */

    if (pdgid!=25) continue;
    //std::cout << "ADRI " << count << " " << vcont->pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;

    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++){
      parent = vcont->parent(par_index);
      if(!parent) continue;
      int parent_pdgId=parent->pdgId();
      //std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
      if (fabs(parent_pdgId)==6) {
        found_top = true;
        break;
      }
    }
    if (found_top) break;
  }

  if (found_top){
    //std::cout << "Found parent " << parent->pdgId() << " " << parent->pt() << std::endl;
    for(unsigned int child_index=0; child_index < parent->nChildren(); child_index++){
      const xAOD::TruthParticle* child = parent->child(child_index);
      if(!child) continue;
      int child_pdgId=child->pdgId();
     // std::cout << "--> " << child_index << " " << child->pt() << " " << child_pdgId << " " << child->status() << std::endl;
      
      if ( (fabs(child_pdgId) == 1 || fabs(child_pdgId) == 2 || fabs(child_pdgId) == 3 || fabs(child_pdgId) == 4) && child->status()==23){
        float loc_px    =child->px(); float loc_py    =child->py();
        float loc_pz    =child->pz(); float loc_e     =child->e() ;
        tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
        tlp.FourVector = tlv;
        tlp.pdgid = child_pdgId;
        break;
      }
    }    
  }

  if (tlp.pdgid==-99){
    std::cout<<"Expected H sister, aborting"<<std::endl;
    exit(0);
  }
  return tlp;
}

myParticle TruthAna_BaseSelector::QuarkOtherTop(const xAOD::TruthParticleContainer *cont)
{
  myParticle tlp;
  TLorentzVector tlv; tlv.SetPtEtaPhiM(0,0,0,0);
  tlp.FourVector = tlv;
  tlp.pdgid = -99; 
  //int count = -1;
  for(auto vcont : *cont){
    //count++;  

    int Higgsinblood = -1;
    int pdgid = vcont->auxdata<int>("pdgId");

    /*
    int apdgid = fabs(vcont->auxdata<int>("pdgId"));
    if( ! ( vcont->auxdata<int>("status")== 22 || vcont->auxdata<int>("status")== 21 || vcont->auxdata<int>("status")== 23 ) ) continue;
    count++;
    std::cout << "ADRI " << count << " " << tlv.Pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;

    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
        std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
      }
    */

    if ( !(( fabs(pdgid) == 1 || fabs(pdgid) == 2 || fabs(pdgid) == 3 || fabs(pdgid) == 4 || fabs(pdgid) == 5 ) && vcont->auxdata<int>("status") == 23) )  continue;
    //std::cout << "ADRI " << count << " " << vcont->pt() << " " << vcont->auxdata<int>("pdgId") << " " << vcont->auxdata<int>("status") << std::endl;

    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++){
      const xAOD::TruthParticle* parent = vcont->parent(par_index);
      if(!parent) continue;
      int parent_pdgId=parent->pdgId();
      //std::cout << "--> " << par_index << " " << parent->pt() << " " << parent_pdgId << " " << parent->status() << std::endl;
      if (fabs(parent_pdgId)==6) {
        Higgsinblood = 0;
        for(unsigned int child_index=0; child_index < parent->nChildren(); child_index++){
          const xAOD::TruthParticle* child = parent->child(child_index);
          if(!child) continue;
          int child_pdgId=child->pdgId();
          //std::cout << "---> " << child_index << " " << child->pt() << " " << child_pdgId << " " << child->status() << std::endl;
          if (child_pdgId==25){
            Higgsinblood = 1;
            break;
          }
        }
      }
     if(Higgsinblood == 1) break;
     else if (Higgsinblood == 0){
        float loc_px    =vcont->px(); float loc_py    =vcont->py();
        float loc_pz    =vcont->pz(); float loc_e     =vcont->e() ;
        tlv.SetPxPyPzE(loc_px,loc_py,loc_pz,loc_e); 
        tlp.FourVector = tlv;
        tlp.pdgid = pdgid;
        break;
     }
    } 
  }

  if (tlp.pdgid==-99){
    std::cout<<"Expected q from othertop, aborting"<<std::endl;
    exit(0);
  }
  return tlp;
}

std::vector<int> TruthAna_BaseSelector::GetHiggsDecays(const xAOD::TruthParticleContainer *cont, int H_pdgid)
{
  std::vector<int> PDGID;

  for(auto vcont : *cont){

    int par_pdgid=vcont->auxdata<int>("pdgId");
    
    if(par_pdgid == H_pdgid) continue;

    bool has_a_higgs_parent=false;
    for(unsigned int par_index=0; par_index < vcont->nParents(); par_index++)
      {
        const xAOD::TruthParticle* parent = vcont->parent(par_index);
        if(!parent) continue;
        int parent_pdgId=parent->pdgId();
	//fabs used for H+
        if( fabs(parent_pdgId) == H_pdgid ) has_a_higgs_parent=true;
									      
      }
    if(!has_a_higgs_parent) continue;
    
    //std::cout<<" Higgs doughter "<<par_pdgid<<std::endl;
    PDGID.push_back(par_pdgid);                                                                   
  }
  
  return PDGID;
}

