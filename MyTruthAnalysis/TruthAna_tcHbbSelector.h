#ifndef MYTRUTHANALYSIS_TRUTHANA_TCHBBSELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_TCHBBSELECTOR_H

#include "MyTruthAnalysis/TruthAna_BaseSelector.h"

class TruthAna_tcHbbSelector : public TruthAna_BaseSelector  
{
 public: 
 TruthAna_tcHbbSelector():TruthAna_BaseSelector(DECAY_H::BB,ASSOCIATED_H::TC){};
  ~TruthAna_tcHbbSelector(){};
  
  inline std::vector<TLorentzVector> 
    GetParticleW(std::vector<TLorentzVector> v_tlv_el,std::vector<TLorentzVector> v_tlv_nu)
  {
    std::vector<TLorentzVector> v_tlv_wboson;
    v_tlv_wboson.push_back(v_tlv_el[0]+v_tlv_nu[0]);
    return v_tlv_wboson;
  }
  
};

#endif 
