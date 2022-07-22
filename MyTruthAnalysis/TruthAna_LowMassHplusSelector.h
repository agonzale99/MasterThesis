#ifndef MYTRUTHANALYSIS_TRUTHANA_LOWMASSHPLUSSELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_LOWMASSHPLUSSELECTOR_H

#include "MyTruthAnalysis/TruthAna_BaseSelector.h"

class TruthAna_LowMassHplusSelector : public TruthAna_BaseSelector  
{
 public: 
 TruthAna_LowMassHplusSelector():TruthAna_BaseSelector(DECAY_H::ALL,ASSOCIATED_H::INCLUSIVE){};
  ~TruthAna_LowMassHplusSelector(){};

  enum class ANGLE {LEP_TOP_DPHI,LEP_TOP_DR,NU_TOP_DPHI,NU_TOP_DR,
      LEP_TOP_DPHI_TOPCM,NU_TOP_DPHI_TOPCM};

  //Study the top spin correlations
  //The first element of the pair is the TLV of the parton, the second element is the pdgid
  std::vector< std::pair< TLorentzVector,int > > GetTopDecayChain(const xAOD::TruthParticleContainer *cont);
  float CalculateAngle(std::vector< std::pair< TLorentzVector,int > > particle_list, ANGLE angle);
  TLorentzVector GetTLorentzVector(const xAOD::TruthParticle* particle);
  std::vector<TLorentzVector> GetParticleInTopDecays(std::vector< std::pair< TLorentzVector,int > > particle_list, std::vector<int> pdgid_list);
  float CalculateCosThetaStar(const xAOD::TruthParticleContainer *cont);

};

#endif 
