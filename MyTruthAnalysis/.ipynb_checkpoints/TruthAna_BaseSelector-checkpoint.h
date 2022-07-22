#ifndef MYTRUTHANALYSIS_TRUTHANA_BASESELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_BASESELECTOR_H

#include <vector>
#include <TLorentzVector.h>

#include "xAODTruth/TruthParticleContainer.h"
#include "MyTruthAnalysis/MyClasses.h"

class TruthAna_BaseSelector
{
 public: 
  TruthAna_BaseSelector();
  ~TruthAna_BaseSelector();
  //super cool c++11 feature
  enum class DECAY_H {ALL,TT,BB,TB,TC};
  enum class ASSOCIATED_H {INCLUSIVE,TT,BB,TB,TC};
  //Selector methods
  std::vector<TLorentzVector> DecayH(const xAOD::TruthParticleContainer *cont, float pt, float eta, DECAY_H decay_h);
  std::vector<TLorentzVector> AssociatedToH(const xAOD::TruthParticleContainer *cont, float pt, float eta, ASSOCIATED_H associated_h);
  std::vector<myParticle> AssociatedToH_P(const xAOD::TruthParticleContainer *cont, float pt, float eta, ASSOCIATED_H associated_h);
  myParticle SisterH(const xAOD::TruthParticleContainer *cont); //finds higgs parent and returns daughter that is not the higgs 
  myParticle QuarkOtherTop(const xAOD::TruthParticleContainer *cont); //finds quark not related to the higgs from top 
  std::vector<int> GetHiggsDecays(const xAOD::TruthParticleContainer *cont, int H_pdgid);

 private:
  const int top_pid=6;
  const int bottom_pid=5;
  const int charm_pid=4;
  DECAY_H ini_decay_h;
  ASSOCIATED_H ini_associated_h;
  
 protected:
 TruthAna_BaseSelector(DECAY_H ini_decay_h,ASSOCIATED_H ini_associated_h) : 
  ini_decay_h(ini_decay_h),ini_associated_h(ini_associated_h){};
  
};



#endif 
