#ifndef MYTRUTHANALYSIS_TRUTHANA_ZBSELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_ZBSELECTOR_H

#include <vector>
#include <TLorentzVector.h>

#include "xAODTruth/TruthParticleContainer.h"

class TruthAna_ZbSelector
{
 public: 
  TruthAna_ZbSelector();
  ~TruthAna_ZbSelector();

  //Selector methods

  std::vector<TLorentzVector> GetZboson(std::vector<TLorentzVector> leptons);
  std::vector<TLorentzVector> GetBJets(std::vector<TLorentzVector> bpartons, std::vector<TLorentzVector> jets);
    
};


#endif 
