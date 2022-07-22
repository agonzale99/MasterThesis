#ifndef MYTRUTHANALYSIS_TRUTHANA_BBHTTSELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_BBHTTSELECTOR_H

#include "MyTruthAnalysis/TruthAna_BaseSelector.h"

class TruthAna_bbHttSelector : public TruthAna_BaseSelector  
{
 public: 
 TruthAna_bbHttSelector():TruthAna_BaseSelector(DECAY_H::TT,ASSOCIATED_H::BB){};
  ~TruthAna_bbHttSelector(){};
};

#endif 
