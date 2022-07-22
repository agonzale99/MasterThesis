#ifndef MYTRUTHANALYSIS_TRUTHANA_TTHTTSELECTOR_H
#define MYTRUTHANALYSIS_TRUTHANA_TTHTTSELECTOR_H

#include "MyTruthAnalysis/TruthAna_BaseSelector.h"

class TruthAna_ttHttSelector : public TruthAna_BaseSelector  
{
 public: 
 TruthAna_ttHttSelector():TruthAna_BaseSelector(DECAY_H::TT,ASSOCIATED_H::TT){};
  ~TruthAna_ttHttSelector(){};
};

#endif 
