#ifndef MYTRUTHANALYSIS_TRUTHANA_CUTFLOWTOOLS_H
#define MYTRUTHANALYSIS_TRUTHANA_CUTFLOWTOOLS_H

#include <string>
#include <map>
#include <EventLoop/Algorithm.h>
#include <TH1F.h>
#include <iostream>
#include <sstream>  

class TruthAna_CutFlowTools
{
 
 public:
  TruthAna_CutFlowTools();
  ~TruthAna_CutFlowTools();

  std::map<std::string, std::vector<std::string> > m_vecCutFlowID; //!                      
  std::map<std::string, TH1F*> m_vecCutFlows; //!                                                
  EL::StatusCode addCutFlow(std::string cutFlowID, std::string sCutIDs);
  EL::StatusCode callAddCutFlow(std::string cutFlowID, std::string sCutIDs);
  EL::StatusCode cutFlow(std::string cutFlowID, std::string cutID, double weight);
  EL::StatusCode printCutFlow(std::string cutFlowID);
  void printAllCutFlow();
  std::string slength(std::string input, int size);
  std::string slength(double d, int size);
  
};

#endif 
