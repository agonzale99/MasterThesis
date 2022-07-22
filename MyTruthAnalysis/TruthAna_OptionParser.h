#ifndef MYTRUTHANALYSIS_TRUTHANA_OPTIONPARSER_H
#define MYTRUTHANALYSIS_TRUTHANA_OPTIONPARSER_H

#include <string>

class TruthAna_OptionParser
{
 
 public:
  TruthAna_OptionParser();
  ~TruthAna_OptionParser();

  inline std::string Process() const { return m_Process; }
  inline std::string InputList() const { return m_InputList; }  
  inline std::string Derivation() const { return m_Derivation; }  

  enum PROCESS{BBH_BBTT,TTH_TTTT,TT_TCH};
  enum DERIVATION{TRUTH0,TRUTH1};
  
  PROCESS proc;
  DERIVATION deriv;
  
  void ReadOptions();
 
 private:
  std::string m_Process;  
  std::string m_InputList;
  std::string m_Derivation;

};

#endif 
