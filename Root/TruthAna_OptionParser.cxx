#include <MyTruthAnalysis/TruthAna_OptionParser.h>

#include <string>

TruthAna_OptionParser::TruthAna_OptionParser():
  m_Process("ttH"),
  m_InputList(""),
  m_Derivation("TOPQ1")
{}


TruthAna_OptionParser::~TruthAna_OptionParser(){}


void TruthAna_OptionParser::ReadOptions()
{
  return;
}
