#include "MyTruthAnalysis/TruthAna_ZbSelector.h"


TruthAna_ZbSelector::TruthAna_ZbSelector(){}


TruthAna_ZbSelector::~TruthAna_ZbSelector(){}


std::vector<TLorentzVector> TruthAna_ZbSelector::GetZboson(std::vector<TLorentzVector> leptons)
{
  TLorentzVector Zboson=leptons[0]+leptons[1];
  //dirty trick .. need to fix it
  std::vector<TLorentzVector> v_Zboson; v_Zboson.push_back(Zboson);

  return v_Zboson;
}


std::vector<TLorentzVector> TruthAna_ZbSelector::GetBJets(std::vector<TLorentzVector> bpartons, std::vector<TLorentzVector> jets)
{
  std::vector<TLorentzVector> tlv_vec_bjets;
  
  unsigned int jets_n = jets.size();
  unsigned int partons_n = bpartons.size();

  std::cout<<"I'm here"<<std::endl;
  
  for(unsigned int ijet=0; ijet<jets_n; ijet++)
    {

      //progress_jet:
      
      for(unsigned int ipartons=0; ipartons<partons_n; ipartons++)
	{
	  std::cout<<"I'm here 2"<<std::endl;
	  
	  //float DeltaR= jets[ijet].DeltaR(bpartons[ipartons]);
	  float DeltaR= 0.1;
	  if(DeltaR<0.3)
	    {

	      std::cout<<"I'm here 3"<<std::endl;
	      
	      tlv_vec_bjets.push_back(jets[ijet]);
	      //goto progress_jet;
	    }
	}
    }
  
  return tlv_vec_bjets;
}


