#ifndef MYTRUTHANALYSIS_TRUTHANA_CLASSES_H
#define MYTRUTHANALYSIS_TRUTHANA_CLASSES_H
class myParticle{
  public:
    TLorentzVector FourVector;
    Int_t pdgid;
  
  myParticle(){ pdgid=-99;};
  ~myParticle(){};
};

struct ComparePt{
  bool operator()(TLorentzVector a, TLorentzVector b) {
    return (a.Pt() > b.Pt());
  }
};

#endif