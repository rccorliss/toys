#ifndef __QPILEUP_H__
#define __QPILEUP_H__
//#include "TH3F.h"
//#include "TVector3.h"
#include "TH3F.h"
#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"

class QPileUp {
 public:
  QPileUp();
  void SetDebugLevel(int n) {fDebug=n;};
  void Make();

 protected:
  void InitMaps();
  void SaveMaps();
  int fDebug;

  //TH3F *fRho;
};

#endif /* __QPILEUP_H__ */
