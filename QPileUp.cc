#include "QPileUp.h"
#include "Constants.h"
#include <string>

#include "TFile.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH3F.h"


using namespace Constants;

//=====================
QPileUp::QPileUp() {
  fDebug = 0;
  //fRho = NULL;
}
//=====================
//QPileUp::~QPileUp() {
//}
//=====================
void QPileUp::Make() {
  // InitMaps();
  if(fDebug>0) printf("QPileUp is initializing rho map... ");
  TH3F *fRho = new TH3F("rho","ChargeDensity [fC/cm^3];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
			kNAzimuthalSteps,0,TMath::TwoPi(),
			kNRadialSteps,kInnerRadius,kOutterRadius,
		      kNLongitudinalSteps,0,+kHalfLength);
  if(fDebug>0) printf("[DONE]\n");
  
  //------------
  //STEER
  if(fDebug>0) printf("QPileUp is being computed from TOY MODEL... \n");
  float e0 = 8.854187817e-12*1e+9; //[C]/[Vm]*1e+9
  float gas = 1.0/76628.0; //[Vs]
  float mult = 400.0; //950.0;
  float rate = 5e+4; //[Hz]
  double a=mult*rate*e0*gas; // fC/cm;
  float b=100.0/Constants::kHalfLength; //[1/m]
  float c=2.0/3.0*20.0;
  if(fDebug>1) printf("a = %f\n",a);
  if(fDebug>1) printf("b = %f\n",b);
  if(fDebug>1) printf("c = %f\n",c);
  for(int r=0; r!=Constants::kNRadialSteps; ++r) {
    float dr = fRho->GetXaxis()->GetBinCenter( r+1 )/100.0; //[m]
    for(int p=0; p!=Constants::kNAzimuthalSteps; ++p) {
      float dp = fRho->GetYaxis()->GetBinCenter( p+1 );
      for(int z=0; z!=Constants::kNLongitudinalSteps; ++z) {
	float dz = fRho->GetZaxis()->GetBinCenter( z+1 )/100.0; //[m]
	float dRho = a*(1-b*TMath::Abs(dz)+c)/(dr*dr); //fC/cm^3
	fRho->SetBinContent(p+1,r+1,z+1,dRho);
	if(fDebug>2) printf("@{Ir,Ip,Iz}={%d (%f),%d (%f),%d (%f)}, rho %f\n",r,dr,p,dp,z,dz,dRho);
      }
    }
  }
  if(fDebug>0) printf("[DONE]\n");
  //------------

    const char *outputfile= Form("%s_0.root",Constants::kFileNameRoot.data());
  if(fDebug>0) printf("QPileUp saving distortion maps... ");
  TFile *ofile = new TFile(outputfile,"RECREATE");
  ofile->WriteObject(fRho,"rho");
  ofile->Close();
  if(fDebug>0) printf("[DONE]\n");
  
  //SaveMaps();
}
//=====================
void QPileUp::InitMaps() {
  /*
  if(fDebug>0) printf("QPileUp is initializing rho map... ");
  fRho = new TH3F("rho","ChargeDensity [fC/cm^3];Radial [cm];Azimuthal [rad];Longitudinal [cm]",
		      kNRadialSteps,kInnerRadius,kOutterRadius,
		      kNAzimuthalSteps,0,TMath::TwoPi(),
		      kNLongitudinalSteps,0,+kHalfLength);
  if(fDebug>0) printf("[DONE]\n");
  */
}
//=====================
void QPileUp::SaveMaps() {
  /*
  const char *outputfile= Form("%s_0.root",Constants::kFileNameRoot.data());
  if(fDebug>0) printf("QPileUp saving distortion maps... ");
  TFile *ofile = new TFile(outputfile,"RECREATE");
  ofile->WriteObject(fRho,"rho");
  ofile->Close();
  if(fDebug>0) printf("[DONE]\n");
  */
}
