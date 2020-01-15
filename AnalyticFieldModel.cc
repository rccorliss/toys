 
#include "TVector3.h"
#include "TFormula.h"
#include "AnalyticFieldModel.h"

AnalyticFieldModel::AnalyticFieldModel(float _ifc_radius, float _ofc_radius, float _z_max, float scalefactor){
  double ifc_radius=_ifc_radius;
  double ofc_radius=_ofc_radius;
  double tpc_halfz=_z_max;
  

  double sum=ifc_radius+ofc_radius;//338 in ALICE, [3] in args
  double prod=ifc_radius*ofc_radius;//21250.75 in ALICE [4] in args
  double diff=ofc_radius-ofc_radius;
  
  double a = ofc_radius * ofc_radius;
  a *= (diff);
  a *= (diff);
  a =  (1000.0 / a);
  double b = 0.5;
  double c = 1.0 / (((tpc_halfz) / 2.0) * ((tpc_halfz) / 2.0));
  double d=sum;
  double e=prod;
  vTestFunction1.SetParameters(a, b, c, d, e);
  rhoTestFunction1.SetParameters(a, b, c, d, e);

  erTestFunction1.SetParameters(-a, b, c, d, e);
  ePhiTestFunction1.SetParameters(-a, b, c, d, e);
  ezTestFunction1.SetParameters(-a, b, c, d, e);
  intErDzTestFunction1.SetParameters(-a, b, c, d, e);
  intEPhiDzTestFunction1.SetParameters(-a, b, c, d, e);
  intEzDzTestFunction1.SetParameters(-a, b, c, d, e);
  return;
}

TVector3 AnalyticFieldModel::E(TVector3 pos){//field as a function of position
  //in rhat phihat zhat coordinates: (at phi=0, phi is the +Y position, Perp is the +X direction and Z is Z)
  TVector3 ret(erTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z()),
	       ePhiTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z()),
	       ezTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z()));
  //now rotate this to the position we evaluated it at, to match the global coordinate system.
  ret.RotateZ(pos.Phi());
  return ret;
}
	       
  
double AnalyticFieldModel::rho(TVector3 pos){//charge density as a function of position
  //at phi=0, phi is the +Y position, Perp is the +X direction and Z is Z.
  return rhoTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z());
}
TVector3 AnalyticFieldModel::Eint(float zfinal, TVector3 pos){//field integral from 'pos' to z-position zfinal.
  //in rhat phihat zhat coordinates: (at phi=0, phi is the +Y position, Perp is the +X direction and Z is Z)
  TVector3 eintI(intErDzTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z()),
	       intEPhiDzTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z()),
	       intEzDzTestFunction1.Eval(pos.Perp(),pos.Phi(),pos.Z()));
  TVector3 eintF(intErDzTestFunction1.Eval(pos.Perp(),pos.Phi(),zfinal),
	       intEPhiDzTestFunction1.Eval(pos.Perp(),pos.Phi(),zfinal),
	       intEzDzTestFunction1.Eval(pos.Perp(),pos.Phi(),zfinal));

  TVector3 ret=eintF-eintI;
    //now rotate this to the position we evaluated it at, to match the global coordinate system.
  ret.RotateZ(pos.Phi());
  return ret;
}
