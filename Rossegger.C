#include "Rossegger.h"

#include <iostream>
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include <string>

#include "TH2.h"
#include "TH3.h"

using namespace std;
using namespace TMath;


/*
  This is a renamed copy of Carlos and Tom's Spacecharge class, with phi terms added.
 */

//  Here we create the interface to the fortran code we got off the web...
extern"C"{
  void dkia_(int *IFAC, double *X, double *A, double *DKI, double *DKID, int *IERRO);
  void dlia_(int *IFAC, double *X, double *A, double *DLI, double *DLID, int *IERRO);
}


Rossegger::Rossegger(std::string filename) {
  TFile *file = new TFile(filename.data());
  if (!file)
    {
      cout << "Could Not Open Greene Function File: " << filename.data() << endl;
      return;
    }
  fByFile = true;

  Tags = (TH2*)file->Get("Tags");
  if (!Tags)
    {
      cout << "Could not find Tags histogram in the Greenes function file." << endl;
      return;
    }

  int Nx = Tags->GetXaxis()->GetNbins();
  int Ny = Tags->GetYaxis()->GetNbins();
  char name[500];
  for (int i=1; i<Nx+1; i++)
    {
      for (int j=1; j<Ny+1; j++)
	{
	  int index = Tags->GetBinContent(i,j);
	  sprintf(name,"Gr%08d",index);
	  if (verbosity) cout << " " << name;
	  Grid[name] = (TH3*)file->Get(name);
	  if (verbosity) cout << " " << Grid[name]->GetName();
	  sprintf(name,"Gz%08d",index);
	  if (verbosity) cout << " " << name;
	  Grid[name] = (TH3*)file->Get(name);
	  if (verbosity) cout << " " << Grid[name]->GetName();
	  if (verbosity) cout << endl;
	}
    }
  
  MinimumDR = 2.0*Tags->GetXaxis()->GetBinWidth(1);
  MinimumDZ = 2.0*Tags->GetYaxis()->GetBinWidth(1);

  a = Tags->GetXaxis()->GetBinCenter(1);
  b = Tags->GetXaxis()->GetBinCenter(Nx);
  L = Tags->GetYaxis()->GetBinCenter(Ny);

  verbosity = 0;
  tweak = true;
  pi = 2.0 * asin(1.0);
  cout << pi << endl;

  FindBetamn(0.0001);
  FindMunk(0.0001);

  cout << "Rossegger object initialized as follows:" << endl;
  cout << "  Inner Radius = " << a << " cm." << endl;
  cout << "  Outer Radius = " << b << " cm." << endl;
  cout << "  Half  Length = " << L << " cm." << endl;


}

Rossegger::Rossegger(double InnerRadius, double OuterRadius, double HalfLength)
{
  a = InnerRadius;
  b = OuterRadius;
  L = HalfLength;

  verbosity = 0;
  tweak = true;
  pi = 2.0 * asin(1.0);
  cout << pi << endl;

  FindBetamn(0.0001);
  FindMunk(0.0001);

  cout << "Rossegger object initialized as follows:" << endl;
  cout << "  Inner Radius = " << a << " cm." << endl;
  cout << "  Outer Radius = " << b << " cm." << endl;
  cout << "  Half  Length = " << L << " cm." << endl;

  fByFile = false;
  return ;
}


void Rossegger::FindBetamn(double epsilon)
{
  cout << "Now filling the Beta[m][n] Array..."<<endl;
  double l = a/b;
  for (int m=0; m<NumberOfOrders; m++)
    {
      if (verbosity) cout << endl << m;
      int n=0;  //  !!!  Off by one from Rossegger convention  !!!
      double x=epsilon;
      double previous = jn(m,x)*yn(m,l*x) - jn(m,l*x)*yn(m,x);
      while (n < NumberOfOrders)
	{
	  //  Rossegger equation 5.12
	  double value = jn(m,x)*yn(m,l*x) - jn(m,l*x)*yn(m,x);
	  //if (verbosity) cout << " " << value;
	  if (value == 0) cout << "FFFFFFUUUUUUUUUUU......" << endl;
	  if ( (value == 0) || (value<0 && previous>0) || (value>0 && previous<0))
	    {
	      double slp = (value-previous)/epsilon;
	      double cpt =  value - slp*x;
	      double x0  = -cpt/slp;

	      Betamn[m][n] = x0/b;
	      if (verbosity) cout << " " << x0 << "," << Betamn[m][n];
	      n++;
	    }
	  previous = value;
	  x+=epsilon;
	}
    }


  //  Now fill in the N2mn array...
  for (int m=0; m<NumberOfOrders; m++)
    {
      for (int n=0; n<NumberOfOrders; n++)
	{
	  //  Rossegger Equation 5.17
	  //  N^2_mn = 2/(pi * beta)^2 [ {Jm(beta a)/Jm(beta b)}^2 - 1 ]
	  N2mn[m][n]  = 2/(pi*pi*Betamn[m][n]*Betamn[m][n]);
	  N2mn[m][n] *= (jn(m,Betamn[m][n]*a)*jn(m,Betamn[m][n]*a)/jn(m,Betamn[m][n]*b)/jn(m,Betamn[m][n]*b) - 1.0);
	  if (verbosity) cout << "m: " << m << " n: " << n << " N2[m][n]: " << N2mn[m][n];
	  double integral=0.0;
	  double step = 0.01;
	  for (double r=a; r<b; r+=step)
	    {
	      integral += Rmn(m,n,r)*Rmn(m,n,r)*r*step;
	    }
	  if (verbosity) cout << " Int: " << integral << endl;
	  //N2mn[m][n] = integral;
	}
    }

  cout << "Done." << endl;
}


void Rossegger::FindMunk(double epsilon)
{
  cout << "Now filling the Mu[n][k] Array..."<<endl;
  cout << "Done." << endl;
}


double Rossegger::Rmn(int m, int n, double r)
{
  if (verbosity) cout << "Determine Rmn("<<m<<","<<n<<","<<r<<") = ";

  //  Check input arguments for sanity...
  int error=0;
  if (m<0 || m>NumberOfOrders) error=1;
  if (n<0 || n>NumberOfOrders) error=1;
  if (r<a || r>b)              error=1;
  if (error)
    {
      cout << "Invalid arguments Rmn("<<m<<","<<n<<","<<r<<")" << endl;;
      return 0;
    }

  //  Calculate the function using C-libraries from math.h
  //  Rossegger Equation 5.11:
  //         Rmn(r) = Ym(Beta_mn a)*Jm(Beta_mn r) - Jm(Beta_mn a)*Ym(Beta_mn r)
  double R=0;
  R = yn(m,Betamn[m][n]*a)*jn(m,Betamn[m][n]*r) - jn(m,Betamn[m][n]*a)*yn(m,Betamn[m][n]*r);

  if (verbosity) cout << R << endl;
  return R;
}

double Rossegger::Rmn1(int m, int n, double r)
{
 //  Check input arguments for sanity...
  int error=0;
  if (m<0 || m>NumberOfOrders) error=1;
  if (n<0 || n>NumberOfOrders) error=1;
  if (r<a || r>b)              error=1;
  if (error)
    {
      cout << "Invalid arguments Rmn1("<<m<<","<<n<<","<<r<<")" << endl;;
      return 0;
    }

  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.32
  //         Rmn1(r) = Km(BetaN a)Im(BetaN r) - Im(BetaN a) Km(BetaN r)
  double R=0;
  double BetaN = (n+1)*pi/L;
  R = BesselK(m,BetaN*a)*BesselI(m,BetaN*r)-BesselI(m,BetaN*a)*BesselK(m,BetaN*r);

  return R;
}

double Rossegger::Rmn2(int m, int n, double r)
{
 //  Check input arguments for sanity...
  int error=0;
  if (m<0 || m>NumberOfOrders) error=1;
  if (n<0 || n>NumberOfOrders) error=1;
  if (r<a || r>b)              error=1;
  if (error)
    {
      cout << "Invalid arguments Rmn2("<<m<<","<<n<<","<<r<<")" << endl;;
      return 0;
    }

  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.33
  //         Rmn2(r) = Km(BetaN b)Im(BetaN r) - Im(BetaN b) Km(BetaN r)
  double R=0;
  double BetaN = (n+1)*pi/L;
  R = BesselK(m,BetaN*b)*BesselI(m,BetaN*r)-BesselI(m,BetaN*b)*BesselK(m,BetaN*r);

  return R;
}

double Rossegger::RPrime(int m, int n, double ref, double r)
{
 //  Check input arguments for sanity...
  int error=0;
  if (m<0   || m>NumberOfOrders) error=1;
  if (n<0   || n>NumberOfOrders) error=1;
  if (ref<a || ref>b)            error=1;
  if (r<a   || r>b)              error=1;
  if (error)
    {
      cout << "Invalid arguments RPrime("<<m<<","<<n<<","<<ref<<","<<r<<")" << endl;;
      return 0;
    }

  double R=0;
  //  Calculate using the TMath functions from root.
  //  Rossegger Equation 5.65
  //         Rmn2(ref,r) = BetaN/2* [   Km(BetaN ref) {Im-1(BetaN r) + Im+1(BetaN r)}
  //                                  - Im(BetaN ref) {Km-1(BetaN r) + Km+1(BetaN r)}  ]
  //  NOTE:  K-m(z) = Km(z) and I-m(z) = Im(z)... 
  //
  //
  double BetaN = (n+1)*pi/L;
  double term1 = BesselK(m,BetaN*ref)*( BesselI(abs(m-1),BetaN*r) + BesselI(m+1,BetaN*r) );
  double term2 = BesselI(m,BetaN*ref)*( BesselK(abs(m-1),BetaN*r) + BesselK(m+1,BetaN*r) );
  R = BetaN/2.0*(term1 + term2);

  return R;
}

double Rossegger::Rnk(int n, int k, double r)
{
 //  Check input arguments for sanity...
  int error=0;
  if (n<0   || n>NumberOfOrders) error=1;
  if (k<0   || k>NumberOfOrders) error=1;
  if (r<a   || r>b)              error=1;
  if (error)
    {
      cout << "Invalid arguments Rnk("<<n<<","<<k<<","<<r<<")" << endl;;
      return 0;
    }

  //  Rossegger Equation 5.45
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN r) - Kimu_nk(BetaN a) Limu_nk (BetaN r)
  double BetaN = (n+1)*pi/L;

  int IFAC=1;
  double A=Munk[n][k];
  double DLI_1=0; 
  double DKI_1=0; 
  double DLI_2=0; 
  double DKI_2=0; 
  double DERR=0;
  int IERRO=0;

  double X=BetaN*a;
  dlia_( &IFAC, &X, &A, &DLI_1, &DERR, &IERRO);

  X=BetaN*r;
  dkia_( &IFAC, &X, &A, &DKI_1, &DERR, &IERRO);

  X=BetaN*a;
  dkia_( &IFAC, &X, &A, &DKI_2, &DERR, &IERRO);

  X=BetaN*r;
  dlia_( &IFAC, &X, &A, &DLI_2, &DERR, &IERRO);

  double R= DLI_1*DKI_1 - DKI_2*DLI_2;

  return R;

}

double Rossegger::ByFileEZ(double r, double phi, double z, double r1, double phi1, double z1)
{
  int index = Tags->GetBinContent(Tags->FindBin(r1,z1));
  char name[500];

  sprintf(name,"Gz%08d",index);

  //  "tweak" means that we shall move the r point so that r-r' better matches the grid...
  if (tweak)
    {
      int xbin,ybin,zbin;
      Tags->GetBinXYZ(Tags->FindBin(r1,z1),xbin,ybin,zbin);
      double xc = Tags->GetXaxis()->GetBinCenter(xbin);
      double yc = Tags->GetYaxis()->GetBinCenter(ybin);
      r = xc + (r-r1);
      z = yc + (z-z1);

      double epsilon = 0.01;
      if (r<a+epsilon) r=a+epsilon;
      if (r>b-epsilon) r=b-epsilon;
      if (z<0+epsilon) z=0+epsilon;
      if (z>L-epsilon) z=L-epsilon;
    }

  double dphi = phi-phi1;
  while (dphi<-pi) {dphi += 2.0*pi;}
  while (dphi> pi) {dphi -= 2.0*pi;}

  double answer = Grid[name]->Interpolate(r,dphi,z);
  return answer;
}

double Rossegger::ByFileER(double r, double phi, double z, double r1, double phi1, double z1)
{
  int index = Tags->GetBinContent(Tags->FindBin(r1,z1));
  char name[500];

  sprintf(name,"Gr%08d",index);

  if (tweak)
    {
      int xbin,ybin,zbin;
      Tags->GetBinXYZ(Tags->FindBin(r1,z1),xbin,ybin,zbin);
      double xc = Tags->GetXaxis()->GetBinCenter(xbin);
      double yc = Tags->GetYaxis()->GetBinCenter(ybin);
      r = xc + (r-r1);
      z = yc + (z-z1);

      double epsilon = 0.01;
      if (r<a+epsilon) r=a+epsilon;
      if (r>b-epsilon) r=b-epsilon;
      if (z<0+epsilon) z=0+epsilon;
      if (z>L-epsilon) z=L-epsilon;
    }

  double dphi = phi-phi1;
  while (dphi<-pi) {dphi += 2.0*pi;}
  while (dphi> pi) {dphi -= 2.0*pi;}

  double answer = Grid[name]->Interpolate(r,dphi,z);
  return answer;
}

double Rossegger::Ez(double r, double phi, double z, double r1, double phi1, double z1)
{
  if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileEZ(r,phi,z,r1,phi1,z1);
  //  Check input arguments for sanity...
  int error=0;
  if (r<a    || r>b)         error=1;
  if (phi<0  || phi>2*pi)    error=1;
  if (z<0    || z>L)         error=1;
  if (r1<a   || r1>b)        error=1;
  if (phi1<0 || phi1>2*pi)   error=1;
  if (z1<0   || z1>L)        error=1;
  if (error)
    {
      cout << "Invalid arguments Ez(";
      cout <<r<<",";
      cout <<phi<<",";
      cout <<z<<",";
      cout <<r1<<",";
      cout <<phi1<<",";
      cout <<z1;
      cout <<")" << endl;;
      return 0;
    }

  double G=0;
  for (int m=0; m<NumberOfOrders; m++)
    {
      if (verbosity) cout << endl << m;
      for (int n=0; n<NumberOfOrders; n++)
	{
	  if (verbosity) cout << " " << n;
	  double term = 1/(2.0*pi);
	  if (verbosity) cout << " " << term; 
	  term *= (2 - ((m==0)?1:0))*cos(m*(phi-phi1));
	  if (verbosity) cout << " " << term; 
	  term *= Rmn(m,n,r)*Rmn(m,n,r1)/N2mn[m][n];
	  if (verbosity) cout << " " << term; 
	  if (z<z1)
	    {
	      term *=  cosh(Betamn[m][n]*z)*sinh(Betamn[m][n]*(L-z1))/sinh(Betamn[m][n]*L);
	    }
	  else
	    {
	      term *= -cosh(Betamn[m][n]*(L-z))*sinh(Betamn[m][n]*z1)/sinh(Betamn[m][n]*L);;
	    }
	  if (verbosity) cout << " " << term; 
	  G += term;
	  if (verbosity) cout << " " << term << " " << G << endl;
	}
    }
  if (verbosity) cout << "Ez = " << G << endl;

  return G;
}


double Rossegger::Er(double r, double phi, double z, double r1, double phi1, double z1)
{
  if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileER(r,phi,z,r1,phi1,z1);
  //  Check input arguments for sanity...
  int error=0;
  if (r<a    || r>b)         error=1;
  if (phi<0  || phi>2*pi)    error=1;
  if (z<0    || z>L)         error=1;
  if (r1<a   || r1>b)        error=1;
  if (phi1<0 || phi1>2*pi)   error=1;
  if (z1<0   || z1>L)        error=1;
  if (error)
    {
      cout << "Invalid arguments Er(";
      cout <<r<<",";
      cout <<phi<<",";
      cout <<z<<",";
      cout <<r1<<",";
      cout <<phi1<<",";
      cout <<z1;
      cout <<")" << endl;;
      return 0;
    }

  double G=0;
  for (int m=0; m<NumberOfOrders; m++)
    {
      for (int n=0; n<NumberOfOrders; n++)
	{
	  double term = 1/(L*pi);

	  term *= (2 - ((m==0)?1:0))*cos(m*(phi-phi1));

	  double BetaN = (n+1)*pi/L;
	  term *= sin(BetaN*z)*sin(BetaN*z1);

	  if (r<r1)
	    {
	      term *= RPrime(m,n,a,r)*Rmn2(m,n,r1);
	    }
	  else
	    {
	      term *= Rmn1(m,n,r1)*RPrime(m,n,b,r);
	    }

	  term /= BesselI(m,BetaN*a)*BesselK(m,BetaN*b)-BesselI(m,BetaN*b)*BesselK(m,BetaN*a);

	  G += term;
	}
    }

  if (verbosity) cout << "Er = " << G << endl;

  return G;
}

double Rossegger::Ephi(double r, double phi, double z, double r1, double phi1, double z1)
{
  //  Check input arguments for sanity...
  int error=0;
  if (r<a    || r>b)         error=1;
  if (phi<0  || phi>2*pi)    error=1;
  if (z<0    || z>L)         error=1;
  if (r1<a   || r1>b)        error=1;
  if (phi1<0 || phi1>2*pi)   error=1;
  if (z1<0   || z1>L)        error=1;
  if (error)
    {
      cout << "Invalid arguments Ephi(";
      cout <<r<<",";
      cout <<phi<<",";
      cout <<z<<",";
      cout <<r1<<",";
      cout <<phi1<<",";
      cout <<z1;
      cout <<")" << endl;;
      return 0;
    }

  double G=0;
  //cout << "Ephi = " << G << endl;
  return G;
}

