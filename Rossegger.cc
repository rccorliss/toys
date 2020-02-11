#include "Rossegger.h"

#include <iostream>
#include <fstream>
#include <math.h>
#include "TMath.h"
#include "TFile.h"
#include <string>
#include <sstream>

//#include "/usr/local/include/complex_bessel.h"
#include <boost/math/special_functions.hpp> //covers all the special functions.

#include "TH2D.h"
#include "TH3.h"

using namespace std;
using namespace TMath;
//using namespace boost::math::special_functions

//Bessel Function J_n(x):
#define jn(order,x) boost::math::cyl_bessel_j(order,x)
//Bessel (Neumann) Function Y_n(x):
#define yn(order,x) boost::math::cyl_neumann(order,x)
//Modified Bessel Function of the first kind I_n(x):
#define in(order,x) boost::math::cyl_bessel_i(order,x)
//Modified Bessel Function of the second kind K_n(x):
#define kn(order,x) boost::math::cyl_bessel_k(order,x)
#define limu(im_order,x) Rossegger::Limu(im_order,x)
#define kimu(im_order,x) Rossegger::Kimu(im_order,x)

/*
  This is a modified/renamed copy of Carlos and Tom's "Spacecharge" class, modified to use boost instead of fortran routines, and with phi terms added.
 */


Rossegger::Rossegger(double InnerRadius, double OuterRadius, double Rdo_Z)
{
  a = InnerRadius;
  b = OuterRadius;
  L = Rdo_Z;

  verbosity =1;
  pi = 2.0 * asin(1.0);
  cout << pi << endl;

 
  char limufile[]="limu_table.csv";
  char kimufile[]="kimu_table.csv";
  LoadCsvToHist(&hLimu,limufile);
  LoadCsvToHist(&hKimu,kimufile);
  cout << "hKimu is real, see:" <<  hKimu->GetXaxis()->GetNbins() << endl;
  FindMunk(0.01);
  FindBetamn(0.0001);
 

  cout << "Rossegger object initialized as follows:" << endl;
  cout << "  Inner Radius = " << a << " cm." << endl;
  cout << "  Outer Radius = " << b << " cm." << endl;
  cout << "  Half  Length = " << L << " cm." << endl;
  cout << "  Limu Dataset = " << limufile << endl;

  return ;
}

double Rossegger::FindNextZero(double xstart,double epsilon,int order, double (Rossegger::*func)(int, double)){
    
  double previous=(this->*func)(order,xstart);
  double x=xstart+epsilon;
  double value=previous;
	 
  while (! (  (value == 0) || (value<0 && previous>0) || (value>0 && previous<0)) ){
    	  //  Rossegger equation 5.12
    value = (this->*func)(order, x);
    if (value == 0) cout << "hit it exactly!  Go buy a lottery ticket!" << endl;
    if ( (value == 0) || (value<0 && previous>0) || (value>0 && previous<0))
      {
	//when we go from one sign to the other, we have bracketed the zero
	//the following is mathematically equivalent to finding the delta
	//between the x of our previous value and the point where we cross the axis
	//then returning x0=x_old+delta
	double slope = (value-previous)/epsilon;
	double intercept =  value - slope*x;
	double x0  = -intercept/slope;
	if (verbosity>1) cout << " " << x0 << "," << endl;
	double n0=  (this->*func)(order, x-epsilon);
	double n1= (this->*func)(order, x+epsilon);
	if ((n0<0 && n1<0)||(n0>0 && n1>0)){
	  printf("neighbors on both sides have the same sign.  Check your function resolution!\n");
	}
 
	return x0;
      }
    previous = value;
    x+=epsilon;
  }
  cout <<"logic break!\n";
  assert(1==2);
  return 0;

}


  
void Rossegger::FindBetamn(double epsilon)
{
  cout << "Now filling the Beta[m][n] Array..."<<endl;
  for (int m=0; m<NumberOfOrders; m++)
    {
      if (verbosity) cout << "Filling Beta["<<m<<"][n]..." << endl;

      double x=epsilon;
      for (int n=0;n<NumberOfOrders;n++){//  !!!  Off by one from Rossegger convention  !!!
	x=FindNextZero(x,epsilon,m,&Rossegger::Rmn_for_zeroes);
	Betamn[m][n]=x/b;
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
	  //N2mn[m][n] *= (jn(m,Betamn[m][n]*a)*jn(m,Betamn[m][n]*a)/jn(m,Betamn[m][n]*b)/jn(m,Betamn[m][n]*b) - 1.0);
	  double jna_over_jnb=jn(m,Betamn[m][n]*a)/jn(m,Betamn[m][n]*b);
	  N2mn[m][n] *= (jna_over_jnb*jna_over_jnb-1.0);
	  //rcc note!  in eq 5.17, N2nm is set with betamn[m][n].  this is reversed here.  Not sure if important.
	  if (verbosity>1) cout << "m: " << m << " n: " << n << " N2[m][n]: " << N2mn[m][n];
	  double integral=0.0;
	  double step = 0.01;
	  if (verbosity>1)
	    {
	      for (double r=a; r<b; r+=step){
		integral += Rmn(m,n,r)*Rmn(m,n,r)*r*step;
	      }
	      cout << " Int: " << integral << endl;
	    }
	  //N2mn[m][n] = integral;
	}
    }

  cout << "Done." << endl;
}


void Rossegger::FindMunk(double epsilon)
{
  cout << "Now filling the Mu[n][k] Array..."<<endl;
  // We're looking for the zeroes of Rossegger eqn. 5.46:
  // R_nk(mu_nk;a,b)=Limu(Beta_n*a)Kimu(Beta_n*b)-Kimu(Beta_n*a)Limu(Beta_n*b)=0
  // since a and b are fixed, R_nk is a function solely of mu_nk and n.
  // for each 'n' we wish to find the a set of k mu_n's that result in R_nk=0

  for (int n=0; n<NumberOfOrders; n++)//  !!!  Off by one from Rossegger convention  !!!
    {
      if (verbosity) cout << "Filling Mu["<<n<<"][k]..." << endl;
      double x=epsilon;
      for (int k=0;k<NumberOfOrders;k++){
	x=FindNextZero(x,epsilon,n,&Rossegger::Rnk_for_zeroes);
	Munk[n][k]=x;
	x+=epsilon;
	if (verbosity>0) {
	  printf("Mu[%d][%d]=%E\n",n,k,Munk[n][k]);	  
	  printf("adjacent values are Rnk[mu-epsilon]=%E\tRnk[mu+epsilon]=%E\n",
		 Rnk_for_zeroes(n,x-epsilon),Rnk_for_zeroes(n,x+epsilon));
	  printf("values of argument to limu and kimu are %f and %f\n",
		 (n+1)*pi/L*a,(n+1)*pi/L*b);
	}
      }
    }
  cout << "Done." << endl;
  return;
}

 double Rossegger::Limu(double mu, double x){
   //defined in Rossegger eqn 5.44, also a canonical 'satisfactory companion' to Kimu.
   //could use Griddit?
   if (verbosity>1) {
     printf("Limu::mu=%f,x=%f\n",mu,x);
     printf("limu axis bounds= mu:(%f to %f) x:(%f to %f)\n",
	    hLimu->GetXaxis()->GetXmin(),hLimu->GetXaxis()->GetXmax(),
	    hLimu->GetYaxis()->GetXmin(),hLimu->GetYaxis()->GetXmax());
   }
   int xbin=hLimu->GetXaxis()->FindBin(mu);

   if (xbin<1 || xbin>hLimu->GetXaxis()->GetNbins()) {
     printf("Limu::mu=%f,x=%f is out of bounds\n",mu,x);
     printf("Limu axis bounds= mu:(%f to %f) x:(%f to %f)\n",
	    hLimu->GetXaxis()->GetXmin(),hLimu->GetXaxis()->GetXmax(),
	    hLimu->GetYaxis()->GetXmin(),hLimu->GetYaxis()->GetXmax());
     assert(1==3);
   }
   int ybin=hLimu->GetYaxis()->FindBin(x);
   if (ybin<1 || ybin>hLimu->GetYaxis()->GetNbins()) {
     printf("Limu::mu=%f,x=%f is out of bounds\n",mu,x);
     printf("Limu axis bounds= mu:(%f to %f) x:(%f to %f)\n",
	    hLimu->GetXaxis()->GetXmin(),hLimu->GetXaxis()->GetXmax(),
	    hLimu->GetYaxis()->GetXmin(),hLimu->GetYaxis()->GetXmax());
     assert(1==3);
   }

   return hLimu->GetBinContent(xbin,ybin); //this should really be interpolating, but we'll deal with that later.
 }
 double Rossegger::Kimu(double mu, double x){
   //could use Griddit?
   if (verbosity>1) {
     printf("Kimu::mu=%f,x=%f\n",mu,x);
   }
   int xbin=hKimu->GetXaxis()->FindBin(mu);
   if (verbosity>1) {
	printf("Kimu:: ybin=%d\n",xbin);
   }
   if (xbin<1 || xbin>hKimu->GetXaxis()->GetNbins()) {
     printf("Kimu::mu=%f,x=%f is out of bounds\n",mu,x);
     printf("Kimu axis bounds= mu:(%f to %f) x:(%f to %f)\n",
	    hKimu->GetXaxis()->GetXmin(),hKimu->GetXaxis()->GetXmax(),
	    hKimu->GetYaxis()->GetXmin(),hKimu->GetYaxis()->GetXmax());
     assert(1==3);
   }
   int ybin=hKimu->GetYaxis()->FindBin(x);
   if (ybin<1 || ybin>hKimu->GetYaxis()->GetNbins()) {
     printf("Kimu::mu=%f,x=%f is out of bounds\n",mu,x);
     printf("Kimu axis bounds= mu:(%f to %f) x:(%f to %f)\n",
	    hKimu->GetXaxis()->GetXmin(),hKimu->GetXaxis()->GetXmax(),
	    hKimu->GetYaxis()->GetXmin(),hKimu->GetYaxis()->GetXmax());
     assert(1==3);
   }
   return hKimu->GetBinContent(xbin,ybin); //this should really be interpolating, but we'll deal with that later.
 }

 double Rossegger::Rmn_for_zeroes(int m, double x){
   double lx = a*x/b;
   return jn(m,x)*yn(m,lx) - jn(m,lx)*yn(m,x);
 }

double Rossegger::Rmn(int m, int n, double r){
  if (verbosity>1) cout << "Determine Rmn("<<m<<","<<n<<","<<r<<") = ";

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

  //  Calculate the function using C-libraries from boost
  //  Rossegger Equation 5.11:
  //         Rmn(r) = Ym(Beta_mn a)*Jm(Beta_mn r) - Jm(Beta_mn a)*Ym(Beta_mn r)
  double R=0;
  R = yn(m,Betamn[m][n]*a)*jn(m,Betamn[m][n]*r) - jn(m,Betamn[m][n]*a)*yn(m,Betamn[m][n]*r);

  if (verbosity>1) cout << R << endl;
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
  R = kn(m,BetaN*a)*in(m,BetaN*r)-in(m,BetaN*a)*kn(m,BetaN*r);

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
  R = kn(m,BetaN*b)*in(m,BetaN*r)-in(m,BetaN*b)*kn(m,BetaN*r);

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
  //  NOTE:  K-m(z) = Km(z) and I-m(z) = Im(z)... though boost handles negative orders.
  //
  // with: s -> ref,  t -> r, 
  //  NOTE:  BetaN is defined near Rossegger Equation 5.27, 5.43, and other places: BetaN= n*pi / L
  //  to match our change in definition of order 0 in " !!!  Off by one from Rossegger convention  !!!", we move n->n+1
  double BetaN = (n+1)*pi/L;
  double term1 = kn(m,BetaN*ref)*( in(m-1,BetaN*r) + in(m+1,BetaN*r) );
  double term2 = in(m,BetaN*ref)*( kn(m-1,BetaN*r) + kn(m+1,BetaN*r) );
  R = BetaN/2.0*(term1 + term2);

  return R;
}

double Rossegger::Rnk_for_zeroes(int n, double mu){
  if (verbosity>1) printf("Rnk_for_zeroes called with n=%d,mu=%f\n",n,mu);
  double BetaN=(n+1)*pi/L;
    //  Rossegger Equation 5.46
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN b) - Kimu_nk(BetaN a) Limu_nk (BetaN b)
  
  return limu(mu,BetaN*a)*kimu(mu,BetaN*b)- kimu(mu,BetaN*a)*limu(mu,BetaN*b);
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
  double BetaN=(n+1)*pi/L;
   //  Rossegger Equation 5.45
  //       Rnk(r) = Limu_nk (BetaN a) Kimu_nk (BetaN r) - Kimu_nk(BetaN a) Limu_nk (BetaN r)

  return limu(Munk[n][k],BetaN*a)*kimu(Munk[n][k],BetaN*r)- kimu(Munk[n][k],BetaN*a)*limu(Munk[n][k],BetaN*r);

}


double Rossegger::Ez(double r, double phi, double z, double r1, double phi1, double z1)
{
  //if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileEZ(r,phi,z,r1,phi1,z1);
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
  //field at r, phi, z due to unit charge at r1, phi1, z1;
  //if(fByFile && fabs(r-r1)>MinimumDR && fabs(z-z1)>MinimumDZ) return ByFileER(r,phi,z,r1,phi1,z1);
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


void Rossegger::LoadCsvToHist(TH2** hist, char* sourcefile){
  cout << "trying to load " << sourcefile << endl;
  std::ifstream inf(sourcefile);
  string line;
  string token;
  printf("opened.\n");
  printf("status? %d\n",inf.is_open()?1:0);
  //    inf>>line;
  //  printf("%s\n",line.c_str());

  double val[3];//={-12345,-12345,-12345};
  double mu, x, f;
  
  //delete hist; // get rid of it if we already had it around.

  //read the histogram bounds from the file:
  string axisname[2];
  double firstelement[2];
  double lastelement[2];
  int nsteps[2];
  double stepsize[2];
  double lbound[2];
  double ubound[2];
  for (int i=0;(i<2 && !inf.eof());i++){
    inf>>line;
    cout << line << endl;
    //printf("%s\n",line.c_str());
    stringstream lineparser(line);
    getline(lineparser,token,',');
    axisname[i]=token;
    getline(lineparser,token,',');
    firstelement[i]=std::stod(token);
    getline(lineparser,token,',');
    lastelement[i]=std::stod(token);
    getline(lineparser,token,',');
    nsteps[i]=std::stoi(token);
    stepsize[i]=(lastelement[i]-firstelement[i])/(nsteps[i]-1);
    lbound[i]=firstelement[i]-stepsize[i]/2;
    ubound[i]=lastelement[i]+stepsize[i]/2;
    printf("loaded axis=%s,limits=%f,%f,steps=%d from file\n",axisname[i].c_str(),firstelement[i],lastelement[i],nsteps[i]);

  }

  

  
  *hist=new TH2D(sourcefile,Form("%s;%s;%s",sourcefile,axisname[0].c_str(),axisname[1].c_str()),
		nsteps[0],lbound[0],ubound[0],
		nsteps[1],lbound[1],ubound[1]);
  printf("hist has name %s, xaxis bins=%d\n",(*hist)->GetName(),(*hist)->GetXaxis()->GetNbins());
  int nlines=0;
  while (!inf.eof()){
    inf>>line;//mu>>x>>val;
    //printf("%s\n",line.c_str());
    stringstream lineparser(line);
    for (int i=0;i<3;i++){
      getline(lineparser,token,',');
      //printf("token=%s\n",token.c_str());
      val[i]=std::stod(token);
      //printf("double=%f\n",val[i]);
    }
    mu=val[0];
    x=val[1];
    f=val[2];
    //printf("loaded %f,%f,%f from file\n",mu,x,f);
    (*hist)->Fill(mu,x,f);
    nlines++;
    //if(nlines>100) assert(1==2);
  }
    
  return;
}
