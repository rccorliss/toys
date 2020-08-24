
void HeuristicCharge(){

  int kNAzimuthalSteps=360;
  int kNRadialSteps=159;
  int nphi=360;
  int kNLongitudinalSteps=62*2;
  float z_rdo=105.5;
  float rmin=20;
  float rmax=78;
  
  //these settings are for 50kHz Au+Au collisions?
  float e0 = 8.854187817e-12; //[C]/[Vm]
  float gas = 1.0/76628.0; //[Vs] -- something about the ionization per unit volume and the time it takes to clear the gas, in units I don't grok.
  float mult = 400.0; //950.0;
  float rate = 5e+4; //[1/s]
  double a=mult*rate*e0*gas; // C/m;
  float b=1.0/z_rdo; //[1/cm]
  float c=2.0/3.0*20.0;
  a=a*1e15; //fC/m
  a=a/100;//fC/cm
  //rho= a/r^2*(1-bz+c)
  //so c is the gain times the backflow percentage -- the number of ions per electron.
  //in the plot I am trying to match, we have gain = 2000 and IBF=0.3% --> c=2000*0.3/100=6, which is actually half of the listed value.
  
  const double proton_charge=1.6e-4;//proton charge in fC
  double a_ions=a/proton_charge;//ion number density per cm?


  TFile *outfile=TFile::Open("HeuristicSc_test.root","RECREATE");

  TH3D *hCharge=new TH3D("sphenix_minbias_average",Form("Heuristic SC per cm^3 (a=%2.2E, b=%2.2E,c=%2.2E);phi (rad);r (cm);z (cm)",a_ions,b,c),kNAzimuthalSteps,0,6.28319,kNRadialSteps,rmin,rmax,kNLongitudinalSteps,0,z_rdo);
  
  printf("a_ions = %f\n",a_ions);
  printf("b = %f\n",b);
  printf("c = %f\n",c);
  for(int r=0; r!=kNRadialSteps; ++r) {
    float dr = hCharge->GetYaxis()->GetBinCenter( r+1 ); //[cm]
    for(int p=0; p!=kNAzimuthalSteps; ++p) {
      float dp = hCharge->GetXaxis()->GetBinCenter( p+1 );
      for(int z=0; z!=kNLongitudinalSteps; ++z) {
	float dz = hCharge->GetZaxis()->GetBinCenter( z+1 ); //[cm]
	float dRho = a_ions*(1-b*TMath::Abs(dz)+c)/(dr*dr); //ions/cm^3 
	hCharge->SetBinContent(p+1,r+1,z+1,dRho);
	//if(fDebug>2) printf("@{Ir,Ip,Iz}={%d (%f),%d (%f),%d (%f)}, rho %f\n",r,dr,p,dp,z,dz,dRho);
      }
    }
  }
  hCharge->Write();
  outfile->Close();
  return;
}
