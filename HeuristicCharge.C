
void HeuristicCharge(){

  int kNAzimuthalSteps=360;
  int kNRadialSteps=159;
  int nphi=360;
  int kNLongitudinalSteps=62*2;
  float z_rdo=105.5;
  float rmin=20;
  float rmax=78;
  
  
  float e0 = 8.854187817e-12*1e+9; //[C]/[Vm]*1e+9
  float gas = 1.0/76628.0; //[Vs]
  float mult = 400.0; //950.0;
  float rate = 5e+4; //[Hz]
  double a=mult*rate*e0*gas; // fC/cm;
  float b=100.0/z_rdo; //[1/m]
  float c=2.0/3.0*20.0;

  const double proton_charge=1.6e-4;//proton charge in fC
  double a_ions=a/proton_charge;//ion number density per cm^3


  TFile *outfile=TFile::Open("HeuristicSc.root","RECREATE");

  TH3D *hCharge=new TH3D("sphenix_minbias_average",Form("Heuristic SC per cm^3 (a=%2.2E, b=%2.2E,c=%2.2E);phi (rad);r (cm);z (cm)",a_ions,b,c),kNAzimuthalSteps,0,6.28319,kNRadialSteps,rmin,rmax,kNLongitudinalSteps,0,z_rdo);
  
  printf("a_ions = %f\n",a_ions);
  printf("b = %f\n",b);
  printf("c = %f\n",c);
  for(int r=0; r!=kNRadialSteps; ++r) {
    float dr = hCharge->GetYaxis()->GetBinCenter( r+1 )/100.0; //[m]
    for(int p=0; p!=kNAzimuthalSteps; ++p) {
      float dp = hCharge->GetXaxis()->GetBinCenter( p+1 );
      for(int z=0; z!=kNLongitudinalSteps; ++z) {
	float dz = hCharge->GetZaxis()->GetBinCenter( z+1 )/100.0; //[m]
	float dRho = a_ions*(1-b*TMath::Abs(dz)+c)/(dr*dr); //ions/cm^3 --- rccsays: but dr is in m, not cm!
	hCharge->SetBinContent(p+1,r+1,z+1,dRho);
	//if(fDebug>2) printf("@{Ir,Ip,Iz}={%d (%f),%d (%f),%d (%f)}, rho %f\n",r,dr,p,dp,z,dz,dRho);
      }
    }
  }
  hCharge->Write();
  outfile->Close();
  return;
}
