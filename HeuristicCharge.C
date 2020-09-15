
struct model_params{
  //detector string:
  char name[100];
  //geometry
  float z;
  float rmin;
  float rmax;

  //gas parameters for a*(1-b*TMath::Abs(dz)+c)/(pow(dr,d))
  float a;//to get charge per volume
  float aprime; //to get current per area on a surface
  float b;
  float c;
  float d;
};

model_params Sphenix2020ModelParams(){
    float z_rdo=105.5;
  float rmin=20;
  float rmax=78;

  
  float e0 = 8.854187817e-12; //[C]/[Vm]

  //these settings are for Ne:CF4 50:50
  
  //note: gas [V*s]=empirical_factor*ions_per_cm*drift_time
  float empirical=1.0/431538.0;//[V*cm] empirical factor derived from ALICE's derivation based on STAR.  See http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-017.pdf
  float empirical_alt=1.0/17094.0;//[V*cm/keV] empirical factor derived from ALICE's derivation based on STAR.  See http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-017.pdf
  
  float ion_mobility=1.65;//[cm/s*cm/V] ion mobility in 50:50 Ne:CF4 gas, not easily computed from first principles.
  float drift_v=400;//[V/cm] drift voltage
  float ion_velocity= ion_mobility*drift_v;// ion velocity in [cm/s]
  float drift_time=z_rdo/ion_velocity;//[s] time charge deposited at z=outer_edge will be in the tpc.

  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm

  float ions_per_cm=0.50 * Ne_NTotal+0.50 * CF4_NTotal;//[ions/cm] approximate ions per cm.  71.5 in 50:50
  float kev_per_cm=0.50 * Ne_dEdx + 0.50 * CF4_dEdx;///[kev/cm] mip energy loss in 4.28 in 50:50

  //gas=empirical*ions_per_cm*drift_time;//[V*s] More plausible scaling argument, but smaller overall yield than alt.
  float gas_alt=empirical_alt*kev_per_cm*drift_time;//[V*s]
  double mult = 450.0; //950.0;
  double rate = 5e+4; //[1/s]
  double a=mult*rate*gas_alt*e0; // C/m;
  float b=1.0/z_rdo; //[1/cm]

  double gain=1410;//[#] electrons on rdo per electron arriving. originally 2000, scaled so that total gain in the old 90:10 is approximately the same as total gain here.;
  double ibf_per_gain=0.0067;//[#] 0.67% IBF, 2 ions escaping back per every 300 electrons in the shower.
   float c=gain*ibf_per_gain;//[#] ions flowing back per electron arriving.  
  a=a*1e15; //fC/m
  a=a/100;//fC/cm
  double aprime=a*ion_velocity*1e-6;//[nC/s=nA]

  //rho= a/r^2*(1-bz+c)

  printf("ion_velocity=%E cm/s\n",ion_velocity);
  printf("drift_time=%E s\n",drift_time);
  printf("gas_alt=%E V*s\n",gas_alt);
  printf("a=%E fC/cm\n",a);
  
  model_params m;
  sprintf(m.name,"sPHENIX2020");
  m.z=z_rdo;
  m.rmin=rmin;
  m.rmax=rmax;
  m.a=a;
  m.aprime=aprime;
  m.b=b;
  m.c=c;
  m.d=2;
  return m;
}

model_params Sphenix2020ElectronsOnRDOSpecialParams(){
  float z_rdo=105.5;//cm
  float rmin=20;//cm
  float rmax=78;//cm

  
  float e0 = 8.854187817e-12; //[C]/[Vm]

  //these settings are for Ne:CF4 50:50
  
  //note: gas [V*s]=empirical_factor*ions_per_cm*drift_time
  float empirical=1.0/431538.0;//[V*cm] empirical factor derived from ALICE's derivation based on STAR.  See http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-017.pdf
  float empirical_alt=1.0/17094.0;//[V*cm/keV] empirical factor derived from ALICE's derivation based on STAR.  See http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-017.pdf
  
  // float ion_mobility=1.65;//[cm/s*cm/V] ion mobility in 50:50 Ne:CF4 gas, not easily computed from first principles.
  // float drift_v=400;//[V/cm] drift voltage
  float ele_velocity= 8e6;// electron velocity in [cm/s] at 400V.
  float drift_time=z_rdo/ele_velocity;//[s] time charge deposited at z=outer_edge will be in the tpc.

  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm

  float ions_per_cm=0.50 * Ne_NTotal+0.50 * CF4_NTotal;//[ions/cm] approximate ions per cm.  71.5 in 50:50
  float kev_per_cm=0.50 * Ne_dEdx + 0.50 * CF4_dEdx;///[kev/cm] mip energy loss in 4.28 in 50:50

  //gas=empirical*ions_per_cm*drift_time;//[V*s] More plausible scaling argument, but smaller overall yield than alt.
  float gas_alt=empirical_alt*kev_per_cm*drift_time;//[V*s]
  double mult = 400.0; //950.0;
  double rate = 1e+5; //[1/s]
  double a=mult*rate*gas_alt*e0; // C/m;
  float b=1.0/z_rdo; //[1/cm]

  //all together:
  a=rate*drift_time*mult*empirical_alt*e0*kev_per_cm; //primaries per cm in the radial direction.
  
  double gain=2000;//[#] electrons on rdo per electron arriving. originally 2000, scaled so that total gain in the old 90:10 is approximately the same as total gain here.;
  double ibf_per_gain=0.0067;//[#] 0.67% IBF, 2 ions escaping back per every 300 electrons in the shower.
   float c=gain;//[#] electrons on rdo per electron arriving.
  a=a*1e15; //fC/m
  a=a/100;//fC/cm
  //rho= a/r^2*(1-bz+c) -- For current, we will generally evaluate this at z=0 or z=rdo.


  double aprime=a*ele_velocity*1e-6;//[fC/s*1e-6=nA]


  printf("ele_velocity=%E cm/s\n",ele_velocity);
  printf("f_E=%E C/cm\n",empirical_alt*e0*kev_per_cm/100);
  printf("drift_time=%E s\n",drift_time);
  printf("gas_alt=%E V*s\n",gas_alt);
  printf("a=%E fC/cm\n",a);
  printf("a'=%E nA\n",aprime);
  
  model_params m;
  sprintf(m.name,"sPHENIX2020electrons");
  m.z=z_rdo;
  m.rmin=rmin;
  m.rmax=rmax;
  m.a=a;
  m.aprime=aprime;
  m.b=b;
  m.c=c;
  m.d=2;
  return m;
}

model_params Sphenix2020RccTestModelParams(){
    float z_rdo=105.5;
  float rmin=20;
  float rmax=78;

  
  float e0 = 8.854187817e-12; //[C]/[Vm]

  //these settings are for Ne:CF4 50:50
  
  //note:  gas*eps0 [Cs/m]=[C/m]*[s]=pseudorapidity_averaged_primaries_per_cm*drift_length/velocity
  float ion_mobility=1.65;//[cm/s*cm/V];
  float drift_v=400;//[V/cm] drift voltage
  float ion_velocity= ion_mobility*drift_v;// ion velocity in [cm/s]
  float drift_time=z_rdo/ion_velocity;//[s] time charge deposited at z=outer_edge will be in the tpc.
  float ions_per_cm=71.5;//[ions/cm] approximate ions per radial cm -- ions per charge * (1+sqrt(2))/2
  float proton_charge=1.6e-4;//[C/ion] proton charge in fC

  float gas_e0=ions_per_cm*proton_charge*drift_time;//[fC*s/cm]
  printf("gas_e0 calc: %1.9f=(1/%5.2f)\n",gas_e0,1/gas_e0);

  //annoying numerology:
  // 1/(  (50*1.2)*(1.6e-4)*(105.5/(3.3*400))/8.85e-12 )= 1.15e-8
  // 1/76628*8.85e-12 = 1.15e-16
  //in order to get this right, I really want to match the term that Carlos was using.
  //this suggests that 1/gas=ions_per_cm*proton_charge*(z_rdo/ion_velocity)/eps0
  //but this would imply that the gas factor goes /down/ as ions_per_cm or charge_per_ion goes up, which is clearly wrong.  Similarly, we want gas factor to go up as ion_velocity goes down, which this does not present.

  double mult = 400.0; //950.0;
  double rate = 5e+4; //[1/s]
  double a=mult*rate*gas_e0; // C/m;
  float b=1.0/z_rdo; //[1/cm]
   float c=(2.0/3.0*1/100)*2000.0;//gain of 2k, 0.67% IBF per unit gain.
  a=a*1e15; //fC/m
  a=a/100;//fC/cm

  double aprime=a*ion_velocity*1e-6;//[nC/s=nA]
  //rho= a/r^2*(1-bz+c)
  
  model_params m;
  sprintf(m.name,"sPHENIX2020RccTest");
  m.z=z_rdo;
  m.rmin=rmin;
  m.rmax=rmax;
  m.a=a;
  m.b=b;
  m.c=c;
  m.d=2;
  return m;
}


model_params Sphenix2018ModelParams(){
    float z_rdo=105.5;
  float rmin=20;
  float rmax=78;

  //these settings are for 50kHz Au+Au collisions?
  float e0 = 8.854187817e-12; //[C]/[Vm]
  float gas = 1.0/76628.0; //[Vs] -- something about the ionization per unit volume and the time it takes to clear the gas, in units I don't grok.  This is from a fit STAR performs, scaled to ALICE values in http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-017.pdf,
  
  double mult = 450.0; //950.0;
  double rate = 5e+4; //[1/s]
  double a=mult*rate*e0*gas; // C/m;
  float b=1.0/z_rdo; //[1/cm]
   float c=(0.3/100.0)*2000.0;//gain of 2k, 0.3% IBF per unit gain.
  a=a*1e15; //fC/m
  a=a/100;//fC/cm

  //don't have this available.   double aprime=a*ion_velocity*1e-6;//[nC/s=nA]

  //rho= a/r^2*(1-bz+c)
  //so c is the gain times the backflow percentage -- the number of ions per electron.
  //in the plot I am trying to match, we have gain = 2000 and IBF=0.3% --> c=2000*0.3/100=6, which is actually half of the listed value.
  model_params m;
  sprintf(m.name,"sPHENIX2018");
  m.z=z_rdo;
  m.rmin=rmin;
  m.rmax=rmax;
  m.a=a;
  m.aprime=0;
  m.b=b;
  m.c=c;
  m.d=2;
  return m;
}

model_params AliceModelParams(){
  float z_rdo=249.7;
  float rmin=83.5;
  float rmax=254.5;
  float ion_velocity=z_rdo/0.156;//[cm/s]
  float e0 = 8.854187817e-12; //[C]/[Vm]
  float gas = 1.0/76628.0; //[Vs] -- something about the ionization per unit volume and the time it takes to clear the gas, in units I don't grok.  
  double mult = 900.0; //950.0;
  double rate = 5e+4; //[1/s]
  double a=mult*rate*e0*gas; // C/m;
  float b=1.0/z_rdo; //[1/cm]
   float c=(1.0/100.0)*2000.0;//gain of 2k, 1% IBF per unit gain.
  a=a*1e15; //fC/m
  a=a/100;//fC/cm

  double aprime=a*ion_velocity*1e-6;//[nC/s=nA]
  model_params m;
    sprintf(m.name,"ALICE");
  m.z=z_rdo;
  m.rmin=rmin;
  m.rmax=rmax;
  m.a=a;
  m.aprime=aprime;
  m.b=b;
  m.c=c;
  m.d=2;
  return m;
}

model_params AliceAltModelParams(){
  float z_rdo=249.7;
  float rmin=83.5;
  float rmax=254.5;
  float ion_velocity=z_rdo/0.156;//[cm/s]
  float e0 = 8.854187817e-12; //[C]/[Vm]
  float gas = 1.0/76628.0; //[Vs] -- something about the ionization per unit volume and the time it takes to clear the gas, in units I don't grok.  
  double mult = 900.0; //950.0;
  double rate = 5e+4; //[1/s]
  double a=mult*rate*e0*gas; // C/m;
  float b=1.0/z_rdo; //[1/cm]
  float c=10;//per ALICE description of 'eps 10'.
  a=a*1e15; //fC/m
  a=a/100;//fC/cm

  double aprime=a*ion_velocity*1e-6;//[nC/s=nA]
  model_params m;
    sprintf(m.name,"ALICEeps10");
  m.z=z_rdo;
  m.rmin=rmin;
  m.rmax=rmax;
  m.a=a;
  m.aprime=aprime;
  m.b=b;
  m.c=c;
  m.d=2;
  return m;
}


double CalcCharge(model_params m,double r,double z){
  return m.a*(1-m.b*TMath::Abs(z)+m.c)/pow(r,m.d);
}
double CalcCurrent(model_params m,double r,double z){
  return m.aprime*(1-m.b*TMath::Abs(z)+m.c)/pow(r,m.d);
}

void HeuristicCharge(){

  int kNAzimuthalSteps=360;
  int kNRadialSteps=159;
  int kNLongitudinalSteps=62*2;

  //model_params par=Sphenix2020ModelParams();
  model_params par=Sphenix2020ElectronsOnRDOSpecialParams();
  //model_params par=Sphenix2018ModelParams();
  //model_params par=AliceModelParams();


  TFile *outfile=TFile::Open(Form("HeuristicSc_%s.root",par.name),"RECREATE");

  TH3D *hCharge=new TH3D("heuristic",Form("Heuristic Spacecharge (fC per cm^3) (a=%2.2E, b=%2.2E,c=%2.2E, d=%2.1f);phi (rad);r (cm);z (cm)",par.a,par.b,par.c, par.d),kNAzimuthalSteps,0,6.28319,kNRadialSteps,par.rmin,par.rmax,kNLongitudinalSteps,0,par.z);
  
  printf("a = %f\n",par.a);
  printf("b = %f\n",par.b);
  printf("c = %f\n",par.c);
  for(int r=0; r!=kNRadialSteps; ++r) {
    float dr = hCharge->GetYaxis()->GetBinCenter( r+1 ); //[cm]
    for(int p=0; p!=kNAzimuthalSteps; ++p) {
      float dp = hCharge->GetXaxis()->GetBinCenter( p+1 );
      for(int z=0; z!=kNLongitudinalSteps; ++z) {
	float dz = hCharge->GetZaxis()->GetBinCenter( z+1 ); //[cm]
	float dRho = CalcCharge(par,dr,dz); //ions/cm^3 
	hCharge->SetBinContent(p+1,r+1,z+1,dRho);
	//if(fDebug>2) printf("@{Ir,Ip,Iz}={%d (%f),%d (%f),%d (%f)}, rho %f\n",r,dr,p,dp,z,dz,dRho);
      }
    }
  }

  TH2D *hCurrent=new TH2D("current",Form("Heuristic Spacecharge Current(nA per cm^2) @CM (a'=%2.2E, b=%2.2E,c=%2.2E, d=%2.1f);phi (rad);r (cm)",par.aprime,par.b,par.c, par.d),kNAzimuthalSteps,0,6.28319,kNRadialSteps,par.rmin,par.rmax);
  for(int r=0; r!=kNRadialSteps; ++r) {
    float dr = hCharge->GetYaxis()->GetBinCenter( r+1 ); //[cm]
    for(int p=0; p!=kNAzimuthalSteps; ++p) {
      float dp = hCharge->GetXaxis()->GetBinCenter( p+1 );
      float dJ = CalcCurrent(par,dr,0); //ions/cm^3 
      hCurrent->SetBinContent(p+1,r+1,dJ);
      //if(fDebug>2) printf("@{Ir,Ip,Iz}={%d (%f),%d (%f),%d (%f)}, rho %f\n",r,dr,p,dp,z,dz,dRho);
    }
    
  }
  hCharge->Write();
  hCurrent->Write();
  outfile->Close();
  printf("Wrote HeuristicSc_%s.root\n",par.name);
  return;
}
