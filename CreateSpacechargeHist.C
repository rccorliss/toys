#if 1
#include "/osx_sphenix/coresoftware/simulation/g4simulation/g4main/PHG4HitContainer.h"
#include "/osx_sphenix/coresoftware/simulation/g4simulation/g4main/PHG4Hit.h"
#else
#include <g4main/PHG4HitContainer.h>
#include <g4main/PHG4Hit.h>
#endif
R__LOAD_LIBRARY(libg4testbench.so)
R__LOAD_LIBRARY(libphg4hit.so)

//need to include something to read DST.TPC.G4HIT_TPC
//and maybe 'G4HIT_ABSORBER_TPC' as well?

void CreateSpacechargeHist(const char *dirname, const char *filename, int istart=0, bool saveTree=false){
  printf("are you running with the sphenix env?  This probably doesn't work without that!\n");

  gSystem->Load("libg4testbench.so");
  gSystem->Load("libphg4hit.so");

  
  // TFile *file=TFile::Open("/osx_sphenix/G4Hits_sHijing_9-11fm_09900_10000.root","R");
  TFile *file=TFile::Open(Form("%s%s",dirname,filename),"R");
  assert(file.IsOpen());
  TTree *T=(TTree*)file->Get("T");
  int neve=T->GetEntries();

  
  PHG4HitContainer *eleHits=new PHG4HitContainer();
  printf("eleHits pointer is %p\n",(void *)eleHits);
  T->SetBranchAddress("DST.TPC.G4HIT_TPC",&eleHits);
  printf("after set-address, eleHits pointer is %p\n",(void *)eleHits);




  float us=1.0,ms=us*1e3,s=ms*1e3;
  float um=1e-4, mm=um*1e3, cm=mm*10,m=mm*1e3; //changed to make 'cm' 1.0, for convenience.
  float Hz=1/s, kHz=1/ms, MHz=1/us;
  //used two ways:  1) to apply units to variables when defined
  //                2) to divide by certain units so that those variables are expressed in those units.
  
  float vIon=16.0*um/us;
  float ampGain=1e4;
  float ampIBFfrac=0.01;
  float ionsPerEle=ampGain*ampIBFfrac;
  float mbRate=22*kHz;
  float z_rdo=105.5*cm;
  float rmin=20*cm;
  float rmax=78*cm;

  double Ne_dEdx = 1.56;   // keV/cm
  double CF4_dEdx = 7.00;  // keV/cm
  double Ne_NTotal = 43;    // Number/cm
  double CF4_NTotal = 100;  // Number/cm
  double Tpc_NTot = 0.90 * Ne_NTotal + 0.10 * CF4_NTotal;
  double Tpc_dEdx = 0.90 * Ne_dEdx + 0.10 * CF4_dEdx;
  double Tpc_ElectronsPerKeV = Tpc_NTot / Tpc_dEdx;
  double Tpc_ElectronsPerGeV = Tpc_NTot / Tpc_dEdx*1e6 /(cm); //input as electrons per gev per cm.

  
  TFile *outfile=TFile::Open(Form("%s.rcc_sc.hist.root",filename),"RECREATE");
  int nr=159;
  int nphi=360;
  int nz=110;
  TH3D *hCharge=new TH3D("sphenix_minbias_charge","SC (ions) per cm^3;phi (rad);r (cm);z (cm)",nphi,0,6.28319,nr,rmin/cm,rmax/cm,110,0,z_rdo/cm);
  double hrstep=(rmax-rmin)/cm/nr;
  double hphistep=6.28319/nphi;
  double hzstep=z_rdo/cm/nz;
  TTree *rawHits;


  //to fully populate the detector half, we need events covering t0=0 to t0=z_rdo/vIon.  These occur at mbRate, so:
  printf("vIon=%f cm/s\tz=%f cm\t rate=%f kHz\n ==> need z/v*r=%f events to cover the detector. T has %d entries.\n",vIon/(cm/s),z_rdo/(cm),mbRate/(kHz),(z_rdo/vIon*mbRate), neve);
  
  float x,y,z,zibf,zprim;
  float r,phi;
  float ne;

  float driftedZ;
  int testi=5;
  int i;
  if (saveTree){
    rawHits=new TTree("hTree","tpc hit tree for ionization");
    rawHits->Branch("x",&x);
    rawHits->Branch("y",&y);
    rawHits->Branch("zorig",&z);
    rawHits->Branch("zibf",&zibf);
    rawHits->Branch("zprim",&zprim);
    rawHits->Branch("r",&r);
    rawHits->Branch("phi",&phi);
    rawHits->Branch("ev",&i);
    rawHits->Branch("ne",&ne);
  }

  
  T->GetEntry(testi);
  printf("after GetEntry(%d), eleHits pointer is %p\n",testi,(void *)eleHits);

  //return;

  for (i=1;i<neve;i++){
    T->GetEntry(i);
    //eventually I ought to roll a random number for each bunch crossing and use that to detemrine the number of events, but for now I just take one per minbiasRate.
    //must be sure to convert the incoming positions into local units.
    //load prim*
    //printf("loading eve=%d\n",i);
    float t0=(i+istart)/mbRate;
    driftedZ=t0*vIon;//drift position in local units
    PHG4HitContainer::ConstRange range=eleHits->getHits();
    assert(range);

    float f=0.5;;//for now, just pick the middle of the hit.  Do better later.

    
    for (PHG4HitContainer::ConstIterator hiter=range.first;hiter!=range.second;hiter++)
      {
	ne=hiter->second->get_eion()*Tpc_ElectronsPerGeV;
	//load the three coordinates in with units of cm.
	x = (hiter->second->get_x(0) + f * (hiter->second->get_x(1) - hiter->second->get_x(0)))*(cm);
	y = (hiter->second->get_y(0) + f * (hiter->second->get_y(1) - hiter->second->get_y(0)))*(cm);
	z = (hiter->second->get_z(0) + f * (hiter->second->get_z(1) - hiter->second->get_z(0)))*(cm);
	if (z<0) continue;
	r=sqrt(x*x+y*y);
	phi=atan2(x,y);
	zprim=z-driftedZ;
	zibf=z_rdo-driftedZ;
	if (phi<0) phi+=6.28319;
	//compute the bin volume:
	int bin=hCharge->FindBin(phi,r/(cm),zprim/(cm));
	double hr=hCharge->GetYaxis()->GetBinLowEdge(bin);
	double vol=(hzstep*hphistep*(hr+hrstep*0.5)*hrstep)/cm/cm/cm;

	hCharge->Fill(phi,r/(cm),zprim/(cm),ne/vol); //primary ion, drifted by t0, in cm
	hCharge->Fill(phi,r/(cm),zibf/(cm),ne*ionsPerEle/vol); //amp ion, drifted by t0, in cm
	if (saveTree){
	  rawHits->Fill();
	}
    }

    
  }
  outfile->cd();
  hCharge->Write();
  if (saveTree){
    rawHits->Write();
  }
  outfile->Close();
  return;
}