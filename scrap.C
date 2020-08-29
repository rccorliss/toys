//#include "FieldSim.h"

/*
digital_current_macro_alice started as code to model the alice-specific TPC properties, to match to work the ALICE group did.  
it turned out to be a convenient point to compare sPHENIX performance, so it became ill-named.

This code loads an AnnularFieldSim model of the TPC, computing spacecharge distortions and using those to propagate particles from arbitrary points in the volume to a specified z-plane.  It assumes propagation or back-propagation based on the sign of (destination-start).  Details of the simulation are set in the first section, and details of the particles to be propagated in the second.

 */


#include "AnnularFieldSim.h"
#include "Rossegger.h"

R__LOAD_LIBRARY(.libs/libfieldsim)


void proveRosseggersAreEqual();
void sumRosseggerPhiLoops();

void scrap(){

  // proveRosseggersAreEqual(); return;
  sumRosseggerPhiLoops(); return;
  
  const float tpc_rmin=20.0;
  const float tpc_rmax=78.0;
  const float tpc_z=105.5;
  const float tpc_driftVolt=-400*tpc_z; //V =V_CM-V_RDO -- volts per cm times the length of the drift volume.
  const float tpc_driftVel=8.0*1e6;//cm per s  -- 2019 nominal value
  const float tpc_magField=1.4;//T -- 2019 nominal value
  const double tpc_chargescale=-1.6e-19;  //multiple charge hist contents by this to get coulombs.
 
  int nr=4;//10;//24;//159;//159 nominal
  int nr_roi_min=0;
  int nr_roi=nr;//10;
  int nr_roi_max=nr_roi_min+nr_roi;
  int nphi=4;//38;//360;//360 nominal
  int nphi_roi_min=0;
  int nphi_roi=nphi;//38;
  int nphi_roi_max=nphi_roi_min+nphi_roi;
  int nz=4;//62;//62 nominal
  int nz_roi_min=0;
  int nz_roi=nz;
  int nz_roi_max=nz_roi_min+nz_roi;

  AnnularFieldSim *tpc=
    new  AnnularFieldSim(tpc_rmin,tpc_rmax,tpc_z,
			 nr, nr_roi_min,nr_roi_max,1,2,
			 nphi,nphi_roi_min, nphi_roi_max,1,2,
			 nz, nz_roi_min, nz_roi_max,1,2,
			 tpc_driftVel, AnnularFieldSim::PhiSlice, AnnularFieldSim::FromFile);
  char lookup_string[200];
  sprintf(lookup_string,"ross_phislice_lookup_r%dxp%dxz%d",nr,nphi,nz);
  char lookupFilename[200];
  sprintf(lookupFilename,"%s.root",lookup_string);

  tpc->load_rossegger();
  printf("loaded rossegger greens functions. (phi set to zero)\n");
  tpc->populate_lookup();
  tpc->save_phislice_lookup(lookupFilename);
  printf("populated lookup.\n");
  return;
}

void proveRosseggersAreEqual(){
  int squaresteps=15;
  float rb[]={20,78};
  float pb[]={0,6.28};
  TH2F *hCompare=new TH2F("hCompare","Comparison of streamlined and non-streamlined values",squaresteps,rb[0],rb[1],squaresteps,pb[0],pb[1]);
  TH2F *hCompare2=new TH2F("hCompare2","Comparison of streamlined and non-streamlined values",squaresteps,rb[0],rb[1],squaresteps,pb[0],pb[1]);
  TH1F *hFrac=new TH1F("hFrac","Ratio-1 of streamlined and non-streamlined values",100,-1e-10,1e-10);
  Rossegger *ro=new Rossegger(20,78,105.5);
  Rossegger *ro2=new Rossegger(20,78,105.5);
  for (int i=0;i<squaresteps;i++){
    float r=rb[0]+(rb[1]-rb[0])/squaresteps*(i+0.5);
    for (int j=0;j<squaresteps;j++){
      float p=pb[0]+(pb[1]-pb[0])/squaresteps*(j+0.5);
      for (int k=0;k<squaresteps;k++){
	float r1=rb[0]+(rb[1]-rb[0])/squaresteps*(k+0.5);
	for (int l=0;l<squaresteps;l++){
	  float p1=pb[0]+(pb[1]-pb[0])/squaresteps*(l+0.5);
	  double opt=ro->Er(r,p,52.75,r1,p1,40);
	  double bare=ro2->Er_(r,p,52.75,r1,p1,40);
	  double ratio=opt/bare;
	  hFrac->Fill(ratio-1);
	  if (abs(ratio-1)>1E-4)
	    hCompare2->Fill(r,p,ratio);

	  hCompare->Fill(r,p,ratio);
	}
      }
    }
  }
  hFrac->Draw();  return;
  hCompare->Draw("colz");
  return;
}

void sumRosseggerPhiLoops(){
  int squaresteps=15;
  int loopsteps=50;
  float rb[]={20,78};
  float pb[]={0,6.28};
  float zb[]={0,105.5};

  int phibin=0;
  float r=40,p=pb[0]+(pb[1]-pb[0])/loopsteps*(phibin+0.5),z=27;
  TH2F *hLoopSum=new TH2F("hLoopSum",Form("Sum of field at (r%2.1f,%2.2f,%2.1f) from %d-step loops in phi at r-z coordinates",r,p,z,loopsteps),squaresteps,rb[0],rb[1],squaresteps,zb[0],zb[1]);
  //TH2F *hCompare2=new TH2F("hCompare2","Comparison of streamlined and non-streamlined values",squaresteps,rb[0],rb[1],squaresteps,pb[0],pb[1]);
  // TH1F *hFrac=new TH1F("hFrac","Ratio-1 of streamlined and non-streamlined values",100,-1e-10,1e-10);

  Rossegger *ro=new Rossegger(20,78,105.5);
  for (int i=0;i<squaresteps;i++){
    float r1=rb[0]+(rb[1]-rb[0])/squaresteps*(i+0.5);
    for (int j=0;j<squaresteps;j++){
      float z1=zb[0]+(zb[1]-zb[0])/squaresteps*(j+0.5);
      for (int k=0;k<loopsteps;k++){
	float p1=pb[0]+(pb[1]-pb[0])/loopsteps*(k+0.5);
	double opt=ro->Ephi(r,p,z,r1,p1,z1);
	hLoopSum->Fill(r1,z1,opt/loopsteps);
      }
    }
    
  }
  hLoopSum->Draw("colz");
  return;
}
