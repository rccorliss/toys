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
void showRosseggerDivergencePhi();
void showRosseggerDivergenceR();


void studyRosseggerBehavior(){

  // proveRosseggersAreEqual(); return;
  //sumRosseggerPhiLoops(); return;
  //showRosseggerDivergenceR(); return;
  
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
  double pb[]={0,4*acos(0)};
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
	if(1 || k!=phibin){
	  float p1=pb[0]+(pb[1]-pb[0])/loopsteps*(k+0.5);
	  double opt=ro->Ephi(r,p,z,r1,p1,z1);
	  hLoopSum->Fill(r1,z1,opt/loopsteps);
	}
      }
    }
    
  }
  hLoopSum->Draw("colz");
  return;
}


void showRosseggerDivergencePhi(){
  double r=35, r1g=45;
  double z=45, z1g=55;
  int squaresteps=20;
  int loopsteps=1000;
  float rb[]={20,78};
  double pb[]={0,4*acos(0)};
  float zb[]={0,105.5};

  int phibin=loopsteps/2;
  double p=pb[0]+(pb[1]-pb[0])/loopsteps*(phibin+0.5);
  double r1,z1;//for the graph.

  double ephi[loopsteps], phi1[loopsteps];
  double er[loopsteps], ez[loopsteps];
  //do a detailed loop over one point
  Rossegger *ro=new Rossegger(20,78,105.5);
  for (int k=0;k<loopsteps;k++){
    phi1[k]=pb[0]+(pb[1]-pb[0])/loopsteps*(k+0.5);
    ephi[k]=ro->Ephi(r,p,z,r1g,phi1[k],z1g);
    ez[k]=ro->Er(r,p,z,r1g,phi1[k],z1g);
    er[k]=ro->Ez(r,p,z,r1g,phi1[k],z1g);
  }

  //show the same-phi contribution at all points in a slice:
  TH1F *hSame[3];
  hSame[0]=new TH1F("hSamePhi",Form("Ephi component from same phi positions for all rphiz sources to rz=(%2.1f,%2.1f)",r,z),100,-0.1,0.1);
  hSame[1]=new TH1F("hSameR",Form("Er component from same phi positions for all rphiz sources to rz=(%2.1f,%2.1f)",r,z),100,-0.1,0.1);
  hSame[2]=new TH1F("hSameZ",Form("Ez component from same phi positions for all rphiz sources to rz=(%2.1f,%2.1f)",r,z),100,-0.1,0.1);
  TH2F *hSameMap[3];
  hSameMap[0]=new TH2F("hSameMapPhi",Form("Ephi component at (r%2.1f,%2.2f,%2.1f) for all source rz;r;z",r,p,z),
		       squaresteps,rb[0],rb[1],squaresteps,zb[0],zb[1]);
  hSameMap[1]=new TH2F("hSameMapR",Form("Er component at (r%2.1f,%2.2f,%2.1f) for all source rz;r;z",r,p,z),
		       squaresteps,rb[0],rb[1],squaresteps,zb[0],zb[1]);
  hSameMap[2]=new TH2F("hSameMapZ",Form("Ez component at (r%2.1f,%2.2f,%2.1f) for all source rz;r;z",r,p,z),
		       squaresteps,rb[0],rb[1],squaresteps,zb[0],zb[1]);
  for (int i=0;i<squaresteps;i++){
    r1=rb[0]+(rb[1]-rb[0])/squaresteps*(i+0.5);
    for (int j=0;j<squaresteps;j++){
      z1=zb[0]+(zb[1]-zb[0])/squaresteps*(j+0.5);
      for (int k=0;k<1;k++){//squaresteps;k++){
	float p1=pb[0]+(pb[1]-pb[0])/squaresteps*(k+0.5);
	float tempephi=ro->Ephi(r,p1,z,r1,p1,z1);
	hSame[0]->Fill(tempephi);
	hSameMap[1]->Fill(r1,z1,ro->Er(r,p1,z,r1,p1,z1));
	hSameMap[2]->Fill(r1,z1,ro->Ez(r,p1,z,r1,p1,z1));
	hSameMap[0]->Fill(r1,z1,tempephi);
      }
    }
  }


  TCanvas *c=new TCanvas("cshowRosseggerDivergence","showRosseggerDivergence",800,600);
  c->Divide(3,2);
  
  c->cd(1);
  TGraph *gPhiComp=new TGraph(loopsteps,phi1,ephi);
  gPhiComp->SetTitle(Form("E_phi(phi1) at (r%2.1f,%2.2f,%2.1f) from %d-step loops in phi at rz (r%2.1f,phi,z%2.1f);source phi (field calc'd at %2.4f=phi);field value (arb)",r,p,z,loopsteps,r1g,z1g,p));
  gPhiComp->Draw();
  TGraph *gPhiPoint=new TGraph(1,&(phi1[phibin]),&(ephi[phibin]));
  gPhiPoint->SetMarkerColor(kRed);
  gPhiPoint->Draw("*");
  c->cd(2);
  TGraph *gRComp=new TGraph(loopsteps,phi1,er);
  gRComp->SetTitle(Form("E_r(phi1) at (r%2.1f,%2.2f,%2.1f) from %d-step loops in phi at rz=(r%2.1f,phi,z%2.1f);source phi (field calc'd at %2.4f=phi);field value (arb)",r,p,z,loopsteps,r1g,z1g,p));
  gRComp->Draw();
  TGraph *gRPoint=new TGraph(1,&(phi1[phibin]),&(er[phibin]));
  gRPoint->SetMarkerColor(kRed);
  gRPoint->Draw("*");
  c->cd(3);
  TGraph *gZComp=new TGraph(loopsteps,phi1,ez);
  gZComp->SetTitle(Form("E_z(phi1) at (r%2.1f,%2.2f,%2.1f) from %d-step loops in phi at r-z coordinates (r%2.1f,phi,z%2.1f);source phi (field calc'd at %2.4f=phi);field value (arb)",r,p,z,loopsteps,r1g,z1g,p));
  gZComp->Draw();
  TGraph *gZPoint=new TGraph(1,&(phi1[phibin]),&(ez[phibin]));
  gZPoint->SetMarkerColor(kRed);
  gZPoint->Draw("*");
  c->cd(4);
  hSameMap[0]->Draw("colz");
 c->cd(5);
  hSameMap[1]->Draw("colz");
 c->cd(6);
  hSameMap[2]->Draw("colz");


  
  return;
}

void showRosseggerDivergenceR(){
 
  int histsteps=20;
  int graphsteps=1000;
  float rb[]={20,78};
  double pb[]={0,4*acos(0)};
  float zb[]={0,105.5};

 double r=35, r1=45;
  double z=45, z1=55;
  int pbin=4;
  double p=pb[0]+(pb[1]-pb[0])/histsteps*(pbin+0.5);
 int p1bin=8;
  double p1=pb[0]+(pb[1]-pb[0])/graphsteps*(p1bin+0.5);

  
  int xbin=graphsteps/2;
  r=rb[0]+(rb[1]-rb[0])/graphsteps*(xbin+0.5);

  double ephi[graphsteps], r1g[graphsteps];
  double er[graphsteps], ez[graphsteps];
  //do a detailed loop over one point
  Rossegger *ro=new Rossegger(20,78,105.5);
  for (int k=0;k<graphsteps;k++){
    r1g[k]=rb[0]+(rb[1]-rb[0])/graphsteps*(k+0.5);
    ephi[k]=ro->Ephi(r,p,z,r1g[k],p1,z1);
    ez[k]=ro->Er(r,p,z,r1g[k],p1,z1);
    er[k]=ro->Ez(r,p,z,r1g[k],p1,z1);
  }

  //show the same-phi contribution at all points in a slice:
  TH1F *hSame[3];
  hSame[0]=new TH1F("hSamePhi",Form("Ephi component from same r positions for all rphiz sources to pz=(%2.1f,%2.1f)",p,z),100,-0.1,0.1);
  hSame[1]=new TH1F("hSameR",Form("Er component from same r positions for all rphiz sources to pz=(%2.1f,%2.1f)",p,z),100,-0.1,0.1);
  hSame[2]=new TH1F("hSameZ",Form("Ez component from same r positions for all rphiz sources to pz=(%2.1f,%2.1f)",p,z),100,-0.1,0.1);
  TH2F *hSameMap[3];
  hSameMap[0]=new TH2F("hSameMapPhi",Form("Ephi component at (r%2.1f,%2.2f,%2.1f) for all source phiz;r;z",r,p,z),
		       histsteps,pb[0],pb[1],histsteps,zb[0],zb[1]);
  hSameMap[1]=new TH2F("hSameMapR",Form("Er component at (r%2.1f,%2.2f,%2.1f) for all source phiz;r;z",r,p,z),
		       histsteps,pb[0],pb[1],histsteps,zb[0],zb[1]);
  hSameMap[2]=new TH2F("hSameMapZ",Form("Ez component at (r%2.1f,%2.2f,%2.1f) for all source phiz;r;z",r,p,z),
		       histsteps,pb[0],pb[1],histsteps,zb[0],zb[1]);
  for (int i=0;i<1;i++){//histsteps;i++){
    r1=rb[0]+(rb[1]-rb[0])/histsteps*(i+0.5);
    for (int j=0;j<histsteps;j++){
      z1=zb[0]+(zb[1]-zb[0])/histsteps*(j+0.5);
      for (int k=0;k<histsteps;k++){
	p1=pb[0]+(pb[1]-pb[0])/histsteps*(k+0.5);
	float tempephi=ro->Ephi(r,p1,z,r1,p1,z1);
	hSame[0]->Fill(tempephi);
	hSameMap[1]->Fill(r1,z1,ro->Er(r1,p,z,r1,p1,z1));
	hSameMap[2]->Fill(r1,z1,ro->Ez(r1,p,z,r1,p1,z1));
	hSameMap[0]->Fill(r1,z1,tempephi);
      }
    }
  }


  TCanvas *c=new TCanvas("cshowRosseggerDivergence","showRosseggerDivergence",800,600);
  c->Divide(3,2);
  
  c->cd(1);
  TGraph *gPhiComp=new TGraph(graphsteps,r1g,ephi);
  gPhiComp->SetTitle(Form("E_phi(r1) at (r%2.1f,%2.2f,%2.1f) from %d steps in r at phiz (p%2.1f,z%2.1f);source r (field calc'd at %2.4f=r);field value (arb)",r,p,z,graphsteps,p1,z1,r));
  gPhiComp->Draw();
  TGraph *gPhiPoint=new TGraph(1,&(r1g[xbin]),&(ephi[xbin]));
  gPhiPoint->SetMarkerColor(kRed);
  gPhiPoint->Draw("*");
  c->cd(2);
  TGraph *gRComp=new TGraph(graphsteps,r1g,er);
  gRComp->SetTitle(Form("E_r(r1) at (r%2.1f,%2.2f,%2.1f) from %d steps in r at phiz (p%2.1f,z%2.1f);source r (field calc'd at %2.4f=r);field value (arb)",r,p,z,graphsteps,p1,z1,r));
  gRComp->Draw();
  TGraph *gRPoint=new TGraph(1,&(r1g[xbin]),&(er[xbin]));
  gRPoint->SetMarkerColor(kRed);
  gRPoint->Draw("*");
  c->cd(3);
  TGraph *gZComp=new TGraph(graphsteps,r1g,ez);
  gZComp->SetTitle(Form("E_z(r1) at (r%2.1f,%2.2f,%2.1f) from %d steps in r at phiz (p%2.1f,z%2.1f);source r (field calc'd at %2.4f=r);field value (arb)",r,p,z,graphsteps,p1,z1,r));
  gZComp->Draw();
  TGraph *gZPoint=new TGraph(1,&(r1g[xbin]),&(ez[xbin]));
  gZPoint->SetMarkerColor(kRed);
  gZPoint->Draw("*");
  c->cd(4);
  hSameMap[0]->Draw("colz");
 c->cd(5);
  hSameMap[1]->Draw("colz");
 c->cd(6);
  hSameMap[2]->Draw("colz");


  
  return;
}
