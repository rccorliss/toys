//#include "FieldSim.h"
#include "AnnularFieldSim.h"
R__LOAD_LIBRARY(.libs/libfieldsim)

void digital_current_macro_alice(){
  printf("hello\n");
  TTime now, start;
  start=now=gSystem->Now();
  printf("the time is %lu\n",(unsigned long)now);


  //load the ALICE space charge model
  const float alice_rmin=83.5;
  const float alice_rmax=254.5;
  float alice_deltar=alice_rmax-alice_rmin;
  const float alice_z=249.7;
  const float alice_driftVolt=-99930; //V
  const float alice_driftVel=2.58*1e6;//cm per s
  const float alice_chargescale=8.85e-16;//their hist. has charge in units of C/cm^3 /eps0.  This is eps0 in SI+cm units so that I can multiple by the volume in cm^3 to get the right Q.
  //get the ALICE histogram
  TFile *f=TFile::Open("InputSCDensityHistograms_8000events.root");
  TH3F* alice_average=(TH3F*)f->Get("inputSCDensity3D_8000_avg");
  now=gSystem->Now();
  printf("loaded hist.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  AnnularFieldSim *alice=
    new AnnularFieldSim(alice_rmin,alice_rmax,alice_z,9,120,9,alice_driftVel);
 AnnularFieldSim *alice2=
    new AnnularFieldSim(alice_rmin,alice_rmax,alice_z,9,90,9,alice_driftVel);
   
    // dropping half-res for test: new AnnularFieldSim(alice_rmin,alice_rmax,alice_z,53,18,31,alice_driftVel);
    //full resolution is too big:  new AnnularFieldSim(alice_rmin,alice_rmax,alice_z,159,360,62,alice_driftVel);
  now=gSystem->Now();
  printf("created sim obj.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  alice->setFlatFields(0.5, alice_driftVolt/alice_z);
  alice2->setFlatFields(0.5, alice_driftVolt/alice_z);
  now=gSystem->Now();
  printf("set fields.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  alice->load_spacecharge(alice_average,0,alice_chargescale);
  alice2->load_spacecharge(alice_average,0,alice_chargescale);
  now=gSystem->Now();
  printf("loaded spacecharge.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  alice->populate_lookup();
  alice2->populate_lookup();
  now=gSystem->Now();
  printf("populated lookup.  the dtime is %lu\n",(unsigned long)(now-start));
  start=now;
  alice->populate_fieldmap();
  alice2->populate_fieldmap();

  //define a grid of test points:
  const int nparticles=100*100;
  TVector3 testparticle[nparticles];
  TVector3 outparticle[nparticles];
  TVector3 outparticle2[nparticles];
  TVector3 backparticle[nparticles];
  float outx[nparticles], outy[nparticles],outz[nparticles];
  for (int i=0;i<nparticles;i++){
    int rpos=i/100;
    testparticle[i].SetXYZ(alice_deltar*(rpos/100.0)+alice_rmin,0,0);
    testparticle[i].RotateZ(6.28/100.0*(i%100));
  }
  for (int i=0;i<nparticles;i++){
    //printf("test[%d]=(%f,%f,%f)\n",i,testparticle[i].X(),testparticle[i].Y(),testparticle[i].Z());
    outparticle[i]=alice->swimToInSteps(alice_z,testparticle[i],600,true);
    outx[i]=outparticle[i].X();
    outy[i]=outparticle[i].Y();
    outz[i]=outparticle[i].Z();
    outparticle2[i]=alice2->swimToInSteps(alice_z,testparticle[i],600,true);
    //printf("out[%d]=(%f,%f,%f)\n",i,outparticle[i].X(),outparticle[i].Y(),outparticle[i].Z());
    //backparticle[i]=alice->swimToInSteps(testparticle[i].Z(),outparticle[i],600,true);
  }

  //  TH2F* hDeltaForward=new TH2F("hDeltaForward","xy offset of all testpoints;x offset (cm); y offset (cm)",100,-50,50,100,-50,50);
  TH2F* hDeltaForward=new TH2F("hDeltaForward","r offset vs thrown r of all testpoints;r test (cm); r offset (cm)",100,0,350,100,-40,20);
   TH2F* hForward=new TH2F("hForward","xy position of all drifted testpoints;x  (cm); y  (cm)",100,-500,500,100,-500,500);
   TH2F* hBackTest=new TH2F("hBackTest","xy residual after drift and return;x  (cm); y  (cm)",500,-0.04,0.04,500,-0.04,0.04);
    TH2F* hTest=new TH2F("hTest","xy position of all original testpoints;x  (cm); y  (cm)",100,-500,500,100,-500,500);
    TH2F* hLargeDelta=new TH2F("hLargeDelta","xy position of testpoints that have large offsets post-drift;x  (cm); y  (cm)",100,-500,500,100,-500,500);
    TH2F* hSmallDelta=new TH2F("hSmallDelta","xy position of testpoints that have small offsets post-drift;x  (cm); y  (cm)",100,-500,500,100,-500,500);
    TH1F* hDeltaMag=new TH1F("hDeltaMag","magnitude of offset after drift and return;mag (cm)",500,-1,1);
    TH2F*hDeltaR=new TH2F("hDeltaR"," drift distortion vs radius;radius (cm);offset (cm)",100,0,300,100,0,50);
    TH2F*hDeltaR2=new TH2F("hDeltaR2"," reco distortion vs radius;radius (cm);offset (cm)",100,0,300,200,-1,1);
    TH2F* h90_120=new TH2F("h90_120","drifted position offset in 90 vs 120 phi segments;x (cm);y (cm)",400,-1,1,400,-1,1);
    TVector3 dtemp,dtemp2;
  for (int i=0;i<nparticles;i++){
    dtemp=outparticle[i]-testparticle[i];
    dtemp2=outparticle[i]-outparticle2[i];
    //dtemp2=backparticle[i]-testparticle[i];
    hDeltaForward->Fill(testparticle[i].Perp(),outparticle[i].Perp()-testparticle[i].Perp());
    hForward->Fill((outparticle[i]).X(),(outparticle[i]).Y());
    hTest->Fill((testparticle[i]).X(),(testparticle[i]).Y());
    hDeltaMag->Fill(dtemp2.Perp());
    if (dtemp2.Perp()>0.05){
      hLargeDelta->Fill(testparticle[i].X(),testparticle[i].Y());
    } else {
      hSmallDelta->Fill(testparticle[i].X(),testparticle[i].Y());
    }
    hDeltaR->Fill(testparticle[i].Perp(),dtemp.Perp());
    //hDeltaR2->Fill(testparticle[i].Perp(),backparticle[i].Perp()-testparticle[i].Perp());
    hDeltaR2->Fill(testparticle[i].Perp(),dtemp2.Perp());
    h90_120->Fill(dtemp2.X(),dtemp2.Y());

    
  }
  hBackTest->Draw("colz");
  new TCanvas();
  hDeltaMag->Draw("colz");
  new TCanvas();
  // hTest->Draw("colz");
  //new TCanvas();
  h90_120->Draw("colz");

  return;
  

  AnnularFieldSim *test;
  AnnularFieldSim *testgen;
  const TVector3 cyl(60,0,100);
  TVector3 ptemp(12.005,45.005,75.99);
  TVector3 ftemp,btemp;
  TTime t[10];
  Int_t isteps=22;
  int imin=8;//can't do interpolation without at least some divisions, otherwise we necessarily go out of bounds.
  int imax=imin+isteps;
  Int_t tCreate[isteps];
  Int_t tSetField[isteps];
  Int_t tLookup[isteps];
  Int_t tGenMap[isteps];
  Int_t tSwim1k[isteps];
  Int_t steps[isteps];
  float rdiff[isteps];
  float rphidiff[isteps];
  
  testgen=new AnnularFieldSim(cyl.Perp(),2*TMath::Pi(),cyl.Z(),isteps+imin-1,isteps+imin-1,isteps+imin-1,100.0/12e-6);
  testgen->setFlatFields(1.4,200);
  testgen->populate_lookup();//2-3
  testgen->q->Set(isteps+imin-2,0,0,1e-12);///that's 10^7 protons in that box.
  testgen->populate_fieldmap();
  ftemp=ptemp;
  for (int j=0;j<isteps+imin;j++){
    ftemp=testgen->swimTo(ptemp.Z()*(isteps+imin-1-j)/(isteps+imin),ftemp,true);
  }  
  for (int i=0;i<isteps;i++){
    printf("create %d\n",i);
    t[0]=gSystem->Now();
    test=new AnnularFieldSim(cyl.Perp(),2*TMath::Pi(),cyl.Z(),i+imin,i+imin,i+imin,100.0/12e-6);
    t[1]=gSystem->Now();
    printf("setFlat %d\n",i);
    test->setFlatFields(1.4,200);
    t[2]=gSystem->Now();
    printf("populate %d\n",i);
    test->populate_lookup();//2-3
    t[3]=gSystem->Now();
    printf("setQ %d\n",i);
    test->q->Set(i+imin-1,0,0,1e-12);///that's 10^7 protons in that box.
    t[4]=gSystem->Now();
    printf("make fieldmap %d\n",i);
   test->populate_fieldmap();
    t[5]=gSystem->Now();
    printf("swim %d\n",i);
    for (int n=0;n<1000;n++){
    btemp=ftemp;
    for (int j=0;j<i+imin;j++){
      btemp=test->swimTo(ptemp.Z()*((j+1)*1.0/(1.0*(i+imin))),btemp,true);
    }
    }
    t[6]=gSystem->Now();
    tCreate[i]=(long)(t[1]-t[0]);
    tSetField[i]=(long)(t[2]-t[1]);
    tLookup[i]=(long)(t[3]-t[2]);
    tGenMap[i]=(long)(t[5]-t[4]);
    tSwim1k[i]=(long)(t[6]-t[5]);
    steps[i]=i+imin;
    rdiff[i]=(btemp.Perp()-ptemp.Perp())*1e4;
    rphidiff[i]=(btemp.Perp()*btemp.Phi()-ptemp.Perp()*ptemp.Phi())*1e4;

  }

  int i=0;

  TGraph *rPhiDiff[3];
  TMultiGraph *rphi=new TMultiGraph();
  rphi->SetTitle("r and r*phi offset for varying reco grid sizes;r(um);rphi(um)");
  rPhiDiff[i]=new TGraph(isteps,rdiff,rphidiff);
  rPhiDiff[i]->SetTitle("rphi diff;size;t(ms)");
  rPhiDiff[i]->SetMarkerColor(kBlack);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraph(1,rdiff,rphidiff);
  rPhiDiff[i]->SetTitle(Form("rphi diff %d^3 grid;size;t(ms)",steps[0]));
  rPhiDiff[i]->SetMarkerColor(kRed);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraph(1,rdiff+isteps-1,rphidiff+isteps-1);
  rPhiDiff[i]->SetTitle(Form("rphi diff %d^3 grid;size;t(ms)",steps[isteps-1]));
  rPhiDiff[i]->SetMarkerColor(kGreen);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rphi->Draw("AC*");
  
 




 i=0;
  
  TMultiGraph *multi=new TMultiGraph();
  multi->SetTitle("timing for various simulation steps;size;t(ms)");
  TGraph *graph[10];
  graph[i]=new TGraph(isteps,steps,tCreate);
  graph[i]->SetTitle("FieldSim::FieldSim;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
  graph[i]=new TGraph(isteps,steps,tSetField);
  graph[i]->SetTitle("FieldSim::SetField;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
  graph[i]=new TGraph(isteps,steps,tLookup);
  graph[i]->SetTitle("FieldSim::GenLookupTable;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
   graph[i]=new TGraph(isteps,steps,tGenMap);
  graph[i]->SetTitle("FieldSim::GenFieldMap;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
   graph[i]=new TGraph(isteps,steps,tSwim1k);
  graph[i]->SetTitle("FieldSim::Swim*1k;size;t(ms)");
  graph[i]->SetMarkerStyle(kStar);
  graph[i]->SetMarkerColor(i+1);
  multi->Add(graph[i]);
  i++;
  multi->Draw("AC*");  
  return;
}
