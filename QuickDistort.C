
void QuickDistort(){

  //load the distortionmap
  //TFile *dFile=TFile::Open("elevatorpitch/fluct_average.rev3.1side.3d.file0.h_negz.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  TFile *dFile=TFile::Open("elevatorpitch/fluct_single.1side.3d.file0.h_Charge_0.real_B1.4_E-400.0.ross_phi1_sphenix_phislice_lookup_r26xp40xz40.distortion_map.hist.root","READ");
  TH3F *hDistortX=(TH3F*)dFile->Get("hIntDistortionX");
  TH3F *hDistortY=(TH3F*)dFile->Get("hIntDistortionY");
  //there's also a z component, but we'll ignore that for now.  It's small.
  
  //load the chargemap for comparison
  //TFile *cFile=TFile::Open("averages/average.rev3.hist.root","READ");
  TFile *cFile=TFile::Open("testfluctcharge/fluct_outputFile_15kHz_G4Hits_sHijing_0-12fm_000000_012000_bX1508071_bias0.root","READ");
  TH3D* hQ3D=(TH3D*)cFile->Get("h_Charge_0");
  TH2D* hQ2D=(TH2D*)hQ3D->Project3D("yx");
  

  //eventually we should probe at every stripe, but for now we'll go once per r-phi bin in the charge map:
  float phibound[2];
  float rbound[2];
  int nphi, nr;
  float phistep,rstep;
  phibound[0]=hQ2D->GetXaxis()->GetXmin();
  phibound[1]=hQ2D->GetXaxis()->GetXmax();
  nphi=hQ2D->GetXaxis()->GetNbins();
  phistep=(phibound[1]-phibound[0])/(nphi*1.0);
  rbound[0]=hQ2D->GetYaxis()->GetXmin();
  rbound[1]=hQ2D->GetYaxis()->GetXmax();
  nr=hQ2D->GetYaxis()->GetNbins();
  rstep=(rbound[1]-rbound[0])/(nr*1.0);
  TH2D* hCmDistortion=new TH2D("hCmDistortion","Radial Shift (r_f-r_i) of CM hits;phi (rad);r (m)",
			       nphi,phibound[0],phibound[1],
			       nr,rbound[0],rbound[1]);
  TH2D* hChargeLow=new TH2D("hChargeLow","Lower-Res Charge Fluctuations;phi (rad);r (m)",
			       nphi/5,phibound[0],phibound[1],
			       nr/5,rbound[0],rbound[1]);
  TH2D* hSanityCheck=new TH2D("hSanityCheck","Sanity Check of how many times each cell is touched;phi (rad);r (m)",
			       nphi,phibound[0],phibound[1],
			       nr,rbound[0],rbound[1]);
  TH2D* hRatio=new TH2D("hRatio","#Delta R / Qtot in the bin;phi (rad);r (m)",
			       nphi,phibound[0],phibound[1],
			       nr,rbound[0],rbound[1]);
  /* for sets that aren't average-subtracted:
  TH2D* hComparison=new TH2D("hComparison","#Delta R vs Qtot in the bin;#Delta R (cm);Qtot (ions/bin)",
			     200,-0.1,1.1,
			     200,0,2e7);
  */
  TH2D* hComparison=new TH2D("hComparison","#Delta R-ave vs Qtot-ave in the bin;#Delta R (cm);Qtot (ions/bin)",
			     200,-0.004,0.006,
			     200,-1e6,3e6);
  TH2D* hDistComp=new TH2D("hDistComp","Fluctuation #Delta R at cm and halfway;#Delta R@0 (um);#Delta R@50 (um)",
			     100,-50,50,
			     100,-50,50);
  vector<float> rdist,qint,rpos,phipos;//radial distortion and charge integral per bin.
  vector<float> rdisthalf;//radial distortion at z=50cm

  TVector3 posBefore, posAfter;
  for (int i=0;i<nr;i++){
    float r=rbound[0]+rstep*(i+0.5);
    posBefore.SetXYZ(r,0,0);//neglect z coordinate for now.
    for (int j=0;j<nphi;j++){
      float phi=phibound[0]+phistep*(j+0.5);
      posBefore.SetPhi(phi);
      //printf("checking phi=%f, r=%f\n",phi,r);
      float deltaX=hDistortX->Interpolate(phi,r*100.0, 1.0);//map is in cm, but charge is in m, so need to convert r here.  set z=1.0cm
      float deltaY=hDistortY->Interpolate(phi,r*100.0,1.0);//map is in cm, but charge is in m, so need to convert r here.
      posAfter.SetX(posBefore.X()+deltaX);
      posAfter.SetY(posBefore.Y()+deltaY);
      float deltaR=posAfter.Perp()-posBefore.Perp();
      float qIntegral=hQ2D->GetBinContent(hQ2D->FindBin(phi,r));

      deltaX=hDistortX->Interpolate(phi,r*100.0, 50.0);//map is in cm, but charge is in m, so need to convert r here.  set z=1.0cm
      deltaY=hDistortY->Interpolate(phi,r*100.0,50.0);//map is in cm, but charge is in m, so need to convert r here.
      posAfter.SetX(posBefore.X()+deltaX);
      posAfter.SetY(posBefore.Y()+deltaY);
      float deltaR2=posAfter.Perp()-posBefore.Perp();
      
      hCmDistortion->Fill(phi,r,deltaR);
      hChargeLow->Fill(phi,r,qIntegral);
      hSanityCheck->Fill(phi,r,1);
      //hRatio->Fill(phi,r,abs(deltaR)/qIntegral);
      hRatio->Fill(phi,r,deltaR/qIntegral);
      hComparison->Fill(deltaR,qIntegral);
      hDistComp->Fill(deltaR*1e4,deltaR2*1e4);
      rdist.push_back(deltaR);
      rdisthalf.push_back(deltaR2);
      qint.push_back(qIntegral);
      rpos.push_back(r);
      phipos.push_back(phi);
    }
  }

  TCanvas *c=new TCanvas("c","distortion comparison to charge",800,800);
  c->Divide(2,2);
  c->cd(1);
  hChargeLow->Draw("colz");
  c->cd(2);
  hCmDistortion->Draw("colz");
  c->cd(3);//->SetLogz();
  TGraph *g;
  g=new TGraph(rdist.size(),&(rdist[0]),&(rdisthalf[0]));
  g->SetTitle("Distortion at z=0 vs at z=50cm;dist@0;dist@50");
  //g->Draw("*A");
  hDistComp->Draw("colz");
  /*
  g=new TGraph(rdist.size(),&(rdist[0]),&(qint[0]));
  g->SetTitle("Distortion vs integrated charge in column;dist;q");
  g->Draw("*A");
  */
  //hRatio->Draw("colz");
  c->cd(4);
  //this looked fine, but if you mess with stuff, look again: hSanityCheck->Draw("colz");
  g=new TGraph(rdist.size(),&(rpos[0]),&(rdist[0]));
  //g->SetTitle("Distortion vs radial position;r (m);#Delta R");
  //g->Draw("*A");
  hComparison->Draw("colz");
  //hqluusinqqtac iie
   

  
  return;
}
