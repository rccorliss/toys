void PoorManSmooth(){
  TFile *source=TFile::Open("BeamXingNBeams.root","R");
  TFile *output=TFile::Open("Smooth.temp.50kHz.root","RECREATE");

  //load the source histograms:
  TH3F *hIbfIn=(TH3F*)(source->Get("sphenix_minbias_IBF"));
  TH3F *hPrimaryIn=(TH3F*)(source->Get("sphenix_minbias_primary"));
TH3F *hTotalIn=(TH3F*)(source->Get("sphenix_minbias_charge"));
  //compute their dimensions:
  int np=hPrimaryIn->GetNbinsX();
  int nr=hPrimaryIn->GetNbinsY();
  int nz=hPrimaryIn->GetNbinsZ();
  float p1=hPrimaryIn->GetXaxis()->GetXmax();
  float p0=hPrimaryIn->GetXaxis()->GetXmin();
  float r1=hPrimaryIn->GetYaxis()->GetXmax();
  float r0=hPrimaryIn->GetYaxis()->GetXmin();
  float z1=hPrimaryIn->GetZaxis()->GetXmax();
  float z0=hPrimaryIn->GetZaxis()->GetXmin();

  //get the projections:
  TH2F *hPrimaryRZ=(TH2F*)(hPrimaryIn->Project3D("YZ"));
  TH2F *hPrimaryPR=(TH2F*)(hPrimaryIn->Project3D("XY"));
  TH2F *hIbfRZ=(TH2F*)(hIbfIn->Project3D("YZ"));
  TH2F *hIbfPR=(TH2F*)(hIbfIn->Project3D("XY"));
  TH1F *hIbfP=(TH1F*)(hIbfIn->Project3D("X"));
  TH1F *hIbfR=(TH1F*)(hIbfIn->Project3D("Y"));

  //make smooth histograms:
  TH3F *hSmoothIbf=new TH3F("hSmoothIbf","Smoothed IBF dist;phi;r;z",np,p0,p1,nr,r0,r1,nz,z0,z1);
  TH3F *hSmoothPrimary=new TH3F("hSmoothPrimary","Smoothed Primary dist;phi;r;z",np,p0,p1,nr,r0,r1,nz,z0,z1);
  TH3F *hSmoothTotal=new TH3F("sphenix_minbias_average","Smoothed Total SC dist;phi;r;z",np,p0,p1,nr,r0,r1,nz,z0,z1);
  TH3F *hFluct=new TH3F("sphenix_minbias_fluct","SC minus average SC dist;phi;r;z",np,p0,p1,nr,r0,r1,nz,z0,z1);

  hSmoothIbf->Fill(0.,0.,0.,0.);
  hSmoothPrimary->Fill(0.,0.,0.,0.);
  hFluct->Fill(0.,0.,0.,0.);
  hSmoothTotal->Fill(0.,0.,0.,0.);
  for (int i=0;i<np;i++){
    for (int j=0;j<nr;j++){
      for (int k=0;k<nz;k++){
	float ibfcontent=hIbfR->GetBinContent(j+1)/np/nz;
	float primcontent=hPrimaryRZ->GetBinContent(k+1,j+1)/np;
	hSmoothIbf->SetBinContent(i+1,j+1,k+1,ibfcontent);
	hSmoothPrimary->SetBinContent(i+1,j+1,k+1,primcontent);
	hSmoothTotal->SetBinContent(i+1,j+1,k+1,primcontent+ibfcontent);
	hFluct->SetBinContent(i+1,j+1,k+1,hTotalIn->GetBinContent(i+1,j+1,k+1)-primcontent-ibfcontent);
      }
    }
  }

  //get the smoothed projections:
  TH2F *sPrimaryRZ=(TH2F*)(hSmoothPrimary->Project3D("YZ"));
  sPrimaryRZ->SetTitle("Smooth Primary RZ");
  TH2F *sPrimaryPR=(TH2F*)(hSmoothPrimary->Project3D("XY"));
  sPrimaryPR->SetTitle("Smooth Primary PhiR");
  TH2F *sIbfRZ=(TH2F*)(hSmoothIbf->Project3D("YZ"));
  sIbfRZ->SetTitle("Smooth Ibf RZ");
  TH2F *sIbfPR=(TH2F*)(hSmoothIbf->Project3D("XY"));
  sIbfPR->SetTitle("Smooth Ibf PhiR");

  
  TCanvas *c=new TCanvas("c","smoothing charge data",1200,800);
  c->Divide(4,3);
  c->cd(1);
  hPrimaryRZ->Draw("colz");
  c->cd(2);
  hPrimaryPR->Draw("colz");
  c->cd(3);
  hIbfRZ->Draw("colz");
  c->cd(4);
  hIbfPR->Draw("colz");
  c->cd(5);
  sPrimaryRZ->Draw("colz");
  c->cd(6);
  sPrimaryPR->Draw("colz");
  c->cd(7);
  sIbfRZ->Draw("colz");
  c->cd(8);
  sIbfPR->Draw("colz");
  c->cd(9);
  ((TH2F*)(hSmoothTotal->Project3D("YZ")))->Draw("colz");
  c->cd(10);
  ((TH2F*)(hSmoothTotal->Project3D("XY")))->Draw("colz");
   c->cd(11);
  ((TH2F*)(hFluct->Project3D("YZ")))->Draw("colz");
  c->cd(12);
  ((TH2F*)(hFluct->Project3D("XY")))->Draw("colz");
  
  /*
  for (int i=0;i<hDistortR[3]->GetNbinsX();i++){
    for (int j=0;j<hDistortR[3]->GetNbinsY();j++){
      if (hDistortC[0]->GetBinContent(i+1,j+1)>99 && hDistortC[1]->GetBinContent(i+1,j+1)>99){
      hProfileSum[0]->Fill(hDistortR[3]->GetBinContent(i+1,j+1));
      hProfileSum[1]->Fill(hDistortP[3]->GetBinContent(i+1,j+1));
      }
    }
  }
  */
  //source->Close();
  c->SaveAs("PoorManSmooth.Summary.pdf");
 output->cd();
 hSmoothTotal->Write();
 hSmoothPrimary->Write();
 hSmoothIbf->Write();
 hFluct->Write();
 output->Close();
  return;
}
