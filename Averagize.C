//make a smooth average out of multiple components

//void Averagize(char *sourceFileName="evgeny_apr/Summary_Average_AA_events.root", float scale=0.08012, char *ibfName="h_IBFCharge_evt_0", char *primaryName="h_PrimCharge_evt_0"){
void Averagize(char *sourceFileName="evgeny_jan_2022.sum.hist.root", float scale=0.0417, char *ibfName="_h_SC_ibf_0", char *primaryName="_h_SC_prim_0"){
  //scale is because there were 24 input files.
  printf("default sumw2=%d\n",TH1::GetDefaultSumw2());
  //load files
  TFile *source=TFile::Open(sourceFileName,"READ");
  TH3D *hSourceIBF=(TH3D*)(source->Get(ibfName));
  TH3D *hSourcePrimary=(TH3D*)(source->Get(primaryName));

  //set up output:
  // TFile *output=TFile::Open("newAverage.hist.root","RECREATE"); //should make this an argument.
  TFile *output=TFile::Open("newAverage.Apr.2022.hist.root","RECREATE"); //should make this an argument.  Changed to this in Apr 2022
  TH3D *hOutputIBF=new TH3D(*hSourceIBF);
  hOutputIBF->Reset();
  hOutputIBF->SetTitle("Number of IBF (ions) (smoothed)");
  hOutputIBF->Sumw2(false);
  TH3D *hOutputPrimary=new TH3D(*hSourcePrimary);
  hOutputPrimary->Reset();
  hOutputPrimary->SetTitle("Number of Primary (ions) (smoothed)");
  printf("sanity check:  outputname=%s, inputname=%s\n",hOutputPrimary->GetTitle(), hSourcePrimary->GetTitle());
  hOutputPrimary->Sumw2(false);

  //make smooth IBF by averaging over z IN EACH HALF SEPARATELY
  //find where to divide, and how many bins we have on each side:
  int nZbins=hSourceIBF->GetZaxis()->GetNbins();
  int midZbin=  hSourceIBF->GetZaxis()->FindBin(0.0);

  //negative side:
  int nNegBins=midZbin-1;
  hSourceIBF->GetZaxis()->SetRange(1,midZbin-1);
  TH2D *hNeg=(TH2D*)(hSourceIBF->Project3D("yx"));
  for (int r=0;r<hSourceIBF->GetYaxis()->GetNbins();r++){
    for (int p=0;p<hSourceIBF->GetXaxis()->GetNbins();p++){
      double aveIBF=hNeg->GetBinContent(p+1,r+1)/(nNegBins*1.0);      
      for (int z=0;z<nNegBins;z++){
	hOutputIBF->SetBinContent(p+1,r+1,z+1,aveIBF);
      }
    }
  }
  
    //positive side:
  int nPosBins=nZbins-midZbin-1;
  hSourceIBF->GetZaxis()->SetRange(midZbin+1,nZbins-1);
  TH2D *hPos=(TH2D*)(hSourceIBF->Project3D("yx"));
  for (int r=0;r<hSourceIBF->GetYaxis()->GetNbins();r++){
    for (int p=0;p<hSourceIBF->GetXaxis()->GetNbins();p++){
      double aveIBF=hPos->GetBinContent(p+1,r+1)/(nPosBins*1.0);      
      for (int z=midZbin;z<nZbins;z++){
	hOutputIBF->SetBinContent(p+1,r+1,z+1,aveIBF);
      }
    }
  }
  hSourceIBF->GetZaxis()->SetRange();


  //intermediate region? -- this is a reminder that eventually we need to SEPARATE the two halves and treat them independently.
  //
  for (int r=0;r<hSourceIBF->GetYaxis()->GetNbins();r++){
    for (int p=0;p<hSourceIBF->GetXaxis()->GetNbins();p++){
      double aveIBF=0.5*(hOutputIBF->GetBinContent(p+1,r+1,midZbin-1)+hOutputIBF->GetBinContent(p+1,r+1,midZbin+1));
      hOutputIBF->SetBinContent(p+1,r+1,midZbin,aveIBF);
      
    }
  }  
  //

  
  

  //make smooth prim by averaging over phi
  //find where to divide, and how many bins we have on each side:
  int nPbins=hSourcePrimary->GetXaxis()->GetNbins();

  TH2D *hPhi=(TH2D*)(hSourcePrimary->Project3D("yz"));
  for (int r=0;r<hSourcePrimary->GetYaxis()->GetNbins();r++){
    for (int z=0;z<hSourcePrimary->GetZaxis()->GetNbins();z++){
      double avePrim=hPhi->GetBinContent(z+1,r+1)/(nPbins*1.0);      
    for (int p=0;p<hSourcePrimary->GetXaxis()->GetNbins();p++){
	hOutputPrimary->SetBinContent(p+1,r+1,z+1,avePrim);
      }
    }
  }


  if(0){
  //compare to the originals:
  TCanvas *c =new TCanvas("c","c",1400,600);
  c->Divide(5,2);
  c->cd(1);
  hSourcePrimary->SetLineColor(kBlack);
  hSourcePrimary->Project3D("z2")->DrawCopy("P");
  hOutputPrimary->SetLineColor(kRed);
  hOutputPrimary->Project3D("z")->DrawCopy("same");
  c->cd(2);
  hOutputPrimary->Project3D("xy")->DrawCopy("colz");
  c->cd(3);
  hSourcePrimary->Project3D("xy2")->DrawCopy("colz");
  c->cd(4);
  hOutputPrimary->Project3D("yz3")->DrawCopy("colz");
  c->cd(5);
  hSourcePrimary->Project3D("yz4")->DrawCopy("colz");
  c->cd(6);
  hSourceIBF->SetLineColor(kBlack);
  hSourceIBF->Project3D("z4")->DrawCopy("P");
  hOutputIBF->SetLineColor(kRed);
  hOutputIBF->Project3D("z3")->DrawCopy("same");

  c->cd(7);
  hOutputIBF->Project3D("xy5")->DrawCopy("colz");
  c->cd(8);
  hSourceIBF->Project3D("xy6")->DrawCopy("colz");
  c->cd(9);
  hOutputIBF->Project3D("yz7")->DrawCopy("colz");
  c->cd(10);
  hSourceIBF->Project3D("yz8")->DrawCopy("colz");
  }
  
  //reassemble, scale, and save.
  hOutputPrimary->Scale(scale);
  hOutputPrimary->GetSumw2()->Set(0);
  hOutputIBF->Scale(scale);
  hOutputIBF->GetSumw2()->Set(0);
  TH3D *hOutputCharge=new TH3D(*hOutputPrimary);
  hOutputCharge->SetName("h_Charge_evt_0");
  hOutputCharge->SetTitle("Number of SC (ions)");
  hOutputCharge->Add(hOutputIBF);
  hOutputCharge->GetSumw2()->Set(0);
  hOutputCharge->Write();
  hOutputIBF->Write();
  hOutputPrimary->Write();
  output->Close();

  return;
}
