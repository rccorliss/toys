void model_shift_study(){

  //for each ttree:
  //get the 'out' and the 'back' and compute the difference in the r and the phi direction
  //find the center and rms of that distribution
  //add those to an array of those numbers
  //plot those numbers vs reduction.


  int ntests=12;
  int nominal=10;
  int mintest=-10; //relative to nominal.
  int step=2;
  int red[ntests];
  float deltar[ntests];
  float deltarphi[ntests];
  float sigmar[ntests];
  float sigmarphi[ntests];
  float offset[ntests];
  float sigmaoffset[ntests];

  TFile *f;
  TTree *t;
  TH1F *hR;
  TH1F *hRphi;
  for (int i=0;i<ntests;i++){
    int ioffset=nominal+i*step+mintest;
    offset[i]=ioffset-nominal;
    sigmaoffset[i]=step/2;
    if (ioffset==10) {
      f=TFile::Open(Form("full3d_10cm_direct.ttree.root"));
    } else{
    f=TFile::Open(Form("full3d_%dcm_10cm.ttree.root",nominal+i*step+mintest));
    }
    t=(TTree*)f->Get("pTree");
    t->Draw("(back1.Perp()-orig.Perp())>>histR");
    t->Draw("(back1.Phi()-orig.Phi())*orig.Perp()>>histRphi");
    hR=(TH1F*)gDirectory->Get("histR");
    hRphi=(TH1F*)gDirectory->Get("histRphi");

    //multiply by 5 because I'm going over 1/5 the range.
    deltar[i]=abs(5*hR->GetMean());
    sigmar[i]=5*hR->GetRMS();
    deltarphi[i]=abs(5*hRphi->GetMean());
    sigmarphi[i]=5*hRphi->GetRMS();
    red[i]=i;
    f->Close();
  }


  //   inom*step+mintest=0 ==>
  int inom=-mintest/step;
  int i=0;
  TGraphErrors *rDiff[4];
  TMultiGraph *rDiffM=new TMultiGraph();
  rDiffM->SetTitle("r out-and-back residual for varying charge model offsets;model offset (cm);delta r(um)");
  rDiff[i]=new TGraphErrors(ntests,offset,deltar,sigmaoffset,sigmar);
  rDiff[i]->SetTitle("rphi diff");
  rDiff[i]->SetMarkerColor(kBlack);
  rDiff[i]->SetMarkerStyle(kStar);
  rDiffM->Add(rDiff[i]);
  i++;
   rDiff[i]=new TGraphErrors(1,offset,deltar,sigmaoffset,sigmar);
   rDiff[i]->SetTitle(Form("%dcm",mintest));
  rDiff[i]->SetMarkerColor(kRed);
  rDiff[i]->SetMarkerStyle(kStar);
  rDiffM->Add(rDiff[i]);
  i++;
  rDiff[i]=new TGraphErrors(1,offset+inom,deltar+inom,sigmaoffset+inom,sigmar+inom);
  rDiff[i]->SetTitle(Form("nominal"));
  rDiff[i]->SetMarkerColor(kBlue);
  rDiff[i]->SetMarkerStyle(kStar);
  rDiffM->Add(rDiff[i]);
  i++;
  rDiff[i]=new TGraphErrors(1,offset+ntests-1,deltar+ntests-1,sigmaoffset+ntests-1,sigmar+ntests-1);
  rDiff[i]->SetTitle(Form("%dcm",mintest+(ntests-1)*step));
  rDiff[i]->SetMarkerColor(kGreen);
  rDiff[i]->SetMarkerStyle(kStar);
  rDiffM->Add(rDiff[i]);
  i++;
  rDiffM->Draw("AC*");

  new TCanvas();
  i=0;
  TGraphErrors *rPhiDiff[4];
  TMultiGraph *rPhiDiffM=new TMultiGraph();
  rPhiDiffM->SetTitle("r*phi out-and-back residual for varying charge model offsets;model offset (cm);r*delta phi(um)");
  rPhiDiff[i]=new TGraphErrors(ntests,offset,deltarphi,sigmaoffset,sigmarphi);
  rPhiDiff[i]->SetTitle("rphi diff");
  rPhiDiff[i]->SetMarkerColor(kBlack);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rPhiDiffM->Add(rPhiDiff[i]);
  i++;
   rPhiDiff[i]=new TGraphErrors(1,offset,deltarphi,sigmaoffset,sigmarphi);
   rPhiDiff[i]->SetTitle(Form("%dcm",mintest));
  rPhiDiff[i]->SetMarkerColor(kRed);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rPhiDiffM->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraphErrors(1,offset+inom,deltarphi+inom,sigmaoffset+inom,sigmarphi+inom);
  rPhiDiff[i]->SetTitle(Form("nominal"));
  rPhiDiff[i]->SetMarkerColor(kBlue);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rPhiDiffM->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraphErrors(1,offset+ntests-1,deltarphi+ntests-1,sigmaoffset+ntests-1,sigmarphi+ntests-1);
  rPhiDiff[i]->SetTitle(Form("%dcm",mintest+(ntests-1)*step));
  rPhiDiff[i]->SetMarkerColor(kGreen);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rPhiDiffM->Add(rPhiDiff[i]);
  i++;
  rPhiDiffM->Draw("AC*");
  
  return;
}
