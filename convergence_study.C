void convergence_study(){

  //for each ttree:
  //get the 'out' and the 'back' and compute the difference in the r and the phi direction
  //find the center and rms of that distribution
  //add those to an array of those numbers
  //plot those numbers vs reduction.


  int max_reduction=20+1;
  int red[max_reduction];
  float deltar[max_reduction];
  float deltarphi[max_reduction];
  float sigmar[max_reduction];
  float sigmarphi[max_reduction];

  TFile *f;
  TTree *t;
  TH1F *hR;
  TH1F *hRphi;
  for (int i=0;i<max_reduction;i++){
    f=TFile::Open(Form("pre-hybrid_fixed_reduction_%d.ttree.root",i));
    t=(TTree*)f->Get("pTree");
    t->Draw("(back1.Perp()-orig.Perp())>>histR");
    t->Draw("(back1.Phi()-orig.Phi())*orig.Perp()>>histRphi");
    hR=(TH1F*)gDirectory->Get("histR");
    hRphi=(TH1F*)gDirectory->Get("histRphi");

    deltar[i]=hR->GetMean();
    sigmar[i]=hR->GetRMS();
    deltarphi[i]=hRphi->GetMean();
    sigmarphi[i]=hRphi->GetRMS();
    red[i]=i;
    f->Close();
  }
  
  int i=0;
  TGraphErrors *rPhiDiff[3];
  TMultiGraph *rphi=new TMultiGraph();
  rphi->SetTitle("r and r*phi offset for varying reco grid sizes;r(um);rphi(um)");
  rPhiDiff[i]=new TGraphErrors(max_reduction,deltar,deltarphi,sigmar,sigmarphi);
  rPhiDiff[i]->SetTitle("rphi diff;size;t(ms)");
  rPhiDiff[i]->SetMarkerColor(kBlack);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraphErrors(1,deltar,deltarphi,sigmar,sigmarphi);
  rPhiDiff[i]->SetTitle(Form("rphi diff (native)^3 grid;size;t(ms)"));
  rPhiDiff[i]->SetMarkerColor(kRed);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraphErrors(1,deltar+max_reduction-1,deltarphi+max_reduction-1,sigmar+max_reduction-1,sigmarphi+max_reduction-1);
  rPhiDiff[i]->SetTitle(Form("rphi diff (native-%d)^3 grid;size;t(ms)",red[max_reduction-1]));
  rPhiDiff[i]->SetMarkerColor(kGreen);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rphi->Draw("AC*");
  
  return;
}
