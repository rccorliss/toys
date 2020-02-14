void convergence_study(){

  //for each ttree:
  //get the 'out' and the 'back' and compute the difference in the r and the phi direction
  //find the center and rms of that distribution
  //add those to an array of those numbers
  //plot those numbers vs reduction.

  int first_sample=1;
  int last_sample=99;
  int n_samples=last_sample-first_sample+1;
  float red[n_samples];
  float sigmared[n_samples];
  float deltar[n_samples];
  float deltarphi[n_samples];
  float sigmar[n_samples];
  float sigmarphi[n_samples];
  float deltarF[n_samples];
  float deltarphiF[n_samples];

  TFile *f;
  TTree *t;
  TH1F *hR;
  TH1F *hRphi;
  TH1F *hRfrac;
  TH1F *hRphiFrac;
  for (int i=0;i<n_samples;i++){
    f=TFile::Open(Form("output/analytic_fixed_fraction_f%d.ttree.root",i+first_sample));
    printf("loading file %d\n",i);
    assert(f->IsOpen());
    t=(TTree*)f->Get("pTree");
    t->Draw("(out1.Perp()-outa.Perp())>>histR");
    t->Draw("(out1.Phi()-outa.Phi())*orig.Perp()>>histRphi","(out1.Phi()-outa.Phi()<3.14)");
    hR=(TH1F*)gDirectory->Get("histR");
    hRphi=(TH1F*)gDirectory->Get("histRphi");
   t->Draw("(out1.Perp()-outa.Perp())/(outa.Perp()-orig.Perp())>>histRfrac");
    t->Draw("(out1.Phi()-outa.Phi())/(outa.Phi()-orig.Phi())>>histRphiFrac");
    hR=(TH1F*)gDirectory->Get("histR");
    hRphi=(TH1F*)gDirectory->Get("histRphi");
    hRfrac=(TH1F*)gDirectory->Get("histRfrac");
    hRphiFrac=(TH1F*)gDirectory->Get("histRphiFrac");
    deltar[i]=hR->GetMean();
    deltarF[i]=hRfrac->GetMean();
    sigmar[i]=hR->GetRMS();
    deltarphi[i]=hRphi->GetMean();
    deltarphiF[i]=hRphiFrac->GetMean();
    sigmarphi[i]=hRphi->GetRMS();
    red[i]=i+first_sample;
    sigmared[i]=0;
    f->Close();
  }



  
 TGraphErrors *getemp;
 TMultiGraph *mgtemp;
 mgtemp=new TMultiGraph();
  mgtemp->SetTitle("r and r*phi residual for varying reco grid sizes;fraction of nominal grid;residual(um)");
  getemp=new TGraphErrors(n_samples,red,deltar,sigmared,sigmar);
  getemp->SetTitle("r diff (um);fraction of nominal grid;(um)");
  getemp->SetMarkerColor(kRed);
  getemp->SetMarkerStyle(kStar);
  mgtemp->Add(getemp);
  getemp=new TGraphErrors(n_samples,red,deltarphi,sigmared,sigmarphi);
  getemp->SetTitle("r*(phi diff) (um);fraction of nominal grid;(um)");
   getemp->SetMarkerColor(kBlue);
  getemp->SetMarkerStyle(kStar);
  mgtemp->Add(getemp);
  mgtemp->Draw("A*");

 TCanvas *c2 =  new TCanvas();
 c2->cd();
   mgtemp=new TMultiGraph();
  mgtemp->SetTitle("r/driftr and phi/driftphi residual for varying reco grid sizes;fraction of nominal grid;residual(um)");
  getemp=new TGraphErrors(n_samples,red,deltarF,sigmared,sigmar);
  getemp->SetTitle("r diff/r drift;fraction of nominal grid;(um)");
  getemp->SetMarkerColor(kRed);
  getemp->SetMarkerStyle(kStar);
  mgtemp->Add(getemp);
  getemp=new TGraphErrors(n_samples,red,deltarphiF,sigmared,sigmarphi);
  getemp->SetTitle("(phi diff/ phi drift;fraction of nominal grid;(um)");
   getemp->SetMarkerColor(kBlue);
  getemp->SetMarkerStyle(kStar);
  mgtemp->Add(getemp);
  mgtemp->Draw("AC*");
  return;

  int i=0;
  TGraphErrors *rPhiDiff[3];
  TMultiGraph *rphi=new TMultiGraph();
  rphi->SetTitle("r and r*phi offset for varying reco grid sizes;r(um);rphi(um)");
  rPhiDiff[i]=new TGraphErrors(n_samples,deltar,deltarphi,sigmar,sigmarphi);
  rPhiDiff[i]->SetTitle("rphi diff;size;t(ms)");
  rPhiDiff[i]->SetMarkerColor(kBlack);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraphErrors(1,deltar,deltarphi,sigmar,sigmarphi);
  rPhiDiff[i]->SetTitle(Form("rphi diff (full)^3 grid;size;t(ms)"));
  rPhiDiff[i]->SetMarkerColor(kRed);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rPhiDiff[i]=new TGraphErrors(1,deltar+n_samples-1,deltarphi+n_samples-1,sigmar+n_samples-1,sigmarphi+n_samples-1);
  rPhiDiff[i]->SetTitle(Form("rphi diff (full-%1.0f)^3 grid;size;t(ms)",red[n_samples-1]));
  rPhiDiff[i]->SetMarkerColor(kGreen);
  rPhiDiff[i]->SetMarkerStyle(kStar);
  rphi->Add(rPhiDiff[i]);
  i++;
  rphi->Draw("AC*");
  
  return;
}
