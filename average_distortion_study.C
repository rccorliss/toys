void average_distortion_study(){

  //for each ttree:
  //get the 'out' and the 'back' and compute the difference in the r and the phi direction
  //find the center and rms of that distribution
  //add those to an array of those numbers
  //plot those numbers vs reduction.

  int first_sample=1;
  int last_sample=99;
  int n_samples=20;
  float red[n_samples];
  float rpos[n_samples];
  float sigrpos[n_samples];
  float sigmared[n_samples];
  float deltar[n_samples];
  float deltarphi[n_samples];
  float sigmar[n_samples];
  float sigmarphi[n_samples];
  float deltarF[n_samples];
  float deltarphiF[n_samples];

  TFile *f;
  f=TFile::Open(Form("last_macro_output.ttree.root"));
  TTree *t;
  t=(TTree*)f->Get("pTree");

  //automatically get our axes:
  t->Draw("orig.Perp()>>raxis");
  t->Draw("orig.Phi()>>paxis");
  float rmin=20*1e4;//((TH1F*)gDirectory->Get("raxis"))->GetXaxis()->GetXmin();
  float rmax=78*1e4;//((TH1F*)gDirectory->Get("raxis"))->GetXaxis()->GetXmax();
  float rstep=(rmax-rmin)/n_samples;
  float pmin=((TH1F*)gDirectory->Get("paxis"))->GetXaxis()->GetXmin();
  float pmax=((TH1F*)gDirectory->Get("paxis"))->GetXaxis()->GetXmax();
  float pstep=(pmax-pmin)/n_samples;

  printf("Step sizes: rs=%f, ps=%f\n",rstep,pstep);
  TH1F *hR;
  TH1F *hRphi;
  TH1F *hRfrac;
  TH1F *hRphiFrac;
  for (int i=0;i<n_samples;i++){
    float r0=rmin+rstep*i;
    float r1=r0+rstep;
    float p0=pmin+pstep*i;
    float p1=p0+pstep;
    printf("loading step %d, %f<r<%f\n",i,r0,r1);
    assert(f->IsOpen());
    t->Draw("(out1.Perp()-orig.Perp())>>histR(100)",Form("orig.Perp()>%f && orig.Perp()<%f",r0,r1));
    t->Draw("(out1.Phi()-orig.Phi())*orig.Perp()>>histRphi(100)",Form("(out1.Phi()-orig.Phi()<3.14)&&(orig.Perp()>%f && orig.Perp()<%f)",r0,r1));
    hR=(TH1F*)gDirectory->Get("histR");
    hRphi=(TH1F*)gDirectory->Get("histRphi");
    rpos[i]=(r0+rstep/2)/1e4;
    sigrpos[i]=(rstep/2)/1e4;
    deltar[i]=hR->GetMean();
    sigmar[i]=hR->GetRMS();
    deltarphi[i]=hRphi->GetMean();
    sigmarphi[i]=hRphi->GetRMS();
  }
    f->Close();



  
 TGraphErrors *getemp;
 TMultiGraph *mgtemp;
 mgtemp=new TMultiGraph();
  mgtemp->SetTitle("r distortion from static E and B fieldmaps in absence of spacecharge;r position (cm);distortion(um)");
  getemp=new TGraphErrors(n_samples,rpos,deltar,sigrpos,sigmar);
  getemp->SetTitle("r diff (um);fraction of nominal grid;(um)");
  getemp->SetMarkerColor(kRed);
  getemp->SetMarkerStyle(kStar);
  mgtemp->Add(getemp);
  mgtemp->Draw("A*");

  return;
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
