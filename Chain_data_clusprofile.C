//
{
  int sector = 1;
  TCanvas *c = new TCanvas("c","",900,600);
  //  c->Divide(2,1);
  c->SetLogy(0);
  c->SetLogx(0);
  c->Draw();
  c->Update();

  TH1F * clusterSumAdcVsLayerLow = new TH1F("clusterSumAdcVsLayerLow","clusterSumAdcVsLayerLow",113,-56.5,56.5);
  clusterSumAdcVsLayerLow->SetMarkerColor(kBlue);
  clusterSumAdcVsLayerLow->SetMarkerStyle(20);
  clusterSumAdcVsLayerLow->GetYaxis()->CenterTitle();
  clusterSumAdcVsLayerLow->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  clusterSumAdcVsLayerLow->SetXTitle("Layer");
  clusterSumAdcVsLayerLow->SetYTitle("Cluster Charge");
  clusterSumAdcVsLayerLow-> Draw();

  TH1F * clusterMaxAdcVsLayerLow = new TH1F("clusterMaxAdcVsLayerLow","clusterMaxAdcVsLayerLow",113,-56.5,56.5);
  clusterMaxAdcVsLayerLow->SetMarkerColor(kBlue);
  clusterMaxAdcVsLayerLow->SetMarkerStyle(20);
  clusterMaxAdcVsLayerLow->GetYaxis()->CenterTitle();
  clusterMaxAdcVsLayerLow->GetXaxis()->CenterTitle();
  clusterMaxAdcVsLayerLow->SetXTitle("Layer");
  clusterMaxAdcVsLayerLow->SetYTitle("Cluster Charge");
  clusterMaxAdcVsLayerLow-> Draw();



  TH1F *nEventsLow = new TH1F("nEventsLow","nEventsLow",3,-0.5,2.5);
  nEventsLow->SetMarkerColor(kBlue);
  nEventsLow->SetMarkerStyle(20);
  nEventsLow->GetYaxis()->CenterTitle();
  nEventsLow->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  nEventsLow->SetXTitle("Layer");
  nEventsLow->SetYTitle("Cluster Charge");
  nEventsLow-> Draw();

 TH1F *clusterSumAdcVsLayerMid = new TH1F("clusterSumAdcVsLayerMid","clusterSumAdcVsLayerMid",113,-56.5,56.5);
  clusterSumAdcVsLayerMid->SetMarkerColor(kBlack);
  clusterSumAdcVsLayerMid->SetMarkerStyle(20);
  clusterSumAdcVsLayerMid->GetYaxis()->CenterTitle();
  clusterSumAdcVsLayerMid->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  clusterSumAdcVsLayerMid->SetXTitle("Layer");
  clusterSumAdcVsLayerMid->SetYTitle("Cluster Charge");
  clusterSumAdcVsLayerMid-> Draw();

  TH1F *clusterMaxAdcVsLayerMid = new TH1F("clusterMaxAdcVsLayerMid","clusterMaxAdcVsLayerMid",113,-56.5,56.5);
  clusterMaxAdcVsLayerMid->SetMarkerColor(kBlack);
  clusterMaxAdcVsLayerMid->SetMarkerStyle(20);
  clusterMaxAdcVsLayerMid->GetYaxis()->CenterTitle();
  clusterMaxAdcVsLayerMid->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  clusterMaxAdcVsLayerMid->SetXTitle("Layer");
  clusterMaxAdcVsLayerMid->SetYTitle("Cluster Charge");
  clusterMaxAdcVsLayerMid-> Draw();

  TH1F *nEventsMid = new TH1F("nEventsMid","nEventsMid",3,-0.5,2.5);
  nEventsMid->SetMarkerColor(kBlack);
  nEventsMid->SetMarkerStyle(20);
  nEventsMid->GetYaxis()->CenterTitle();
  nEventsMid->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  nEventsMid->SetXTitle("Layer");
  nEventsMid->SetYTitle("Cluster Charge");
  nEventsMid-> Draw();

  TH1F *clusterSumAdcVsLayerHigh = new TH1F("clusterSumAdcVsLayerHigh","clusterSumAdcVsLayerHigh",113,-56.5,56.5);
  clusterSumAdcVsLayerHigh->SetMarkerColor(kRed);
  clusterSumAdcVsLayerHigh->SetMarkerStyle(20);
  clusterSumAdcVsLayerHigh->GetYaxis()->CenterTitle();
  clusterSumAdcVsLayerHigh->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  clusterSumAdcVsLayerHigh->SetXTitle("Layer");
  clusterSumAdcVsLayerHigh->SetYTitle("Cluster Charge");
  clusterSumAdcVsLayerHigh-> Draw();

   TH1F *clusterMaxAdcVsLayerHigh = new TH1F("clusterMaxAdcVsLayerHigh","clusterMaxAdcVsLayerHigh",113,-56.5,56.5);
  clusterMaxAdcVsLayerHigh->SetMarkerColor(kRed);
  clusterMaxAdcVsLayerHigh->SetMarkerStyle(20);
  clusterMaxAdcVsLayerHigh->GetYaxis()->CenterTitle();
  clusterMaxAdcVsLayerHigh->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  clusterMaxAdcVsLayerHigh->SetXTitle("Layer");
  clusterMaxAdcVsLayerHigh->SetYTitle("Cluster Charge");
  clusterMaxAdcVsLayerHigh-> Draw();

  TH1F *nEventsHigh = new TH1F("nEventsHigh","nEventsHigh",3,-0.5,2.5);
  nEventsHigh->SetMarkerColor(kRed);
  nEventsHigh->SetMarkerStyle(20);
  nEventsHigh->GetYaxis()->CenterTitle();
  nEventsHigh->GetXaxis()->CenterTitle();
  //phires0->SetYTitle("");
  //phires0->SetXTitle("");
  nEventsHigh->SetXTitle("Layer");
  nEventsHigh->SetYTitle("Cluster Charge");
  nEventsHigh-> Draw();


  TChain *clustrlo = new TChain("clustertree","clustrlo");
  TChain *evtrlo = new TChain("eventtree","evtrlo");
  TChain *clustrmid = new TChain("clustertree","clustrmid");
  TChain *evtrmid = new TChain("eventtree","evtrmid");
  TChain *clustrhi = new TChain("clustertree","clustrhi");
  TChain *evtrhi = new TChain("eventtree","evtrhi");

  /* low 18091 */
  clustrlo->Add("/sphenix/tg/tg01/hf/bogui/data/Prod_Tree_52321_*.root");
  evtrlo->Add("/sphenix/tg/tg01/hf/bogui/data/Prod_Tree_52321_*.root");
  
  /* mid 30512 */
  clustrmid->Add("/sphenix/tg/tg01/hf/bogui/data/Prod_Tree_52577_*.root");
  evtrmid->Add("/sphenix/tg/tg01/hf/bogui/data/Prod_Tree_52577_*.root");
  
  /* high 66179 */
  clustrhi->Add("/sphenix/tg/tg01/hf/bogui/data/Prod_Tree_52856_*.root");
  evtrhi->Add("/sphenix/tg/tg01/hf/bogui/data/Prod_Tree_52856_*.root");
  /**/

  clustrlo->Draw("2*(side-0.5)*layer>>clusterMaxAdcVsLayerLow","maxadc*(layer>=7&&layer<55)","");
  clustrlo->Draw("2*(side-0.5)*layer>>clusterSumAdcVsLayerLow","adc*(layer>=7&&layer<55)","");
  evtrlo->Draw("1>>nEventsLow","run>0");

  clusterMaxAdcVsLayerLow->Scale(1.0/nEventsLow->GetEntries());
  clusterSumAdcVsLayerLow->Scale(1.0/nEventsLow->GetEntries());
  clusterSumAdcVsLayerLow->Draw();

  clustrmid->Draw("2*(side-0.5)*layer>>clusterMaxAdcVsLayerMid","maxadc*(layer>=7&&layer<55)","");
  clustrmid->Draw("2*(side-0.5)*layer>>clusterSumAdcVsLayerMid","adc*(layer>=7&&layer<55)","");
  evtrmid->Draw("1>>nEventsMid","run>0");

  clusterMaxAdcVsLayerMid->Scale(1.0/nEventsMid->GetEntries());
  clusterSumAdcVsLayerMid->Scale(1.0/nEventsMid->GetEntries());
  clusterSumAdcVsLayerMid->Draw();

  clustrhi->Draw("2*(side-0.5)*layer>>clusterMaxAdcVsLayerHigh","maxadc*(layer>=7&&layer<55)","");
  clustrhi->Draw("2*(side-0.5)*layer>>clusterSumAdcVsLayerHigh","adc*(layer>=7&&layer<55)","");
  evtrhi->Draw("1>>nEventsHigh","run>0");

  clusterMaxAdcVsLayerHigh->Scale(1.0/nEventsHigh->GetEntries());
  clusterSumAdcVsLayerHigh->Scale(1.0/nEventsHigh->GetEntries());
  clusterSumAdcVsLayerHigh->Draw();

  clusterSumAdcVsLayerHigh->SetMaximum(400000);
  clusterSumAdcVsLayerHigh->Draw();
  clusterSumAdcVsLayerMid->Draw("same");
  clusterSumAdcVsLayerLow->Draw("same");

  TFile *f = new TFile("data_charge_profiles.root","RECREATE");
  clusterSumAdcVsLayerHigh->Write();
  clusterSumAdcVsLayerMid->Write();
  clusterSumAdcVsLayerLow->Write();
  clusterMaxAdcVsLayerHigh->Write();
  clusterMaxAdcVsLayerMid->Write();
  clusterMaxAdcVsLayerLow->Write();
  f->Close();

  /*
  phires0->Draw();
  c->Update();
  c->Print("phires_side0.gif");

  phires1->Draw();
  c->Update();
  c->Print("phires_side1.gif");
  
  zres0->Draw();
  c->Update();
  c->Print("zres_side0.gif");

  zres1->Draw();
  c->Update();
  c->Print("zres_side1.gif");
  */

  TLegend *t3=new TLegend(0.35,0.83,0.75,0.63);
  t3->AddEntry(clusterSumAdcVsLayerHigh,"All clusters, maxadc weighted","");
  t3->AddEntry(clusterSumAdcVsLayerHigh,"run 52856 <nclus> = 66179","P");
  t3->AddEntry(clusterMaxAdcVsLayerMid,"run 52577 <nclus> = 30512","P");
  t3->AddEntry(nhitslaylo,"run 52321 <nclus> = 18091","P");
  //  t3->AddEntry(hcent_4o12,"May Modular","P");
  t3->SetFillColor(0);
  t3->SetBorderSize(0);
  t3->SetFillStyle(0);
  t3->SetTextFont(63);
  t3->SetTextSize(20);
  t3->Draw();
  
  c->Update();
  c->Print("data_charge_profile_max.gif");
}
