int maxpercanvas=5;

void plotSpaceChargeCheck(char * inputfile){
  TFile *infile=TFile::Open(inputfile,"READ");
  TList *keys=infile->GetListOfKeys();
  keys->Print();
  TCanvas *c=new TCanvas("c","charge profiles",1200,800);

  int nKeys=infile->GetNkeys();
  if (nKeys<maxpercanvas){
    c->Divide(2,nKeys);
  } else {
        c->Divide(2,maxpercanvas);
  }

  int drawn=0;
  int nc=0;
  for (int j=0;j<nKeys;j++){
    TObject *tobj=infile->Get(keys->At(j)->GetName());
    //if this isn't a 3d histogram, skip it:
    bool isHist=tobj->InheritsFrom("TH3");
    if (drawn>=maxpercanvas){
      c=new TCanvas(Form("c%d",nc),"charge profiles",1200,800);
      nc++;
      c->Divide(2,maxpercanvas);
      drawn=0;
    }
    c->cd(2*drawn+1);
    ((TH3*)(tobj))->ProjectionZ()->Draw();
    c->cd(2*drawn+2);
    ((TH3*)(tobj))->Project3D("xy")->Draw("colz");
    drawn++;
  }
  return;
}
