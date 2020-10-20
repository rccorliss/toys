
void plotSpaceChargeCheck(char * inputfile){
  TFile *infile=TFile::Open(inputfile,"READ");
  TList *keys=infile->GetListOfKeys();
  keys->Print();
  TCanvas *c=new TCanvas("c","charge profiles",1200,800);

  int nKeys=infile->GetNkeys();
  c->Divide(2,nKeys);


  for (int j=0;j<nKeys;j++){
    TObject *tobj=infile->Get(keys->At(j)->GetName());
    //if this isn't a 3d histogram, skip it:
    bool isHist=tobj->InheritsFrom("TH3");
    c->cd(2*j+1);
    ((TH3*)(tobj))->ProjectionZ()->Draw();
    c->cd(2*j+2);
    ((TH3*)(tobj))->Project3D("xy")->Draw("colz");
  }
  return;
}
