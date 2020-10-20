
void plotSpaceChargeCheck(char * inputfile){
  TFile *infile=TFile::Open(inputfile,"READ");
  TList *keys=infile->GetListOfKeys();
  keys->Print();
  TCanvas *c=new TCanvas("c","charge profiles",800,800);

  int nKeys=infile->GetNkeys();
  c->Divide(1,nKeys);


  for (int j=0;j<nKeys;j++){
    TObject *tobj=infile->Get(keys->At(j)->GetName());
    //if this isn't a 3d histogram, skip it:
    bool isHist=tobj->InheritsFrom("TH3");
    c->cd(j+1);
    ((TH3*)(tobj))->ProjectionZ()->Draw();
  }
  return;
}
