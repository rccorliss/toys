
void subtract_average_th3(const char * inputpattern="./evgeny/*.root", const char *outputdir="./output/", const char *averagefilename="average.out.root"){

  const int nHists=1; //3 differential, 3 integral.
  
  TFile *avefile=TFile::Open(averagefilename,"READ");
  //actually, I have a stack of histograms...  are they in order?
  // int nAveKeys=infile->GetNkeys();

  TH3D *avehist[nHists];
  for (int i=0;i<nHists;i++){
    avehist[i]=(TH3D*)(avefile->Get(avefile->GetListOfKeys()->At(i)->GetName())); //get the ith hist from the average file.
  }

    
  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  printf("found: %s\n",((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetFile());//Title());//Print();

  //return;
  TFile *infile, *outfile;
  TH3D *inhist, *outhist;
  TString filename, basename, outname;
  bool isFirstHist=true;
  int totalHists=0;

  for (int i=0;i<filelist->GetNFiles();i++){
   //for each file, find all histograms in that file.
    filename=((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetUrl();
    //get our basename by tokenizing:
    TString tok;
    Ssiz_t from = 0, lastslash=0;
    while (filename.Tokenize(tok, from, "/")) {
      basename=tok;
      lastslash=from;
    }
    outname=outputdir;
    outname=outname+"/fluct_"+basename;
    
    infile=TFile::Open(filename,"READ");
    outfile=TFile::Open(outname,"RECREATE");
    printf("file %d: basename: %s\n", i,basename.Data());
    TList *keys=infile->GetListOfKeys();
    //keys->Print();
    int nKeys=infile->GetNkeys();
    if (nKeys<nHists) continue;//assert (1==2); //file doesn't have enough histograms in it.  Something's wrong.
    for (int j=0;j<nHists && j<nKeys;j++){
      TObject *tobj=infile->Get(keys->At(j)->GetName());
      printf("  obj %d: getname: %s  inherits from TH3D:%d , matching to %s\n",j,tobj->GetName(),tobj->InheritsFrom("TH3"),avehist[j]->GetName());
      outfile->cd();
      outhist=new TH3D(*(TH3D*)tobj);
      outhist->Add((TH3D*)tobj,avehist[j],1,-1);
      outhist->Write();

    }
    outfile->Close();
    infile->Close();
  }
  return;
  

}
