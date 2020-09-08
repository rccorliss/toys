
void generate_average_th3(const char * inputpattern="./evgeny/*.root", const char *outputfilename="average.out.root"){

  TFile *outfile=TFile::Open(outputfilename,"RECREATE");
  TH3D *outhist;

  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  printf("found: %s\n",((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetUrl());//Title());//Print();

  TFile *infile;
  bool isFirstHist=true;
  int totalHists=0;

  for (int i=0;i<filelist->GetNFiles();i++){
   //for each file, find all histograms in that file.
    infile=TFile::Open(((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetUrl(),"READ");//gross.
    TList *keys=infile->GetListOfKeys();
    keys->Print();
    int nKeys=infile->GetNkeys();

    for (int j=0;j<nKeys;j++){
      TObject *tobj=infile->Get(keys->At(j)->GetName());
      //if this isn't a 3d histogram, skip it:
      bool isHist=tobj->InheritsFrom("TH3");
      if (!isHist) continue;
      
      if (isFirstHist){
	outfile->cd();
	outhist=new TH3D(*(TH3D*)tobj);
	isFirstHist=false;
	totalHists=1;
      } else {
	//todo:  add a check to make sure they have the same dimensions...
	outhist->Add((TH3D*)tobj);
	totalHists++;
      }
      printf("obj %d: getname: %s  inherits from TH3D:%d \n",j,tobj->GetName(),tobj->InheritsFrom("TH3"));
    }
      infile->Close();
  }
  outhist->Scale(1.0/(totalHists*1.0));

  //scale from total charge to charge density as a quick patch for evgeny's old files:
  TH3D *hDensity=new TH3D(*outhist);
  hDensity->SetName("hAverageChargeDensity");
  double hphimin=outhist->GetXaxis()->GetXmin();
  double hphimax=outhist->GetXaxis()->GetXmax();
  int hphibins=outhist->GetXaxis()->GetNbins();
  double hphistep=(hphimax-hphimin)/(1.0*hphibins);

  double hzmin=outhist->GetZaxis()->GetXmin();
  double hzmax=outhist->GetZaxis()->GetXmax();
  int hzbins=outhist->GetZaxis()->GetNbins();
  double hzstep=(hzmax-hzmin)/(1.0*hzbins);

  double hrmin=outhist->GetYaxis()->GetXmin();
  double hrmax=outhist->GetYaxis()->GetXmax();
  int hrbins=outhist->GetYaxis()->GetNbins();
  double hrstep=(hrmax-hrmin)/(1.0*hphibins);
  for (int i=0;i<hphibins;i++){
    for (int j=0;j<hrbins;j++){
      double hr=hrmin+j*hrstep;
      for (int k=0;k<hzbins;k++){
	double vol=hzstep*hphistep*(hr+0.5*hrstep)*hrstep;
	hDensity->SetBinContent(i+1,j+1,k+1,1.0/vol*outhist->GetBinContent(outhist->GetBin(i+1,j+1,k+1)));
      }
    }
  }

  outfile->cd();
  outhist->Write();
  hDensity->Write();
  outfile->Close();
    
  //for each histogram
  //if this is the first histogram, copy it to our output hist
  //if this is not the first histogram, check if it has the same bounds
  //if it does, add it to our output histogram

  //once we're done, divide by the number of histograms (scale by 1/n)

  //save our histogram.
  //
  return;
  


  /*

char *outputname[]={"cmflash/flash_002000_008000","evgeny/mapFile_bX734587_bias0","Smooth.50kHz","Single.50kHz"};
  char *scfilename[]={"outputFile_15kHz_G4Hits_sHijing_0-12fm_002000_008000_bX734587_bias10.root", "evgeny/mapFile_bX734587_bias0.root","Smooth.50kHz.root","BeamXingNBeams.root"};
  char *schistname[]={"h_Charge_0",
		      "h_Charge_1",
		      "h_Charge_2",
		      "h_Charge_3",
		      "h_Charge_4",
		      "h_Charge_5",
		      "h_Charge_6",
		      "h_Charge_7",
		      "h_Charge_8",
		      "h_Charge_9"};
  */
}
