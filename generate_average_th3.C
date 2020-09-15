
void generate_average_th3(const char * inputpattern="./evgeny/*.root", const char *outputfilename="average.out.root"){

  TFile *outfile=TFile::Open(outputfilename,"RECREATE");
  TH3D *outhist;

  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  printf("found: %s\n",((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetFile());//Title());//Print();

  //return;
  TFile *infile;
  bool isFirstHist=true;
  int totalHists=0;

  for (int i=0;i<filelist->GetNFiles();i++){
   //for each file, find all histograms in that file.
    infile=TFile::Open(((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetUrl(),"READ");//gross.
    TList *keys=infile->GetListOfKeys();
    //keys->Print();
    int nKeys=infile->GetNkeys();

    //for (int j=0;j<nKeys;j++){
    for (int j=2;j<3;j++){//just grab the third key -- h_avg_Charge_30;1 etc.
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
      printf("new charge=%E\ttotal charge=%E C\n",1.6e-19*((TH3D*)tobj)->Integral(),1.6e-19 * outhist->Integral());
    }
      infile->Close();
      if (totalHists>4) break;
  }
  outhist->Scale(1.0/(totalHists*1.0));
  printf("Done.  %d hists added. average charge=%E C\n",totalHists,1.6e-19 * outhist->Integral());


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


  //now we have two histograms:
  //hDensity -- 3d density profile from as many events as we found
  //outhit -- 3d charge-per-bin profile of same.
  //from Takao, we expect a certain functional form for the r dependence:
  TF1 * takaoRadial=new TF1("takaoradialscaling","-log(tan(atan2(x,1.05)/2))*[0]/x",0.05,1.0);
  TH1D *radial=hDensity->ProjectionY();
  radial->Scale(1.0/124.0*1.0/360.0);
  radial->Fit(takaoRadial);
  //and the z scaling is straightforward:
  TF1 * longitudinal=new TF1("zscaling","[0]-[1]*x",0,105.5);
  TH1D *longi=hDensity->ProjectionZ();
  longi->Scale(1.0/159.0*1.0/360.0);
  longi->Fit(longitudinal);

  /*
  //for the frames, scaling, we project out z from the output:
  TH1F *framemaskraw=hAverageChargeDensity->Project3D("xy");
  float threshold=framemask->GetMaximum()/2.0;
  //then we'll take everything /above/ half the the average and map that to 1.0, everything below the average and map that to zero to approximate the frame shadows.


  for (int i=0;i<hphibins;i++){
    for (int j=0;j<hrbins;j++){
      double hr=hrmin+j*hrstep;
      float rfactor.... blah;// rcc here.
      for (int k=0;k<hzbins;k++){
	double vol=hzstep*hphistep*(hr+0.5*hrstep)*hrstep;
	hSmoothAverage->SetBinContent(i+1,j+1,k+1,outhist->GetBinContent(outhist->GetBin(i+1,j+1,k+1)));
      }
    }
  }
  
*/

  
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
  

}
