//basic framework to read and compare events in a time ordered distortion file
//these files contain a TTree with xingnum, and three distortion map branches
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"

void writeTimeOrderedDistortions(char *filename="/sphenix/user/rcorliss/distortion_maps/2022.07/TimeOrderedDistortions.root", char *inputpattern="/sphenix/user/rcorliss/distortion_maps/2022.07/*.distortion_map.hist.root"){

  TFile *treefile=TFile::Open(filename,"RECREATE");
  //TTree *tree=(TTree*)(file->Get("TimeDists"));

  TH3F *basehist[6];
  TH3F *temphist[6];
  int xingnum=0;
  std::string histname[]={"IntDistortionP_negz",
			  "IntDistortionR_negz",
			  "IntDistortionZ_negz",
			  "IntDistortionP_posz",
			  "IntDistortionR_posz",
			  "IntDistortionZ_posz"};

  TTree *tree=new TTree("TimeDists", "TimeDists");
  tree->Branch("xingnum",&xingnum);
  for (int i=0;i<6;i++){
    temphist[i]=new TH3F(Form("temphist%d",i),Form("temphist%d",i),10,0,10,20,0,20,30,0,30);
    tree->Branch(histname[i].c_str(),&(temphist[i]));
  }
  printf("histograms built and branched.\n");

  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  printf("found %d files like: %s\n",filelist->GetNFiles(),((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetFile());//Title());//Print();

  //return;
  TFile *infile;
  bool isFirstHist=true;
  int nMaps=0;
  

  for (int i=0;i<filelist->GetNFiles();i++){
    //for each file, find all histograms in that file.
    infile=TFile::Open(((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetUrl(),"READ");//gross.
    if (!infile->IsOpen()) continue; //file didn't open right.  move on to the next one.
    TList *keys=infile->GetListOfKeys();
    for (int i=0;i<6;i++){
      temphist[i]=infile->Get<TH3F*>(histname[i].c_str());
      if (!temphist[i]){
	fileIsValid=false; //histogram doesn't exist.  don't bother loading the other hists.
	break;
      }
    }
    if (!fileIsValid) continue; //didn't get all our hists.  move on to the next file.
    xingnum=i;//temporary fix to paste something in there.

    tree->Fill();
    nMaps++;
    infile->Close();
  }
      
  treefile->cd();
  tree->Write();
  treefile->Close();
  return;

  
