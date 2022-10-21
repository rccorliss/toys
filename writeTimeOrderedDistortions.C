//basic framework to read and compare events in a time ordered distortion file
//these files contain a TTree with xingnum, and three distortion map branches
#include "TTree.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TCanvas.h" //this prevents a lazy binding issue and/or is a magic spell.
#include "TH3F.h"
#include "TH2F.h"
#include "TH1F.h"

void writeTimeOrderedDistortions(bool subtractFirst=false, char *filename="/sphenix/user/rcorliss/distortion_maps/2022.07/TimeOrderedDistortions.root", char *inputpattern="/sphenix/user/rcorliss/distortion_maps/2022.07/*.distortion_map.hist.root"){

  TFile *treefile=TFile::Open(filename,"RECREATE");
  //TTree *tree=(TTree*)(file->Get("TimeDists"));

  TH3F *basehist[6];
  TH3F *temphist[6];
  int xingnum=0;
  std::string branchname[]={"hIntDistortionP_negz",
			  "hIntDistortionR_negz",
			  "hIntDistortionZ_negz",
			  "hIntDistortionP_posz",
			  "hIntDistortionR_posz",
			  "hIntDistortionZ_posz"};
  std::string histname[]={"hIntDistortionP_negz",
			  "hIntDistortionR_negz",
			  "hIntDistortionZ_negz",
			  "hIntDistortionP_posz",
			  "hIntDistortionR_posz",
			  "hIntDistortionZ_posz"};
  TTree *tree=new TTree("TimeDists", "TimeDists");
  tree->Branch("xingnum",&xingnum);
  for (int i=0;i<6;i++){
    temphist[i]=new TH3F(Form("temphist%d",i),Form("temphist%d",i),10,0,10,20,0,20,30,0,30);
    basehist[i]=new TH3F(Form("basehist%d",i),Form("basehist%d",i),10,0,10,20,0,20,30,0,30);
    tree->Branch(branchname[i].c_str(),&(temphist[i]));
  }
  printf("histograms built and branched.\n");

  //find all files that match the input string
  TFileCollection *filelist=new TFileCollection();
  filelist->Add(inputpattern);
  filelist->Print();
  printf("found %d files like: %s\n",filelist->GetNFiles(),((TFileInfo*)(filelist->GetList()->At(0)))->GetCurrentUrl()->GetFile());//Title());//Print();

  //return;
  TFile *infile;
  bool fileIsValid=true;
  bool isFirst=true;
  int nMaps=0;
  

  for (int i=0;i<filelist->GetNFiles() && i<100;i++){
    //for each file, find all histograms in that file.
    infile=TFile::Open(((TFileInfo*)(filelist->GetList()->At(i)))->GetCurrentUrl()->GetUrl(),"READ");//gross.
    fileIsValid=true;
    printf("=====> Trying File %d\n",i);
    if (!infile->IsOpen()) continue; //file didn't open right.  move on to the next one.
    TList *keys=infile->GetListOfKeys();
    for (int j=0;j<6;j++){
      temphist[j]=infile->Get<TH3F>(histname[j].c_str());
      if (!temphist[j]){
	fileIsValid=false; //histogram doesn't exist.  don't bother loading the other hists.
	break;
      }
    }
    if (!fileIsValid) {
      infile->Close();
	continue; //didn't get all our hists.  move on to the next file.
    }
    printf("=====> File %d is valid\n",i);
    xingnum=i;//temporary fix to paste something in there.
    if(subtractFirst){
      if (isFirst){
	for (int j=0;j<6;j++){
	  treefile->cd();
	  basehist[j]=(TH3F*)(temphist[j]->Clone());
	  //temphist[j]->Copy(*(basehist[j]));
	}
	isFirst=false;
      }
      for (int j=0;j<6;j++){
	printf("ptr j=%d:  b:%p\tt:%p\n",j,basehist[j],temphist[j]);
	temphist[j]->Add(basehist[j],-1);
      }
    }
    
    tree->Fill();
    nMaps++;
    infile->Close();
  }
  printf("nMaps=%d\n",nMaps);
      
  treefile->cd();
  tree->Write();
  treefile->Close();
  return;
}
