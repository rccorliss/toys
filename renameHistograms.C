/*
This is a little toy to rename histograms and only save a subset from one file to another
 */

//std::string infileName="";
//std::string outfileName="";


vector<std::pair<std::string,std::string>> fileName;
vector<std::pair<std::string,std::string>> histName;

int main(){
  renameHistograms();
  return 0;
}

void renameHistograms(){
  histName.push_back(std::make_pair("hIntDistortionNegR","hIntDistortionR_negz"));
  histName.push_back(std::make_pair("hIntDistortionNegP","hIntDistortionP_negz"));
  histName.push_back(std::make_pair("hIntDistortionNegZ","hIntDistortionZ_negz"));
  histName.push_back(std::make_pair("hIntDistortionPosR","hIntDistortionR_posz"));
  histName.push_back(std::make_pair("hIntDistortionPosP","hIntDistortionP_posz"));
  histName.push_back(std::make_pair("hIntDistortionPosZ","hIntDistortionZ_posz"));

  fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/workfest2021/full_distortion.workfest2021.distortion_map.hist.root",
				    "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/actual_event.distortion_map.hist.root"));
  fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/workfest2021/empty_distortion.workfest2021.distortion_map.hist.root",
				    "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/static_only.distortion_map.hist.root"));
  fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/workfest2021/average_distortion.workfest2021.distortion_map.hist.root",
				    "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/average_event.distortion_map.hist.root"));
  /*
    fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/workfest2021/",
				    "/star/u/rcorliss/sphenix/workfest2021/"));
  */

  
  TFile *infile;=TFile::Open(infileName,"READ");
  TFile *outfile;=TFile::Open(outfileName,"RECREATE");

  TH1* hin,hout;

  for (int f=0;f<fileName.length(),f++){

    infile=TFile::Open(fileName[f].first(),"READ");
    outfile=TFile::Open(fileName[f].second(),"RECREATE");
    
    for (int i=0;i<histName.length();i++){
      hin=(TH1*)infile->Get(histName[i].first());
      hout=(TH1*)hin->Clone(histName[i].second());
    }
    
    outfile->cd();  
    for (int i=0;i<histName.length();i++){
      hout->Write();
    }
    outfile->Close();
    infile->Close();
  }
    return;
}
