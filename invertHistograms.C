
vector<std::pair<std::string,std::string>> fileName;
vector<std::string> histName;

const int resample_factor=3; //how many samples, evenly spaced but centered in the bin boundaries so that they're never on the edge unless n=inf, to resample.  may have aliasing issues if we're not careful.

void Resample(std::vector<TH3*> hin, std::vector<TH3*> hout);
void CheckClosure(std::vector<TH3*> hdistort, std::vector<TH3*> hcorrect);

void invertHistograms(){
  histName.push_back("hIntDistortionR_negz");
  histName.push_back("hIntDistortionP_negz");
  histName.push_back("hIntDistortionZ_negz");
  histName.push_back("hIntDistortionR_posz");
  histName.push_back("hIntDistortionP_posz");
  histName.push_back("hIntDistortionZ_posz");

  fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/trackingStudySampleNov2021/actual_event.distortion_map.hist.root",
				    "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/actual_event_invert.distortion_map.hist.root"));
  fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/trackingStudySampleNov2021/static_only.distortion_map.hist.root",
				    "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/static_only_invert.distortion_map.hist.root"));
  fileName.push_back(std::make_pair("/star/u/rcorliss/sphenix/trackingStudySampleNov2021/average_event.distortion_map.hist.root",
				    "/star/u/rcorliss/sphenix/trackingStudySampleNov2021/average_event_invert.distortion_map.hist.root"));

  
  TFile *infile;
  TFile *outfile;

  std::vector<TH3*> hin;
  std::vector<TH3*> hout;

  for (int f=0;f<fileName.size();f++){

    infile=TFile::Open(fileName[f].first.data(),"READ");
    outfile=TFile::Open(fileName[f].second.data(),"RECREATE");
    outfile->cd();  

    for (int i=0;i<3;i++){
      hin.push_back((TH3*)infile->Get(histName[i].data()));
      hout.push_back((TH3*)hin[i]->Clone(histName[i].data()));
    }
    Resample(hout,hin);
    for (int i=0;i<3;i++){
      hout[i]->Write();
    }
    CheckClosure(hin,hout);

    hin.clear();
    hout.clear();
    
    for (int i=3;i<6;i++){
      hin.push_back((TH3*)infile->Get(histName[i].data()));
      hout.push_back((TH3*)hin[i]->Clone(histName[i].data()));
    }
    Resample(hout,hin);
    for (int i=0;i<3;i++){
      hout[i]->Write();
    }

    CheckClosure(hin,hout);
    
    outfile->Close();
    infile->Close();
  }
    return;
}

void Resample(std::vector<TH3*> hin, std::vector<TH3*> hout){
  TH3* hhits=(TH3*)hin[0]->Clone("hhits"); //number of elements in each output bin, for normalization purposes.

  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hhits->GetXaxis();
  ax[1]=hhits->GetYaxis();
  ax[2]=hhits->GetZaxis();

  int nbins[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();
  }

  float low[3],high[3], step[3], sample_pos[3];
  float distortion[3], distorted_pos[3];
  int a=0;
  for (int i=0;i<nbins[0];i++){
    a=0;
    low[a]=ax[a]->GetBinLowEdge(i+1);
    high[a]=ax[a]->GetBinUpEdge(i+1);
    step[a]=(high[a]-low[a])/(1.*resample_factor);
    for (int isub=0;isub<resample_factor;isub++){
      a=0;
      sample_pos[a]=step[a]*(isub+0.5);
      for (int j=0;j<nbins[1];j++){
	a=1;
	low[a]=ax[a]->GetBinLowEdge(j+1);
	high[a]=ax[a]->GetBinUpEdge(j+1);
	step[a]=(high[a]-low[a])/(1.*resample_factor);
	for (int jsub=0;jsub<resample_factor;jsub++){
	  a=1;
	  sample_pos[a]=step[a]*(jsub+0.5);
	  for (int k=0;k<nbins[2];k++){
	    a=2;
	    low[a]=ax[a]->GetBinLowEdge(k+1);
	    high[a]=ax[a]->GetBinUpEdge(k+1);
	    step[a]=(high[a]-low[a])/(1.*resample_factor);
	    for (int ksub=0;ksub<resample_factor;ksub++){
	      a=2;
	      sample_pos[a]=step[a]*(ksub+0.5);

	      //get the distorted position
	      for (int m=0;m<3;m++){
		distortion[m]=hin[m]->GetBinContent(hin[m]->FindBin(sample_pos[0],sample_pos[1],sample_pos[2]));
		distorted_pos[m]=sample_pos[m]+distortion[m];
	      }

	      //histogram the distortion in the distorted position.
	      hhits->Fill(distorted_pos[0],distorted_pos[1],distorted_pos[2],1);
	      for (int m=0;m<3;m++){
		hout[m]->Fill(distorted_pos[0],distorted_pos[1],distorted_pos[2],distortion[m]);
	      }
	    }
	  }
	}
      }
    }
  }


  //normalize based on how many hits we put into each.  yell if we have less than 1/2 the sampling.
  for (int i=0;i<nbins[0];i++){
    a=0;
    distorted_pos[a]=ax[a]->GetBinCenter(i+1);
    for (int j=0;j<nbins[1];j++){
      a=1;
      distorted_pos[a]=ax[a]->GetBinCenter(j+1);
      for (int k=0;k<nbins[2];k++){
	a=2;
	distorted_pos[a]=ax[a]->GetBinCenter(k+1);
	//histogram the distortion in the distorted position.
	int global_bin=hin[0]->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]);
	float global_content=hhits->GetBinContent(global_bin);
	if (global_content<1.0){
	  printf("(%2.2f,%2.2f,%2.2f) has %1.2f entries\n",distorted_pos[0],distorted_pos[1],distorted_pos[2],global_content);
	}
	for (int m=0;m<3;m++){
	  hout[m]->SetBinContent(global_bin,hout[m]->GetBinContent(global_bin)/global_content);
	}     
      }
    }
  }
  return;
}



void CheckClosure(std::vector<TH3*> hdistort, std::vector<TH3*> hcorrect){
  TH3* hclosure[3];
  TH1F* hresidual[3];
  static int ncalls=0;
  ncalls++;
  for(int i=0;i<3;i++){
    hclosure[i]=(TH3*)hdistort[0]->Clone(Form("hclosure%d_%d",i,ncalls));
    hresidual[i]=new TH1F(Form("hresidual%d_%d",i,ncalls),Form("residual in axis %d",i),200,-0.1,0.1);
  }
  TH3* hhits=(TH3*)hdistort[0]->Clone("hhits"); //number of elements in each output bin, for normalization purposes.

  TAxis *ax[3]={nullptr,nullptr,nullptr};
  ax[0]=hdistort[0]->GetXaxis();
  ax[1]=hdistort[0]->GetYaxis();
  ax[2]=hdistort[0]->GetZaxis();

  int nbins[3];
  for (int i=0;i<3;i++){
    nbins[i]=ax[i]->GetNbins();
  }

  float low[3],high[3], step[3], sample_pos[3];
  float distortion[3], distorted_pos[3];
  float correction[3], corrected_pos[3];
  int a=0;
  for (int i=0;i<nbins[0];i++){
    a=0;
    low[a]=ax[a]->GetBinLowEdge(i+1);
    high[a]=ax[a]->GetBinUpEdge(i+1);
    step[a]=(high[a]-low[a])/(1.*resample_factor);
    for (int isub=0;isub<resample_factor;isub++){
      a=0;
      sample_pos[a]=step[a]*(isub+0.5);
      for (int j=0;j<nbins[1];j++){
	a=1;
	low[a]=ax[a]->GetBinLowEdge(j+1);
	high[a]=ax[a]->GetBinUpEdge(j+1);
	step[a]=(high[a]-low[a])/(1.*resample_factor);
	for (int jsub=0;jsub<resample_factor;jsub++){
	  a=1;
	  sample_pos[a]=step[a]*(jsub+0.5);
	  for (int k=0;k<nbins[2];k++){
	    a=2;
	    low[a]=ax[a]->GetBinLowEdge(k+1);
	    high[a]=ax[a]->GetBinUpEdge(k+1);
	    step[a]=(high[a]-low[a])/(1.*resample_factor);
	    for (int ksub=0;ksub<resample_factor;ksub++){
	      a=2;
	      sample_pos[a]=step[a]*(ksub+0.5);

	      //get the distorted position
	      for (int m=0;m<3;m++){
		distortion[m]=hdistort[m]->GetBinContent(hdistort[m]->FindBin(sample_pos[0],sample_pos[1],sample_pos[2]));
		distorted_pos[m]=sample_pos[m]+distortion[m];
	      }
	      //get the corrected position
	      for (int m=0;m<3;m++){
		correction[m]=hcorrect[m]->GetBinContent(hcorrect[m]->FindBin(distorted_pos[0],distorted_pos[1],distorted_pos[2]));
		corrected_pos[m]=distorted_pos[m]-correction[m];
	      }

	      //get the position residual
	      for (int m=0;m<3;m++){
		hclosure[m]->Fill(sample_pos[0],sample_pos[1],sample_pos[2],distortion[m]-correction[m]);
		hresidual[m]->Fill(distortion[m]-correction[m]);
	      }
	    }
	  }
	}
      }
    }
  }
  //TFile f=TFile::Open("ClosurePlots.hist.root");
  for (int i=0;i<3;i++){
    hclosure[i]->Write();
    hresidual[i]->Write();
  }
  //f->Close();
  return;
}

int main(){
  invertHistograms();
  return 0;
}
