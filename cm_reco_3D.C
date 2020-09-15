//proof of principle that we can recover a chain of unknown distortions in a simple model

vector<TVector3> ReconstructDistortionSeries(vector<TVector3> recoOffset, int nSets);

void cm_reco_3D(){

  const int nXingPerTpcRefill; //how many bunch clocks it takes to drift completely across the TPC
  const int nXingPerLaserFlash; //how many bunch clocks between consecutive laser flashes.
  const int nSets=nXingPerTpcRefill/nXingPerLaserFlash;//how many distinct subsets we have.
  const int nPads=10;//put the right number here!


  //vector of reco positions for each pad.
  vector<vector<TVector3>> recoPosition;//outer vector is over pad id, inner vector is flash #
  for (int i=0;i<nPads;i++){
    vector<TVector3> padRelativePosition;
    relativePosition.push_back(padRelativePosition);
  }


    //get clusters from events
  //sort clusters by most likely pad
  //compute delta between true pad xyz and reco cluster xy0
  PopulateRecoPositions(&relativePosition); //something opaque that generates the relative positions for hits associated with each pad.

  //now we have sets of clusters with fixed spacing (eventually, we should try to generalize this)

  //remove the static distortion from electric/magnetic fields

  //remove the average distortion from spacecharge

  //per pad, reconstruct the differentials as best we can:

  vector<vector<TVector3>> relativeDistortion;//outer vector is over pad id, inner vector is flash #
  for (int i=0;i<nPads;i++){
    relativeDistortion.push_back(ReconstructDistortionSeries(relativePosition[i]), nSets);
  }

  

    
   TH1F *hMeasRel[nCells];
  for (int offset=0;offset<nCells;offset++){
    hMeasRel[offset]=new TH1F(Form("hMeasRel%d",offset),Form("Relative Distortion offset=%d;distortion",offset),nHistBins*2,-10*fluctScale,10*fluctScale);
    measuredRelativeDistort[offset][0]=0;//the dummy value for the distortion from the 0th refresh.
    for (int j=1;j<nRefreshes;j++){
      measuredRelativeDistort[offset][j]=measuredRelativeDistort[offset][j-1]+xf[(nCells)*(j-1)+offset+1]-xf[(nCells)*(j-1)+offset];
      hMeasRel[offset]->Fill(measuredRelativeDistort[offset][j]);
    }
  }

  //now we know that each of these ought to have the same distribution, since they're drawn from the same sample.
  //we can also generate a proxy for the mean of that sample by looking at the total distortion measured by each test particle
  //and dividing that by the number of cells:
  TH1F *hAveDistortion=new TH1F("hAveDistortion","Average distortion for each refresh;distortion mag.",nHistBins,-3*fluctScale,3*fluctScale);
  for (int i=0;i<nRefreshes;i++){
    hAveDistortion->Fill(xf[i*nCells]/nCells);
  }
  float overallMean=hAveDistortion->GetMean();

  //so we know/suspect the average of each of our measured distortions should be the same as the average of the mean
  //assuming we have enough samples.  This allows us to key in the one missing term -- the initial offset.
  float recoOffset[nCells];
  for (int i=0;i<nCells;i++){
    recoOffset[i]=overallMean-hMeasRel[i]->GetMean();
  }
  //now we compare our guess against reality.

  TH1F *hDistMatch=new TH1F("hDistMatch","reco-true distortion per step;reco-true",nHistBins*4,-3*fluctScale,3*fluctScale);
  TH2F *hDistMatch2D=new TH2F("hDistMatch2D","reco vs true distortion per step;reco;true",nHistBins*4,-3*fluctScale,3*fluctScale,nHistBins*4,-3*fluctScale,3*fluctScale);

  for (int i=0;i<nSteps;i++){
    int offset=i%nCells;
    int refresh=i/nCells;
    hDistMatch->Fill(recoOffset[offset]+measuredRelativeDistort[offset][refresh]-distort[i]);
    hDistMatch2D->Fill(recoOffset[offset]+measuredRelativeDistort[offset][refresh],distort[i]);
  }
  

  

  //all distances are measured in cm

  const float fluctWidth=0.3;//cm fluctuation per cm of drift -- width of that distribution
  const float zLength=105.5;//length of the tpc, used to calculate the fluctuation per cell.
  const int nCells=20;//how many distortion cells are there linearly in z in the model
  const float fluctScale=fluctWidth*zLength/nCells;
  const int nRefreshes=10;//how many times do we completely cycle a new distortion through the model region
  const int nSteps=nCells*nRefreshes;//total number of steps that need to be generated for the full time series
  const float xErrRaw=1e-10;//130e-4;//130um.  the uncertainty in a single x measurement, used to sample from a gaussian with this width as a noise factor added to the true xf.
  const int nLaserShots=1;//number of times we fire the laser.  The effective error will be reduced by sqrt(this)
  const float xErr=xErrRaw/sqrt(nLaserShots);
  
  const int nHistBins=(2*nRefreshes)>100?100:2*nRefreshes;//to try to auto-fit the size of the histograms so they don't over/under segment.

  //the distribution of distortion values in cm:
    TF1 *dDistribution=new TF1("dDistribution","gaus(0)",-5*fluctScale,5*fluctScale);
    dDistribution->SetTitle("Cell Distortion Underlying Distribution");
  dDistribution->SetParameters(1,0,fluctScale);//normalization,mean,width

  //the distribution of measurement smear in cm:
  TF1 *xSmear=new TF1("xSmear","gaus(0)",-1,1);
  xSmear->SetParameters(1,0,xErr);
  

  TH1F *hDist=new TH1F("hDist","true distortion distribution",nHistBins,-3*fluctScale,3*fluctScale);
  float distort[nSteps];//the true time series of distortions
  for (int i=0;i<nSteps;i++){
    distort[i]=dDistribution->GetRandom();
    hDist->Fill(distort[i]);
  }

  //in this simple model, we shoot particles through at high enough speed that the distortions are stationary.
  //we assume in addition that the distortions are uniform for each z slice - no x dependence.  The total
  //distortion is thus the sum of the distortions in each slice, for a given time.  Since we know where we throw,
  //we need only accumulate the distortion, and assume wlog that x0=0 for all of them.
  //we also notably assume that the distortions do not themselves evolve as they migrate.  They merely translate.

  const int nPart=nSteps-nCells;//throw one particle every time a new distortion cycles in.
  float xf[nPart];
  float distSum=0;
  //load the starting distortion into the distortion sum:
  for (int i=0;i<nCells;i++){
    distSum+=distort[i];
  }
  xf[0]=distSum+xSmear->GetRandom();
  for (int i=0;i<nPart;i++){
    distSum+=(-distort[i]+distort[i+nCells]);//remove the oldest distortion from the sum and add the new one
    xf[i+1]=distSum+xSmear->GetRandom();
}

  float xi,yi,zi;
  float xf,yf,zf;

  /* for Sara:
//making a tree using the get-distortions stuff written above
  float partR,partP,partZ;
  float xiSara,xfSara;
  TTree *sTree=new TTree("saraTree","made-up blah.");
  sTree->Branch("xi",&xiSara);
  sTree->Branch("xf",&xfSara);

  for (int i=0;i<inTree->GetEntries();i++){
xiSara=math to get xi;
xfSara=math to get xf;
	sTree->Fill();
}


  //filename comes from henry.
  char *treename="whateverHenryCallsTheTree";
  TFile *input=TFile::Open(fname);
  TTree *inTree=(TTree*)input->Get(treename);
  inTree->SetBranchAddress("xbefore",&xi);
  inTree->SetBranchAddress("xafter",&xf);
  //repeat for all other variable names
  for (int i=0;i<inTree->GetEntries();i++){
    inTree->GetEntry(i);
    printf("i=%d xi=%f\n",i, xi, xf);
    //process this laser flash, store it the way it was being stored before?
  }
    input->Close();
  */
  
  //rcc note:  rethink the ncells size.  I tried to catch all my errors, but it was late.
  
  //now we do the math.  We know there are nCells segments, so:
  // xf[i+1]-xf[i]=distort[i+nCells]-distort[i], which relates distortions in groups spaced by nCells+1:
  //          distort[i+nCells]=xf[i+1]-xf[i]+distort[i]; -- we know everything but distort[i-1]
  //and onward:
  //xf[i+nCells+1]-xf[i+nCells]=distort[i+2*nCells]-distort[i+nCells]
  //                           =distort[i+2*nCells]-(xf[i+1]-xf[i]+distort[i])
  //      distort[i+2*nCells]=xf[i+nCells+1]-xf[i+nCells]+xf[i+1]-xf[i]+distort[i], etc.
  float measuredRelativeDistort[nCells][nRefreshes];//the first refresh is the part we're relative to, dummied out to 0.
  TH1F *hMeasRel[nCells];
  for (int offset=0;offset<nCells;offset++){
    hMeasRel[offset]=new TH1F(Form("hMeasRel%d",offset),Form("Relative Distortion offset=%d;distortion",offset),nHistBins*2,-10*fluctScale,10*fluctScale);
    measuredRelativeDistort[offset][0]=0;//the dummy value for the distortion from the 0th refresh.
    for (int j=1;j<nRefreshes;j++){
      measuredRelativeDistort[offset][j]=measuredRelativeDistort[offset][j-1]+xf[(nCells)*(j-1)+offset+1]-xf[(nCells)*(j-1)+offset];
      hMeasRel[offset]->Fill(measuredRelativeDistort[offset][j]);
    }
  }

  //now we know that each of these ought to have the same distribution, since they're drawn from the same sample.
  //we can also generate a proxy for the mean of that sample by looking at the total distortion measured by each test particle
  //and dividing that by the number of cells:
  TH1F *hAveDistortion=new TH1F("hAveDistortion","Average distortion for each refresh;distortion mag.",nHistBins,-3*fluctScale,3*fluctScale);
  for (int i=0;i<nRefreshes;i++){
    hAveDistortion->Fill(xf[i*nCells]/nCells);
  }
   float overallMean=hAveDistortion->GetMean();

  //so we know/suspect the average of each of our measured distortions should be the same as the average of the mean
  //assuming we have enough samples.  This allows us to key in the one missing term -- the initial offset.
  float recoOffset[nCells];
  for (int i=0;i<nCells;i++){
    recoOffset[i]=overallMean-hMeasRel[i]->GetMean();
  }
  //now we compare our guess against reality.

  TH1F *hDistMatch=new TH1F("hDistMatch","reco-true distortion per step;reco-true",nHistBins*4,-3*fluctScale,3*fluctScale);
  TH2F *hDistMatch2D=new TH2F("hDistMatch2D","reco vs true distortion per step;reco;true",nHistBins*4,-3*fluctScale,3*fluctScale,nHistBins*4,-3*fluctScale,3*fluctScale);

  for (int i=0;i<nSteps;i++){
    int offset=i%nCells;
    int refresh=i/nCells;
    hDistMatch->Fill(recoOffset[offset]+measuredRelativeDistort[offset][refresh]-distort[i]);
    hDistMatch2D->Fill(recoOffset[offset]+measuredRelativeDistort[offset][refresh],distort[i]);
  }

  TCanvas *c=new TCanvas("c","cm_reco_test.C",1000,600);
  c->Divide(3,2);
  c->cd(1);
  dDistribution->Draw();
  TLatex *tex=new TLatex(0.0,0.8,Form("#sigma =%f",dDistribution->GetParameter(2)));
  tex->Draw();
  c->cd(2);
  hDist->Draw();
  c->cd(3);
  hMeasRel[1]->Draw();
  c->cd(4);
  hAveDistortion->Draw();
  c->cd(5);
  hDistMatch->Draw();
  c->cd(6);
  tex=new TLatex(0.0,0.8,Form("N Laser Shots=%d  N Refreshes=%d  N Cells=%d",nLaserShots, nRefreshes, nCells));
  tex->Draw();
  tex=new TLatex(0.0,0.6,Form("#sigma_{distortion} =%1.3f cm",dDistribution->GetParameter(2)));
  tex->Draw();
   tex=new TLatex(0.0,0.5,Form("#sigma_{hit} =%1.3f #mum",xSmear->GetParameter(2)*1e4));
  tex->Draw();
  tex=new TLatex(0.0,0.4,Form("Mean of average distortion per refresh =%1.3f cm",hAveDistortion->GetMean()));
  tex->Draw();
  tex=new TLatex(0.0,0.3,Form("Mean of reco-true =%1.3f #mum",hDistMatch->GetMean()*1e4));
  tex->Draw();
 tex=new TLatex(0.0,0.2,Form("RMS of reco-true =%1.3f #mum",hDistMatch->GetRMS()*1e4));
  tex->Draw();
  c->SaveAs("cm_reco.pdf");
  return;
}

	
vector<TVector3> ReconstructDistortionSeries(vector<TVector3> recoOffset, int nSets){
  //the distortion is a static component, plus a time-varying component.  Assuming that the time-varying component is localized, we can extract it by comparing the delta of when it leaves -- associating it with the new term that just entered.  See the document for more details.

  int nCycles=recoOffset.length()/nSets;


  //determine the time-average integrated distortion.  Not valid unless nCycles>>1.
  TVector3 staticOffset;
  for (int i=0;i<recoOffset.length();i++){
    staticOffset+=recoOffset[i];
  }
  staticOffset*=1.0/(1.0*recoOffset.length());

  
  //divide the distortions into correlated sets, subtracting off the static component.
  TVector3 offset[nSets][nCycles];
  for (int i=0;i<nCycles*nSets;i++){
    offset[i%nSets][i/nSets]=recoOffset[i]-staticOffset;
  }

  //for each correlated set compute the deltas of the time series.
  TVector3 relative[nSets][nCycles];
  for (int i=0;i<nSets;i++){
    relative[i][0].SetXYZ(0,0,0);
    for (int j=1;j<nCycles;j++){
      relative[i][j]=offset[i][j]-offset[i][j-1];
    }
  }
  //each set of relative distortions ought to have the same distribution, As long as the distortions associated with charge don't grow or shrink in any direction as they flow through the tpc (this is not true, but may be approximately true once we subtract off the average integrated distortion behavior).  Since we subtracted off the average behavior, the deltas in each case ought to be centered on zero.

  //I need to think harder on this than I can at this time of night.  Draw pictures, write in the document, grow understanding.

  //but if we take enough total distortions, we can average those, then subtract that off, so that we're left with only the varying parts.  We can learn nothing about the z distribution of the static distortion, but we couldn't do that anyway.
      
 
  //now we do the math.  We know there are nCells segments, so:
  // xf[i+1]-xf[i]=distort[i+nCells]-distort[i], which relates distortions in groups spaced by nCells+1:
  //          distort[i+nCells]=xf[i+1]-xf[i]+distort[i]; -- we know everything but distort[i-1]
  //and onward:
  //xf[i+nCells+1]-xf[i+nCells]=distort[i+2*nCells]-distort[i+nCells]
  //                           =distort[i+2*nCells]-(xf[i+1]-xf[i]+distort[i])
  //      distort[i+2*nCells]=xf[i+nCells+1]-xf[i+nCells]+xf[i+1]-xf[i]+distort[i], etc.
  float measuredRelativeDistort[nCells][nRefreshes];//the first refresh is the part we're relative to, dummied out to 0.
  TH1F *hMeasRel[nCells];
  for (int offset=0;offset<nCells;offset++){
    hMeasRel[offset]=new TH1F(Form("hMeasRel%d",offset),Form("Relative Distortion offset=%d;distortion",offset),nHistBins*2,-10*fluctScale,10*fluctScale);
    measuredRelativeDistort[offset][0]=0;//the dummy value for the distortion from the 0th refresh.
    for (int j=1;j<nRefreshes;j++){
      measuredRelativeDistort[offset][j]=measuredRelativeDistort[offset][j-1]+xf[(nCells)*(j-1)+offset+1]-xf[(nCells)*(j-1)+offset];
      hMeasRel[offset]->Fill(measuredRelativeDistort[offset][j]);
    }
  }

  //now we know that each of these ought to have the same distribution, since they're drawn from the same sample.
  //we can also generate a proxy for the mean of that sample by looking at the total distortion measured by each test particle
  //and dividing that by the number of cells:
  TH1F *hAveDistortion=new TH1F("hAveDistortion","Average distortion for each refresh;distortion mag.",nHistBins,-3*fluctScale,3*fluctScale);
  for (int i=0;i<nRefreshes;i++){
    hAveDistortion->Fill(xf[i*nCells]/nCells);
  }
   float overallMean=hAveDistortion->GetMean();

  //so we know/suspect the average of each of our measured distortions should be the same as the average of the mean
  //assuming we have enough samples.  This allows us to key in the one missing term -- the initial offset.
  float recoOffset[nCells];
  for (int i=0;i<nCells;i++){
    recoOffset[i]=overallMean-hMeasRel[i]->GetMean();
  }
  vector<TVector3> output;

  return output;
}
