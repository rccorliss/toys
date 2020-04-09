//proof of principle that we can recover a chain of unknown distortions in a simple model

void cm_reco_test(){
  TF1 *dDistribution=new TF1("dDistribution","gaus(0)",-1,1);
  dDistribution->SetParameters(1,0,0.3);//normalization,mean,width

  const int nCells=20;//how many distortion cells are there linearly in z in the model
  const int nRefreshes=10;//how many times do we completely cycle a new distortion through the model region
  const int nSteps=nCells*nRefreshes;//total number of steps that need to be generated for the full time series

  TH1F *hDist=new TH1F("hDist","true distortion distribution",nRefreshes,-1,1);
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
  xf[0]=distSum;
  for (int i=0;i<nPart;i++){
    distSum+=(-distort[i]+distort[i+nCells]);//remove the oldest distortion from the sum and add the new one
    xf[i+1]=distSum;
  }

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
    hMeasRel[offset]=new TH1F(Form("hMeasRel%d",offset),Form("Relative Distortion offset=%d;distortion",offset),nRefreshes*10,-10,10);
    measuredRelativeDistort[offset][0]=0;//the dummy value for the distortion from the 0th refresh.
    for (int j=1;j<nRefreshes;j++){
      measuredRelativeDistort[offset][j]=measuredRelativeDistort[offset][j-1]+xf[(nCells)*(j-1)+offset+1]-xf[(nCells)*(j-1)+offset];
      hMeasRel[offset]->Fill(measuredRelativeDistort[offset][j]);
    }
  }

  //now we know that each of these ought to have the same distribution, since they're drawn from the same sample.
  //we can also generate a proxy for the mean of that sample by looking at the total distortion measured by each test particle
  //and dividing that by the number of cells:
  TH1F *hAveDistortion=new TH1F("hAveDistortion","Average distortion for each refresh;distortion mag.",nRefreshes,-1,1);
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

  TH1F *hDistMatch=new TH1F("hDistMatch","reco-true distortion per step;reco-true",nRefreshes,-0.3,0.3);
  for (int i=0;i<nSteps;i++){
    int offset=i%nCells;
    int refresh=i/nCells;
    hDistMatch->Fill(recoOffset[offset]+measuredRelativeDistort[offset][refresh]-distort[i]);
  }
  hDistMatch->Draw();
  return;
}
