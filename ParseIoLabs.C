bool AverageDirectionPositive(vector<double>pos,int i){
  const int depth=6;
  int start=(i-depth);
  if (start<0) start=0;
  double beforelen=i-start;
  double before=0;
  int end=i+depth;
  if (end>pos.size()) end=pos.size();
  double afterlen=end-i;
  double after;

  for (int j=start;j<i;j++){
    before+=pos[j]/beforelen;
  }
  for (int j=i;j<end;j++){
    after+=pos[j]/afterlen;
  }

  if (after>before) return true;
  return false;
}


void ParseIoLabs(const char * title, int wheelNumber, int lightNumber, int surveyMode=0){
  const int date=20210506;
  TString baseName=Form("/Users/rcorliss/Documents/IOLab-WorkFiles/export/%d-",date);
  TString wheelFileName=baseName+wheelNumber+"_Wheel.csv";
  TString lightFileName=baseName+lightNumber+"_Light.csv";
  TTree *wheelTree=new TNtuple("wheelTree","wheelTree","index:frame:sample:time:rawX:calX:rawY:calY:rawZ:calZ");
  TTree *lightTree=new TNtuple("wheelTree","wheelTree","indexL:frameL:sampleL:timeL:rawL:calL");
  /*
  char wheelname[100]="/Users/rcorliss/Documents/IOLab-WorkFiles/export/20210503-110942_Wheel.csv";
    std::ifstream in;
     in.open(wheelname);
     if (!in.good()) {
       printf("file is broken somehow\n");
       return;
     }
     printf("file is fine?\n");
     return;
  */
  //note you have to manually remove the top line.
  
  
  wheelTree->ReadFile(wheelFileName,"",',');
  lightTree->ReadFile(lightFileName,"",',');
  wheelTree->AddFriend(lightTree);
  //TTree *lightTree=TTree::ReadFile(wheelFileName,"",",");
  //wheelTree->Draw("time-timeL","1","colz");
  TGraph *gt;
  TMultiGraph *mgt=new TMultiGraph();
  
  float pos,posBefore;
  bool goingUp,goingUpBefore;
  bool first=true;//
  int nEntries=wheelTree->Draw("calL:calX","1","goff");
  std::vector<double>posVec;
  std::vector<double>lightVec;
  for (int i=0;i<nEntries;i++){
    posVec.push_back(1000.*wheelTree->GetVal(1)[i]);
    lightVec.push_back(wheelTree->GetVal(0)[i]);
  }
  std::vector<int>setStart={0};
  std::vector<int>setLength;

  //TF1 *staticGaus=new TF1("sgaus","gaus(0)+[3]",-50,100);
  TF1 *staticGaus=new TF1("sgaus","[0]*exp(-0.5*(((x-[1])/[2])**2)**[4])+[3]",-50,100);
  const int nParams=5;
  std::vector<double>fitParams[nParams];
  double fitMin[nParams];
  double fitMax[nParams];

  staticGaus->SetParameters(0.05,0,50,0);
  staticGaus->SetParLimits(0,0.,8.);//no inverted gaussians
  staticGaus->SetParLimits(1,-60.,-30.);//center must be on the left somewhere
  staticGaus->SetParLimits(2,0.,80.);//width isn't too wide
  staticGaus->SetParLimits(3,-0.1,0.1);//ambient light isn't too bright.

  //adjust parameter limits based on the tuning for particular passes:
  if (wheelNumber==111314) {
    staticGaus->SetParLimits(0,0.,0.1);//no inverted gaussians
    staticGaus->SetParLimits(1,-20.,-20.);//center must be on the left somewhere
    staticGaus->SetParLimits(2,0.,80.);//width isn't too wide
    staticGaus->SetParLimits(3,-0.1,0.1);//ambient light isn't too bright.
    staticGaus->SetParLimits(4,0.75,1.25);//supergaussian isn't strong.
  }
 if (wheelNumber==111528) {
    staticGaus->SetParLimits(0,0.,0.1);//no inverted gaussians
    staticGaus->SetParLimits(1,0.,40.);//center must be on the left somewhere
    staticGaus->SetParLimits(2,0.,80.);//width isn't too wide
    staticGaus->SetParLimits(3,-0.1,0.1);//ambient light isn't too bright.
    staticGaus->SetParLimits(4,0.85,1.15);//supergaussian isn't strong.
  }
 if (wheelNumber==111743) {
    staticGaus->SetParLimits(0,0.,0.1);//no inverted gaussians
    staticGaus->SetParLimits(1,-20.,40.);//center must be on the left somewhere
    staticGaus->SetParLimits(2,0.,80.);//width isn't too wide
    staticGaus->SetParLimits(3,-0.1,0.1);//ambient light isn't too bright.
    staticGaus->SetParLimits(4,0.85,1.15);//supergaussian isn't strong.
  } 
  if (wheelNumber==122900) {
    staticGaus->SetParLimits(0,0.,0.1);//no inverted gaussians
    staticGaus->SetParLimits(1,-40.,-0.);//center must be on the left somewhere
    staticGaus->SetParLimits(2,0.,80.);//width isn't too wide
    staticGaus->SetParLimits(3,-0.1,0.1);//ambient light isn't too bright.
  }
  if (wheelNumber==123507) {
    staticGaus->SetParLimits(0,0.,0.1);//no inverted gaussians
    staticGaus->SetParLimits(1,-70.,-30.);//center must be on the left somewhere
    staticGaus->SetParLimits(2,0.,90.);//width isn't too wide
    staticGaus->SetParLimits(3,-0.1,0.1);//ambient light isn't too bright.
    staticGaus->SetParLimits(4,0.75,1.25);//supergaussian isn't strong.
  }

  if (date==20210506) {
    staticGaus->SetParLimits(0,0.,9);//no inverted gaussians
    staticGaus->SetParLimits(1,-50.,100.);//center must be in the middle
    staticGaus->SetParLimits(2,0.,800.);//width isn't too wide
    staticGaus->SetParLimits(3,0,0.055);//ambient light isn't negative, nor too bright.
    staticGaus->SetParLimits(4,0.45,3.5);//supergaussian can be strong.
  if (wheelNumber==123507) {
  }

  }  

  
  //break this up into its different passes, based on whether it seems to be going up or down on average
  TCanvas *cp;
  if (surveyMode==1){
    cp=new TCanvas("cp","cp",1400,800);
 cp->Divide(7,4);
  }

  bool keep[100];
  if (surveyMode==2){
    for (int i=0;i<100;i++){
      keep[i]=false;
    }
    //curate which ones to keep:
    if (wheelNumber==111314) {int keeps[]={2,3,6,12,14};for (int i=0;i<5;i++){keep[keeps[i]]=true;}}
    if (wheelNumber==111528) {for (int i=4;i<23;i++){keep[i]=true;}}
    if (wheelNumber==111743) {for (int i=4;i<26;i++){keep[i]=true;}}
    if (wheelNumber==123040) {for (int i=10;i<18;i++){keep[i]=true;}}
    if (wheelNumber==122900) {for (int i=2;i<12;i++){keep[i]=true;}}
    if (wheelNumber==123507) {int keeps[]={2,3,6,7,8,10,11,12,13};for (int i=0;i<9;i++){keep[keeps[i]]=true;}}
  }
    
  
  goingUpBefore=false;//(posVec[1]-posVec[0])>0;
  int goodSetNumber=0;
  for (int i=3;i<posVec.size()-2;i++){
    goingUp=AverageDirectionPositive(posVec,i);
    //    goingUp=(posVec[i-3]+posVec[i-2]+posVec[i-1])<(posVec[i]+posVec[i+1]+posVec[i+2]);
    //goingUp=(posVec[i]-posVec[i-1])>0;
    if (goingUp!=goingUpBefore){
      setLength.push_back(i-setStart.back());
      int setNumber=setLength.size()-1;
      printf("switch directions at i=%d.  set %lu has length %d\n",i,setLength.size()-1,setLength.back());
      gt=new TGraph(setLength.back(),&(posVec[setStart.back()]),&(lightVec[setStart.back()]));
      gt->SetMarkerColor(1+setNumber%8);
      gt->SetLineColor(1+setNumber%8);
      gt->SetTitle(Form("Pass %d (%s);position [mm];light intensity [arb]",setNumber,goingUpBefore?"leftward":"rightward"));
      setStart.push_back(i);
      goingUpBefore=goingUp;

      if(setLength.back()>100){//skip short sets)
      
	if (surveyMode==1){
	  goodSetNumber+=1;
	  gt->Fit(staticGaus);
	  
	  cp->cd(goodSetNumber);
	  gt->SetLineColor(kBlack);
	  gt->SetMarkerColor(kBlack);
	  gt->Draw("A*L");
	}
	//gt->SetLineColor(1+setNumber%8);
	//gt->SetMarkerColor(1+setNumber%8);
	
	if (surveyMode==0 || (surveyMode==2 && keep[setNumber])){
	  mgt->Add(gt);
	  gt->Fit(staticGaus);
	
	  printf("fit params: ");
	  for (int j=0;j<nParams;j++){
	    float fitPar=gt->GetFunction("sgaus")->GetParameter(j);
	    if (j==2 && fitPar<0) fitPar*=-1;//force width to be positive.
	    fitParams[j].push_back(fitPar);
	    printf("%d: %1.2E\t",j,fitPar);
	    //track minimum and maximum values of the parameters
	    if (fitParams[j].size()==1 || fitPar<fitMin[j]){
	      fitMin[j]=fitPar;
	    }
	    if (fitParams[j].size()==1 || fitPar>fitMax[j]){
	      fitMax[j]=fitPar;
	    }
	  }
	  printf("\n");
	}
      }
    }
  }
  if (surveyMode==1) return;

  //take each pass and fit it to a gaussian (done above)
  TH1F *hFitPar[nParams];
  TString fitParName[]={"coeff","mean","sigma","offset","extra exponent"};
  float fitMargin[nParams];//so we don't clip high and low values;
  for (int i=0;i<nParams;i++){
    printf ("Before forcing, %s has min %f and max %f\n",fitParName[i].Data(),fitMin[i],fitMax[i]);
    fitMargin[i]=0.5*(fitMax[i]-fitMin[i]);
    }
  //fitMin[0]=0;
  //fitMax[0]=8;
  //fitMin[1]=-20;
  //fitMax[1]=20;
  //fitMin[2]=0;
  //fitMax[2]=50;
  //fitMin[3]=-0.1;
  //fitMax[3]=0.3;
  
  for (int i=0;i<nParams;i++){
    printf ("%s has min %f and max %f\n",fitParName[i].Data(),fitMin[i],fitMax[i]);
    hFitPar[i]=new TH1F(Form("hFitPar%d",i),fitParName[i],40,fitMin[i]-fitMargin[i],fitMax[i]+fitMargin[i]);
    hFitPar[i]->FillN(fitParams[i].size(),&(fitParams[i][0]),NULL); //fill the histogram with a elements from the array at b, with weights c (or null)
  }

  
  TCanvas *c=new TCanvas("c","c",1400,300);
  gStyle->SetOptStat(111111);
  gStyle->SetStatW(0.4);
  gStyle->SetStatH(0.25);
   c->Divide(1+nParams,1);
  c->cd(1);
  mgt->SetTitle(Form("%s (%lu good passes);position [mm];light [arb]",title,fitParams[0].size()));
  mgt->Draw("A*L");

  if (0){
  float DiffuseBounds=200;
  float cardDist=327;
  int nPhotons=1e6;
  float background=hFitPar[3]->GetMean();
  float center=hFitPar[1]->GetMean();
  float area=sqrt(2*3.14)*hFitPar[0]->GetMean()*hFitPar[2]->GetMean(); //not quite right for supergaussians.
  
 TF1* fEdAngleDeg = new TF1("fEdAngle", "([0]+[1]*exp(-(x-[2])**2/(2*[3]**2)))", -DiffuseBounds, DiffuseBounds);
  fEdAngleDeg->SetParameters(0, 0.9483733002234809, 0, 16.1331011690544950);//test of gaussian fit to image of 40 deg

  TH1F *hEd=new TH1F("hEd","Expected distribution from edmonton fit;pos(mm)",100,-DiffuseBounds,DiffuseBounds);
  for (int i=0;i<nPhotons;i++){
    float angle=fEdAngleDeg->GetRandom();
    float position=center+TMath::Tan(angle*TMath::Pi()/180.)*cardDist;
    hEd->Fill(position,1./(1.*nPhotons));
  }
  //scale it:
  float edArea=hEd->Integral();
  hEd->Scale(area/edArea);
  for (int i=0;i<100;i++){
    hEd->Fill(hEd->GetXaxis()->GetBinCenter(i+1),background);
  }
  hEd->SetLineColor(kBlue);
  hEd->Draw("same,hist");
  }
  

  
  TLegend *leg=c->cd(1)->BuildLegend(0.4,0.6,0.4,0.6);
  //leg->SetMargin(0.05);
  //leg->SetNColumns(2);
  leg->SetEntrySeparation(0.05);
  leg->SetColumnSeparation(0.01);
  for(int i=0;i<nParams;i++){
    c->cd(2+i);
    hFitPar[i]->Draw();
  }
 

  printf("Note:  Some fit limits are hard-coded into this macro, and may need to be adjusted to get the gaussians to behave.  If they come out wrong, adjust the   staticGaus->SetParLimits() calls\n");      
  //chopper:  if position p[i] has: direction[i]=p[i]-p[i-1] : dir[i]!=dir[i-1] chop to a new pass.  
  //return;
}
