//run with no argument to get the default panel of stripe results
//run with a positive argument >2 to get details about a certain petal's stripes
//run with argument==2 to get the default panel of survey mark results
//run with a negative argument to get the details about that petal's survey marks (not yet implemented)
//see ParseOGP() to adjust file and path names.


#include <iostream>
#include <string>




const int MAXOGPSTEP=223;
const bool build_gnuplot=false;
const std::string directory="/Users/rcorliss/Desktop/Output.Files.3/";
string GetHeaderLine(int step){
  string ret;
  //if (step<7 || step>MAXOGPSTEP) ret="FAIL!  Fix your code.";
  //if (step>=7 && step<=MAXOGPSTEP){
    ret=Form("Step %d",step);
    //}
  
  return ret;
}

double GetDatum(string leader, string spacer, int element, fstream *file){
  //this was initially intended to look for tab spacers, but it looks like it's just ' ', not '\t'.
  //I hardcoded in the correct spacing element.
  /*
    the pattern I keyed to was 'find the first instance of 'mm', and know that the numbers begin after that
    with single spaces between them, so if I want the 2nd number, I need to find string that starts
    after the 2nd space after 'mm'.  I then parse that substring as a floating point number, which all are data are.
  */
  float ret=-9000;
  string tp;
  while(getline(*file, tp) ){ //read data from file object and put it into string.
    if (tp.find(leader)!=std::string::npos){
      int pos=tp.find("mm",0);
      //printf("mm at %d, ",pos);
      for (int i=0;i<element;i++){
	pos=tp.find(" ",pos+1);
	//printf("space at %d, ",pos);
      }
      int end=tp.find(" ",pos+1);
      string datum=tp.substr(pos+1, end-pos);
      //printf("datum at %d:%s\n",pos,datum.c_str());

      ret=atof(datum.c_str());
      return ret;
    }
  }
  return -999;
}

void FillStepData(fstream *file, int step, float *val, TVector2* pos ){
  //should check state!
  //switch from OGP step number to logical array number:
  int index=step-11;
  bool isSurvey=false;
  if (index<0){
    //then we're in the survey set.  get the survey index:
    index+=4;
    isSurvey=true;
  }
  TVector2 thisPos;

  
  if (isSurvey){
    //get the diameter, which comes before xy location
    val[index]=GetDatum("Diameter","\t",2, file);
  }
  //get location
  pos[index].SetX(GetDatum("X Location","\t",2, file));
  pos[index].SetY(GetDatum("Y Location","\t",2, file));
  
  if (!isSurvey){
    if (index==2){
      printf("skipping step 13 area because that item is not written out \n");
    }else{
    //get the area, which comes after xy location
    val[index]=GetDatum("Area","\t",2, file);
    }
  }

  return;
}


void AddStepData(fstream *file, int step, float *val, TVector2* pos ){
  //should check state!
  //switch from OGP step number to logical array number:
  int index=step-11;
  bool isSurvey=false;
  if (index<0){
    //then we're in the survey set.  get the survey index:
    index+=4;
    isSurvey=true;
  }
  TVector2 thisPos;

  
  if (isSurvey){
    printf("AddStepData adds another contour's information.  doesn't work for the circles\n");
    //get the diameter, which comes before xy location
    return;
  }
  //get location
  TVector2 newpos;
  newpos.SetX(GetDatum("X Location","\t",2, file));
  newpos.SetY(GetDatum("Y Location","\t",2, file));

  //get area
  float area=0;
  area=GetDatum("Area","\t",2, file);

  //construct the average position, and the sum area
  TVector2 avepos=pos[index]*val[index]+newpos*area;
  val[index]+=area;
  pos[index]=avepos/(val[index]);

  return;
}


void ParseOGP(){

  int petal, index;
  
  float surveyDiameter[4];
  TVector2 surveyPosition[4];

  float stripeArea[213];
  TVector2 stripePosition[213];

  TFile *petalFile=TFile::Open(Form("%sall_petals.hist.root",directory.c_str()),"RECREATE");
  TTree *pTree=new TTree("pTree","petal OGP summary tree");
  pTree->Branch("petal",&petal);
  for (int i=0;i<4;i++){
    pTree->Branch(Form("surveyDiameter%d",i),&(surveyDiameter[i]));
    pTree->Branch(Form("surveyCenter%d",i),&(surveyPosition[i]));
  }
  for (int i=0;i<213;i++){
    pTree->Branch(Form("stripeArea%d",i),&(stripeArea[i]));
    pTree->Branch(Form("stripeCenter%d",i),&(stripePosition[i]));
  }

  TTree *stripeTree=new TTree("stripeTree","petal OGP stripe tree");
  TVector2 tpos;
  float tarea;
  stripeTree->Branch("petal",&petal);
  stripeTree->Branch("i",&index);
  stripeTree->Branch("pos",&tpos);
  stripeTree->Branch("area",&tarea);


  TTree *surveyTree=new TTree("surveyTree","petal OGP survey mark tree");
  float tdiam;
  surveyTree->Branch("petal",&petal);
  surveyTree->Branch("i",&index);
  surveyTree->Branch("pos",&tpos);
  surveyTree->Branch("diam",&tdiam);

  TVector2 empty(0,0);
  
  
  fstream txtfile;
  for (int petalIndex=0;petalIndex<95;petalIndex++){
    if (petalIndex==24) continue;
    txtfile.open(Form("%spetal%d_edges.txt",directory.c_str(),petalIndex),ios::in);
    //txtfile.open("/Users/rcorliss/Desktop/petal55_edges.txt",ios::in); 
    for (int i=0;i<213;i++){
      stripeArea[i]=-900;
      stripePosition[i]=empty;
    }
    for (int i=0;i<213;i++){
      surveyDiameter[i]=-900;
      surveyPosition[i]=empty;
    }
    
     
    if (txtfile.is_open()){   //checking whether the file is open
      //if the file exists, we do the parsing and output.
      printf("opened %d\n", petalIndex);
      int readerstep=7;//first item we look for is step 7
      int maxadjust=0;
      //if (petalIndex==70) maxadjust=1;
      int adjust=0;//special adjustment to the text file contour numbering for certain petals.
      string tp;
      while(getline(txtfile, tp) && readerstep<=MAXOGPSTEP){ //read data from file object and put it into string.
	//printf("readerstep=%d\n",readerstep);
	//cout << tp << "\n"; //print the data of the string
	if (tp.find(GetHeaderLine(readerstep+adjust))!=std::string::npos){
	  if (petalIndex==70) cout << " Found matching string:  " << tp << "\n";
	  int parserStep=readerstep-11; //shift to the zero of the stripes
	  if (parserStep<0) parserStep+=4;//if we're negative, we must've been in the circles
	  if (readerstep<11) {
	    FillStepData(&txtfile,readerstep,
			 surveyDiameter, surveyPosition);
	     if (petalIndex==70) printf("Diameter (%d)=%f, (x,y)=(%f,%f)\n",parserStep,surveyDiameter[parserStep],
	    	   surveyPosition[parserStep].X(),surveyPosition[parserStep].Y());
	  } else {
	    FillStepData(&txtfile,readerstep,
			 stripeArea, stripePosition);
	    if (petalIndex==21 && readerstep==37) adjust++; //step 38 outlines more of stripe 37, but does not write to the text file as of Aug 2 (this is petal 71 after cleaning)
	    if (petalIndex==22 && readerstep==106){ //step 107 outlines more of stripe 106, and DOES write to the text file as of Aug 2 (this is petal 82 after cleaning, with step added)	  
		adjust++;
		AddStepData(&txtfile,readerstep,
			     stripeArea, stripePosition);
	    }
	    if (petalIndex==23 && readerstep==91) adjust++; //step 92 outlines more of stripe 91, but does not write to the text file as of Aug 2 (this is petal 73 after cleaning)
	    if (petalIndex==70 && readerstep==171) adjust++; //step 172 outlines more of stripe 171, but does not write to the text file as of July 29;
	
	    if (petalIndex==71 && readerstep==37) {
	      //step 37 fails, but the record does not have an extra step, so assume step 38 is the usual.  Do nothing (record this for posterity) as of July 29
	    }
	    if (petalIndex==81) printf("Area (%d)=%f, (x,y)=(%f,%f) (readerstep+adjust=%d)\n",parserStep,stripeArea[parserStep],
				       stripePosition[parserStep].X(),stripePosition[parserStep].Y(),readerstep+adjust);
	  }
	  readerstep++;
	}
      }
      txtfile.close(); //close the file object.
      //and now the arrays are filled with stripe and petal data.

      petal=petalIndex;
   
      pTree->Fill();

      for (int i=0;i<213;i++){
	//already done: petal=petalIndex;
	index=i;
	tpos=stripePosition[i];
	tarea=stripeArea[i];
	stripeTree->Fill();
      }

      for (int i=0;i<4;i++){
	index=i;
	tpos=surveyPosition[i];
	tdiam=surveyDiameter[i];
	surveyTree->Fill();
      }
    }
  }

  pTree->Write();
  stripeTree->Write();
  surveyTree->Write();
  petalFile->Close();
  return;
}

void ParseCircles(){
  printf("option <0 ==>  parsing circles\n");
  
  TFile *petalFile=TFile::Open(Form("%sall_petals.hist.root",directory.c_str()),"READ");
  TTree *surveyTree=(TTree*)petalFile->Get("surveyTree");
  
  
  //step 1:  extract arrays for all survey dots (4*(less than 1kb)*(less than 50) = small enough to hold in memory)
  const int MAXPETALS=100;
  const int MAXROWS=4;
  std::array<std::array<float,MAXROWS>,MAXPETALS> surveyDiam;
  std::array<std::array<TVector2,MAXROWS>,MAXPETALS> surveyPos;
  std::array<bool,MAXPETALS> petalExists;

  int nGoodPetals=0;
  for (int i=0;i<MAXPETALS;i++){
    surveyTree->Draw("i:pos.X():pos.Y():diam",Form("petal==%d",i),"goff");
    int nRows=surveyTree->GetSelectedRows();
    //printf("parsing petal %d (nrows=%d)\n",i,nRows);
    petalExists[i]=(nRows>0);
    if (petalExists[i]){
      nGoodPetals++;
    }
    for (int j=0;j<nRows;j++){
      surveyDiam[i][j]=surveyTree->GetVal(3)[j];
      surveyPos[i][j].Set(surveyTree->GetVal(1)[j],surveyTree->GetVal(2)[j]);
    }
  }


  //step 1a:  find the lowest petal number and define this to be our master petal:
  int basepetal=0;
  for (int i=0;i<MAXPETALS;i++){
    if (petalExists[i]) {
      printf("base petal is petal %d\n",i);

      basepetal=i;
      break;
    }
  }

    
    TH2F *hRelPosDiam=new TH2F("hRelPosDiam",Form("Survey Circle Size vs delta pos to petal #%d;size;deltaR",basepetal),50,0,1.4,200,0,3);
   TH1F *hDiam=new TH1F("hDiam",Form("Survey Dot Diameters;diameter"),50,8,12.);
   TH2F *hRelDiam=new TH2F("hRelDiam",Form("Survey Dot Diameters per petal;petal;diam"),50,45.5,95.5,100,8,12.);
   TH2F *hRelDiamn=new TH2F("hRelDiamn",Form("Survey Dot Diameters per corner;corner;diam"),4,-0.5,3.5,100,8,12.);
   TH2F *hRelPosX=new TH2F("hRelPosX",Form("Survey Positions per petal, rotated and relative to petal #%d;petal;deltaX",basepetal),50,45.5,95.5,100,-3,3.);
   TH2F *hRelPosY=new TH2F("hRelPosY",Form("Survey Positions per petal, rotated and relative to petal #%d;petal;deltaY",basepetal),50,45.5,95.5,100,-3,3.);
   TH2F *hRelPosXn=new TH2F("hRelPosXn",Form("Stripe Relative Positions vs stripe number, relative to petal #%d;stripe;deltaX",basepetal),213,-0.5,212.5,100,-3,3.);
   TH2F *hRelPosYn=new TH2F("hRelPosYn",Form("Stripe Relative Positions vs stripe number, relative to petal #%d;stripe;deltaY",basepetal),213,-0.5,212.5,100,-3,3.);

  for (int i=0;i<MAXPETALS;i++){
     if (!petalExists[i]) continue; //skip petals that aren't real;
     for (int j=0;j<MAXROWS;j++){
       hDiam->Fill(surveyDiam[i][j]);
       hRelPosDiam->Fill(surveyDiam[i][j],surveyPos[i][j].Mod());
       hRelDiam->Fill(i,surveyDiam[i][j]);
       hRelDiamn->Fill(j,surveyDiam[i][j]);
       hRelPosX->Fill(i,surveyPos[i][j].X());
       hRelPosY->Fill(i,surveyPos[i][j].Y());
       hRelPosXn->Fill(j,surveyPos[i][j].X());
       hRelPosYn->Fill(j,surveyPos[i][j].Y());
     }
   }


  for (int i=0;i<MAXPETALS;i++){
    int nBadSize=0;
    int nWarnSize=0;
    int nBadPos=0;
    int nWarnPos=0;
    if (!petalExists[i]) continue; //skip petals that aren't real;
    for (int j=0;j<MAXROWS;j++){
      if (j==2) continue; //skip the one that's missing the area
      int badSize=0;
      if (surveyDiam[i][j]<9.8 || surveyDiam[i][j]>10.6){
	badSize=2;
	nBadSize++;
	printf("Petal %d circle %d is bad (diam=%f\n",i,j,surveyDiam[i][j]);
      } else if (surveyDiam[i][j]<9.9 || surveyDiam[i][j]>10.5){
	badSize=1;
	nWarnSize++;
      }
      int badPos=0;
      /*
      if (relRotPos[i][j].Mod()>0.5){
	badPos=2;
	nBadPos++;
	printf("Petal %d stripe %d is bad (delX=%f, delY=%f)\n",i,j,relRotPos[i][j].X(),relRotPos[i][j].Y());

      } else if (relRotPos[i][j].Mod()>0.25){
	badPos=1;
	nWarnPos++;
      }
      hBad2D->Fill(i,badSize>badPos?badSize:badPos,1);
      hBad2DSize->Fill(i,badSize,1);
      hBad2DPos->Fill(i,badPos,1);
      */
    }
    int total=nBadSize+nWarnSize+nBadPos+nWarnPos;
    if (total>0){
      printf("Petal %d:  Tot %d  (BadSize=%d, WarnSize=%d, BadPos=%d, WarnPos=%d)\n",i,total,nBadSize,nWarnSize,nBadPos,nWarnPos);
    }
  }
  
  TCanvas *cSummary=new TCanvas("cSummary","cSummary",1000,1000);
  cSummary->Divide(3,1);
  cSummary->cd(1);
  hDiam->Draw();
    cSummary->cd(2);
   hRelDiam->Draw("colz");
  cSummary->cd(3);
  hRelDiamn->Draw("colz");
  return;
}

void ParseOGP(int q){
  printf("option %d = read in the data file and do more complicated analysis\n",q);
  if (q==2){
    ParseCircles();
    return;
  }
  int onlyDoPetal=0;
  if (q>2) onlyDoPetal=q;
  TFile *petalFile=TFile::Open(Form("%sall_petals.hist.root",directory.c_str()),"READ");
  TTree *stripeTree=(TTree*)petalFile->Get("stripeTree");
  
  /*
    a quick way to pull arrays out of a ttree:
    tree->Draw("v0:v1:...etc","1","goff");
    int nRows=tree->GetSelectedRows(); //gets the number of elements in the arrays.
    float v0_temp=tree->GetVal(0)[i]; //gets the ith element of the array of first argument of the draw command
    float v1_temp=tree->GetVal(1)[i]; //ditto for the second argument, etc.  But note that these get overwritten by the next draw command

    vector<float> v0,v1;
    for (int i=0;i<nRows;i++){
    v0.push_back(tree->GetVal(0)[i]);
    v1.push_back(tree->GetVal(1)[i]);
    }
    //now you have vectors (like arrays) that have the desired data permanently
    */
  
  //step 1:  extract arrays for all stripes (213*(less than 1kb)*(less than 50) = small enough to hold in memory)
  const int MAXPETALS=100;
  const int MAXROWS=213;
  std::array<std::array<float,MAXROWS>,MAXPETALS> stripeArea;
  std::array<std::array<float,MAXROWS>,MAXPETALS> relArea;
  std::array<std::array<TVector2,MAXROWS>,MAXPETALS> stripePos;
  std::array<std::array<TVector2,MAXROWS>,MAXPETALS> relPos,localPos,localRotPos, relRotPos;
  std::array<bool,MAXPETALS> petalExists;

  std::array<float,MAXROWS> meanArea;
  std::array<TVector2,MAXROWS> meanRotPos;
  
  int nGoodPetals=0;
  for (int i=0;i<MAXPETALS;i++){
    stripeTree->Draw("i:pos.X():pos.Y():area",Form("petal==%d",i),"goff");
    int nRows=stripeTree->GetSelectedRows();
    printf("parsing petal %d (nrows=%d)\n",i,nRows);
    petalExists[i]=(nRows>0);
    if (petalExists[i]){
      nGoodPetals++;
    }
    for (int j=0;j<nRows;j++){
      stripeArea[i][j]=stripeTree->GetVal(3)[j];
      stripePos[i][j].Set(stripeTree->GetVal(1)[j],stripeTree->GetVal(2)[j]);
    }
  }


  //step 1a:  find the lowest petal number and define this to be our master petal:
  int basepetal=0;
  for (int i=0;i<MAXPETALS;i++){
    if (petalExists[i]) {
      printf("base petal is petal %d\n",i);

      basepetal=i;
      break;
    }
  }

  //step 2: convert these into 'relative coordinates' by subtracting off the position of the base stripe, and scaling to that area
  for (int i=0;i<MAXPETALS;i++){
    if (!petalExists[i]) continue; //skip petals that aren't real;
    for (int j=0;j<MAXROWS;j++){
      localPos[i][j]=stripePos[i][j]-stripePos[i][0];
      relPos[i][j]=stripePos[i][j]-stripePos[i][0]-(stripePos[basepetal][j]-stripePos[basepetal][0]);
      relArea[i][j]=1;
      if (stripeArea[basepetal][j]>0) relArea[i][j]=stripeArea[i][j]/stripeArea[basepetal][j];
    }
  }

  //step 3: rotate so that the localRotatedPosition (localRotPos) is aligned for all petals:
  TH2F *hArea0to212=new TH2F("hArea0to212","Area comparison;area0;area212",50,7.5,10,50,14,16);
  TH2F *hAreaRatioToAngle=new TH2F("hAreaRatioToAngle","Area Ratio to Angle comparison;ratio;angle",100,6,10,100,13,16);
    if (build_gnuplot){
      printf("#note:  this is saved to my repo for posterity in its current form, but can also be generated by ParseOGP.C if you set the 'build_gnuplot' flag to true\n");
      printf("set terminal pdf size 60,30 linewidth 1 font 'Verdana,100'\n");
      printf("set output 'all_lines_large.pdf'\n");
    }
  for (int i=0;i<MAXPETALS;i++){
    if (!petalExists[i]) continue; //skip petals that aren't real;
    TVector2 axis=localPos[i][212]-localPos[i][0]; //farthest away.
    float phi=localPos[i][212].Phi(); //farthest away.  its angle is angle wrt x, since [i][0] is at (0,0) by construction.
    if (!build_gnuplot){
      //  printf("Petal %d (area0=%2.2f,area212=%2.2f,dist=%2.2f:  dx%d=%2.2f,dy%d=%2.2f,dphi%d=%2.4f\n",
      //	   i,stripeArea[i][0],stripeArea[i][212],axis.Mod(),
      //	   i,-1*stripePos[i][0].X(),i,-1*stripePos[i][0].Y(),i,-1*phi);
    } else{
        printf("area%d_0=%2.2f;area%d_212=%2.2f;dist%d=%2.2f;  dx%d=%2.2f;dy%d=%2.2f;dphi%d=%2.4f\n",
      	   i,stripeArea[i][0],i,stripeArea[i][212],i,axis.Mod(),
      	   i,-1*stripePos[i][0].X(),i,-1*stripePos[i][0].Y(),i,-1*phi);
    }
    if (i>40) hArea0to212->Fill(stripeArea[i][0],stripeArea[i][212]);
    for (int j=0;j<MAXROWS;j++){
      localRotPos[i][j]=localPos[i][j].Rotate(-phi);
      relRotPos[i][j]=localRotPos[i][j]-localRotPos[basepetal][j];
    }
  }
  // if (!build_gnuplot) hArea0to212->Draw("colz");

  if (build_gnuplot){
    for (int i=0;i<MAXPETALS;i++){
      if (!petalExists[i]) continue; //skip petals that aren't real;
      printf("plot  \"<grep 'Circle     7' '%spetal%d_edges.DAT' -A 350000\" u (($1+dx%d)*cos(dphi%d)-($2+dy%d)*sin(dphi%d)):(($2+dy%d)*cos(dphi%d)+($1+dx%d)*sin(dphi%d)) w l title 'petal %d'\n",directory.c_str(),i,
	     i,i,i,i,i,i,i,i,i);
    }
    printf("\n");

    //make one plot with ALL of them overlapped:
    bool firstLine=true;
    for (int i=0;i<MAXPETALS;i++){
      if (!petalExists[i]) continue; //skip petals that aren't real;
      if (firstLine){
	printf("plot ");
	firstLine=false;
      }else{
	printf(",\\\n");
      }
      printf("  \"<grep 'Circle     7' '%spetal%d_edges.DAT' -A 350000\" u (($1+dx%d)*cos(dphi%d)-($2+dy%d)*sin(dphi%d)):(($2+dy%d)*cos(dphi%d)+($1+dx%d)*sin(dphi%d)) w l notitle",directory.c_str(),i,
	     i,i,i,i,i,i,i,i);
    }
    printf("\n");
    exit;
  }

  
  //setp 4:  we know [0] is empty.  let's make that our actual average, now that we've aligned:
  for (int j=0;j<MAXROWS;j++){
    localRotPos[0][j].Set(0.,0.);
    stripeArea[0][j]=0;
  }

  //step 4a:  fill the average data.
  for (int i=0;i<MAXPETALS;i++){
    if (!petalExists[i] || i==70 || i==71) continue; //skip petals that aren't real, skip petal 70 because its positions are suspect, skip #71 because stripe 26 failed.
    for (int j=0;j<MAXROWS;j++){
      localRotPos[0][j]=localRotPos[0][j]+localRotPos[i][j];
      if (stripeArea[i][j]>20) printf("stripe %d in petal %d is suspect:  area=%f\n",j,i,stripeArea[i][j]);
      stripeArea[0][j]=stripeArea[0][j]+stripeArea[i][j];
    }
  }
   for (int j=0;j<MAXROWS;j++){
    localRotPos[0][j]=localRotPos[0][j]/(1.*(nGoodPetals-2));
    stripeArea[0][j]=stripeArea[0][j]/(1.*(nGoodPetals-2));//removing petal 70 and 71 from consideration
   }

 
  

   if(1){//redo our relRotPos and stripeArea using our average instead of #47:
     /*
       TH2F *normal=new TH2F("normal","normal;x;y",100,-100,600,100,-100,300);
       TH2F *mean=new TH2F("mean","mean;x;y",100,-100,600,100,-100,300);
       TCanvas *cHuh=new TCanvas("huh");
       cHuh->Divide(2,1);
     */

     /*
       mean->Fill(localRotPos[0][j].X(),localRotPos[0][j].Y());
       normal->Fill(localRotPos[47][j].X(),localRotPos[47][j].Y());
       if (j>10 && j<30) printf("Stripe %d mean area=%f normal area =%f\n",j,stripeArea[0][j],stripeArea[47][j]);
     */
     /*
       cHuh->cd(1);
       normal->Draw("colz");
       cHuh->cd(2);
       mean->Draw("colz");
     */
    
     basepetal=0;
     for (int i=0;i<MAXPETALS;i++){
       if (!petalExists[i]) continue; //skip petals that aren't real;
       for (int j=0;j<MAXROWS;j++){
	 relRotPos[i][j]=localRotPos[i][j]-localRotPos[basepetal][j];
	 relArea[i][j]=1.9;
	 if (stripeArea[basepetal][j]>0) relArea[i][j]=stripeArea[i][j]/stripeArea[basepetal][j];
       }
     }
   }

    
    TH2F *hRelPosSize=new TH2F("hRelPosSize",Form("Stripe Relative Size vs delta pos to petal #%d;size;deltaR",basepetal),50,0,1.4,200,0,3);
   TH2F *hRelSizes=new TH2F("hRelSizes",Form("Stripe Sizes per petal, scaled to petal #%d;petal;size",basepetal),50,45.5,95.5,100,0,2.);
   TH1F *hRelSize1D=new TH1F("hRelSize1D",Form("Stripe Relative Size, scaled to petal #%d;size;counts",basepetal),150,0,2.);
   TH2F *hRelSizen=new TH2F("hRelSizen",Form("Stripe Sizes, scaled to petal #%d;stripe;size",basepetal),213,-0.5,212.5,100,0,2.);
   TH2F *hRelPosX=new TH2F("hRelPosX",Form("Stripe Positions per petal, rotated and relative to petal #%d;petal;deltaX",basepetal),50,45.5,95.5,100,-3,3.);
   TH2F *hRelPosXY=new TH2F("hRelPosXY",Form("Stripe Positions per petal, rotated and relative to petal #%d;deltaX;deltaY",basepetal),100,-3,3,100,-3,3.);
   TH2F *hRelPosY=new TH2F("hRelPosY",Form("Stripe Positions per petal, rotated and relative to petal #%d;petal;deltaY",basepetal),50,45.5,95.5,100,-3,3.);
   TH2F *hRelPosXn=new TH2F("hRelPosXn",Form("Stripe Relative Positions vs stripe number, relative to petal #%d;stripe;deltaX",basepetal),213,-0.5,212.5,100,-3,3.);
   TH2F *hRelPosYn=new TH2F("hRelPosYn",Form("Stripe Relative Positions vs stripe number, relative to petal #%d;stripe;deltaY",basepetal),213,-0.5,212.5,100,-3,3.);

  for (int i=0;i<MAXPETALS;i++){
     if (!petalExists[i]) continue; //skip petals that aren't real;
     for (int j=0;j<MAXROWS;j++){
       hRelPosSize->Fill(relArea[i][j],relRotPos[i][j].Mod());
       hRelSizes->Fill(i,relArea[i][j]);
       hRelSize1D->Fill(relArea[i][j]);
       hRelSizen->Fill(j,relArea[i][j]);
       hRelPosX->Fill(i,relRotPos[i][j].X());
       hRelPosXY->Fill(relRotPos[i][j].X(),relRotPos[i][j].Y());
       hRelPosY->Fill(i,relRotPos[i][j].Y());
       hRelPosXn->Fill(j,relRotPos[i][j].X());
       hRelPosYn->Fill(j,relRotPos[i][j].Y());
     }
   }


  TH2F *hBad2D=new TH2F("hBad2D","N Bad (2), Worry (1), Okay (0);Badness;Petal",50, 45.5,95.5,2,0.5,2.5);
  TH2F *hBad2DSize=new TH2F("hBad2DSize","NBad (2), Worry (1), Okay (0) vs Petal  (size only )",50, 45.5,95.5,2,0.5,2.5);
  TH2F *hBad2DPos=new TH2F("hBad2DPos","N Bad (2), Worry (1), Okay (0) vs Petal (pos only)",50, 45.5,95.5,2,0.5,2.5);

  std::array<int,213> stripeBadSize;
  for (int i=0;i<213;i++){
    stripeBadSize[i]=0;
  }
  if (1){
  for (int i=0;i<MAXPETALS;i++){
    if (onlyDoPetal>0 && i!=onlyDoPetal) continue;
    int nBadSize=0;
    int nWarnSize=0;
    int nBadPos=0;
    int nWarnPos=0;
    if (!petalExists[i]) continue; //skip petals that aren't real;
    for (int j=0;j<MAXROWS;j++){
      if (j==2) continue; //skip the one that's missing the area
      int badSize=0;
      if (relArea[i][j]<0.9 || relArea[i][j]>1.25){
	badSize=2;
	stripeBadSize[j]++;
	nBadSize++;
	printf("Petal %d stripe %d is bad (relArea=%f, absAre=%f\n",i,j,relArea[i][j],stripeArea[i][j]);
      } else if (relArea[i][j]<0.9 || relArea[i][j]>1.1){
	//	badSize=1;
	//	nWarnSize++;
      }
      int badPos=0;
      if (relRotPos[i][j].Mod()>1.5){
	badPos=2;
	nBadPos++;
	printf("Petal %d stripe %d is bad (delX=%f, delY=%f, absX=%f, absY=%f)\n",i,j,relRotPos[i][j].X(),relRotPos[i][j].Y(), stripePos[i][j].X(),stripePos[i][j].Y());

      } else if (relRotPos[i][j].Mod()>0.25){
	//badPos=1;
	//nWarnPos++;
      }
      hBad2D->Fill(i,badSize>badPos?badSize:badPos,1);
      hBad2DSize->Fill(i,badSize,1);
      hBad2DPos->Fill(i,badPos,1);
    }
    int total=nBadSize+nWarnSize+nBadPos+nWarnPos;
    if (total>0){
      printf("Petal %d:  Tot %d  (BadSize=%d, WarnSize=%d, BadPos=%d, WarnPos=%d)\n",i,total,nBadSize,nWarnSize,nBadPos,nWarnPos);
    }
  }
  }
  
  TCanvas *cSummary=new TCanvas("cSummary","cSummary",1000,1000);
  cSummary->Divide(3,3);
  cSummary->cd(1);
  hRelSizes->Draw("colz");
   cSummary->cd(2);
   hRelSizen->Draw("colz");
  cSummary->cd(4);
  hRelPosX->Draw("colz");
  cSummary->cd(5);
  hRelPosXn->Draw("colz");
  cSummary->cd(7);
  hRelPosY->Draw("colz");
 cSummary->cd(8);
  hRelPosYn->Draw("colz");
  cSummary->cd(3);
 hRelPosSize->Draw("colz");
  cSummary->cd(6);
  hRelPosXY->Draw("colz");
  cSummary->cd(9);
  hRelSize1D->Draw("colz");
  
  return;
}
