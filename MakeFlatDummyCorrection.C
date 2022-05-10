//Makes a file containing the necessary histograms to fake a distortion that is the same value everywhere
//needs bounds such that all valid positions in the tpc are in the center of the 3x3x3 bins so that interpolation works everwhere.

void MakeFlatDummyCorrection(const char* filebase="dummyCorrection", const float rval=0,const float phival=0, const float zval=0){
  printf("Writing correction file %s with correction =(%2E,%2E,%2E)(r,phi,z)\n",filename,rval,phival,zval);
  enum coord{R=0,P=1,Z=2};
  TH3F *hFlat[2][3];
  float val[3];
  val[R]=rval;
  val[P]=phival;
  val[Z]=zval;

  float clow[3],chigh[3];//edges of the center bin, to be safe
  float low[3],high[3]; //calculated edges of the low and high bins.

  clow[R]=18;
  chigh[R]=82;
  
  clow[P]=-3.2;
  chigh[P]=6.3;
  
  clow[Z]=-110;
  chigh[R]=110;

  for (int i=0;i<3;i++){
    float range=chigh[i]-clow[i];
    low[i]=clow[i]-range;
    high[i]=chigh[i]+range;
  }
  
  const std::array<const std::string,2> extension = {{ "_negz", "_posz" }};
  const std::array<const std::string,3> axis = {{ "R", "P", "Z" }};
  for (int i=0;i<2;i++){
    for (int j=0;j<3;j++){
      hFlat[i][j]=new TH3F(Form("hIntDistortion%s%s",axis,extension),
			   Form("Dummy Correction for %s (%s), %2E everywhere",axis[j].c_str(),extension[i].c_str(),val[i]),
			   /*phi in rads.*/3,low[P],high[P],
			   /*r in cm.*/3,low[R],high[R],
			   /*z in cm*/3,low[Z],high[Z]);
	//now fill the bins:
	for (int k=0;k<hFlat[i][j]->GetNcells(),k++){
	  hFlat[i][j]->SetBinContent(k,val[j]);
	}	  
    }
  }

  //now write these to the desired file:
  TFile *output=TFile::Open(Form("%s.hist.root",filename),"RECREATE");
  for (int i=0;i<2;i++){
    for (int j=0;j<3;j++){
      hFlat[i][j]->Write();
    }
  }
}
