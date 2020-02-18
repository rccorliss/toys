#include <iostream>
#include <vector>
#include <fstream>
#include <map>
#include "TIterator.h"
#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TString.h"
#include "TStyle.h"
#include "TObjString.h"
#include "TLegend.h"
#include "Rtypes.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TH2F.h"
#include "TColor.h"

using namespace std;


void CreateExternalFieldMap(){

  float range1 = .001; float range2 = .01; float range3 = .1;
  
  char name[100] = "1050_EF_Ideal";
  char dirname[100]="/Users/rcorliss/Downloads/Ross_Ideal_Distortion_Map/";
 
  char inFileNodes[175];
  sprintf(inFileNodes,"%sNLIST_%s.lis",dirname,name);
  char inFileElements[175];
  sprintf(inFileElements,"%sPRNSOL_%s.lis",dirname,name);

  int nY = 8;
  int nX = 40;
  int totalSamples = (nY-1)*(nX-1);
  float tol =2;

  //field cage width (mm)
  float height = 1090;//z
  float width = 800; //x
  float widthY = 1090*2;//y

  

  int ymin = 0.;
  int ymax = 1090;//560;
  int zmin = 0.;
  int zmax = 1090.;
  int xmin = 0.;
  int xmax = 800;
  printf("trying to open %s\n",inFileNodes);
  ifstream fileN(inFileNodes);
  if(!fileN.is_open()) cerr << "node list file not open" << endl;
  TString lineN;
  vector<int> nodeVec;
  vector<float> xVec;
  vector<float> yVec;
  vector<float> zVec;
  int lineCount = 0;

  map<int,pair<float,float> > nodeMap;
  cout << "checkpoint 1 " << endl;

  // Parse ansys NLIST (node list) txt files to get a list of coordinates for the measurements:
  while (lineN.ReadLine(fileN)) {
    // strip leading spaces and skip comments
    if (lineCount == 0 || lineCount == 1) {++lineCount; continue;}
    lineN.Remove(TString::kBoth, ' ');
    TObjArray *tokensN = lineN.Tokenize(" ");
    TIter myIterN(tokensN);
    int count = 0;
    if ( (int) tokensN->GetEntries() < 2) continue;
    while(TObjString *stN = (TObjString*) myIterN.Next()) {
      TString str=stN->GetString().Strip(TString::kBoth,' ');
      if (str.EqualTo("NODE",TString::kIgnoreCase)) break;
      if (count == 0) nodeVec.push_back( str.Atoi());
      else if (count == 1) xVec.push_back( str.Atof());
      else if (count == 2) yVec.push_back( str.Atof());
      else if (count == 3) zVec.push_back( str.Atof());
      ++count;
    }
    ++lineCount;
  }
  fileN.close();
 
  cout << "checkpoint 2 " << " r vec size " << xVec.size() << " z vec size " << yVec.size() << endl;
  int sampleCount = 0;

  //convert coordinates into coordinate pairs:  x coord== r position, y coord== z position.
  //put them in a map by their node ID.
  for (unsigned int z = 0; z < yVec.size(); ++z) {
    nodeMap[nodeVec[z]]= make_pair(yVec[z],xVec[z]);
    //cout << " node " << nodeVec[z] << endl;
  }


  cout << " X Nodes found " << (int) nodeMap.size() << endl;

  ifstream fileE(inFileElements);
  if(!fileE.is_open()) cerr << "file not open" << endl;
  
  vector<int> nodeVecE;
  vector<float> efXVec;
  vector<float> efYVec;
  vector<float> efZVec;
  vector<float> efSumVec;
  TString lineE;
  cout << "checkpoint 3 " << endl;
  int nodeCount = 0;
  bool startRecording = 0;
  int totalCount = 0;
  int elCount = 0;
  int sumCount = 0;
  int nCount=0;

  
  // read in element list PRNSOL txt file to get field values at coordinates
  while (lineE.ReadLine(fileE)) {
    // strip leading spaces and skip comments
    lineE.Remove(TString::kBoth, ' ');
    // Array of tokens separated by ","
    TObjArray *tokensE = lineE.Tokenize(" ");
    bool negWarning = 0;
    int lineSize = (int) tokensE->GetEntries() ;
    //    cout << "lineSize " << endl;
    TIter myIterE(tokensE);
    int nodeE;
    
 
    int count = 0;

    while(TObjString *stE = (TObjString*) myIterE.Next()) {
      TString str=stE->GetString().Strip(TString::kBoth,' ');

      if (nodeCount != 0 && (str).EqualTo(" " ,TString::kIgnoreCase)) {startRecording = 0; nodeCount = 0; break;}
      else if (nodeCount == 0 && str.EqualTo("NODE",TString::kIgnoreCase)) {startRecording = 1; break;}
      else if (startRecording == 0) break;
      
      if (count == 0 ) {
	 int node = str.Atoi(); 
	 //cout << " node " << node << " linesize " << lineSize << endl;
        if (node == 0) {
          startRecording = 0; nodeCount = 0; ++elCount; break; 
        }

        else if (lineSize == 6) {  nodeVecE.push_back(node);/* cout << "node " << nodeVecE[nCount];*/ ++nCount; ++nodeCount; }
      }
      //rcc efxVec and efyVec contain the x and y components (r and z, really)
      else if (count == 1 && lineSize == 6){ efXVec.push_back(str.Atof()); }
      else if (count == 2 && lineSize == 6 /*&& lineSize > 3*/) {efYVec.push_back(str.Atof());   }
      else if (count == 3 && lineSize == 6 /*&& lineSize > 4*/) {efZVec.push_back(str.Atof());  }
      else if (count == 4 && lineSize == 6 /*|| lineSize < 5*/){ efSumVec.push_back(str.Atof()); ++totalCount;  /*cout << " sum " << efSumVec[sumCount]<< endl;*/  ++sumCount;  }
      ++count;
    }
    if (startRecording == 1 && nodeCount == 36) {startRecording = 0; nodeCount = 0;++elCount; /*cout << " element " << elCount << endl; int check; cin >> check;*/}
  }
  fileN.close();
  
 
  cout << "checkpoint 4 " << " node vec E size " << nodeVecE.size() << endl;
  sampleCount=0;



  //create a tree that matches the format of the B field map:
  TFile output("externalEfield.ttree.root","RECREATE");
  TTree fTree("fTree","external field Tree");
  float rfield,zfield,rcoord,zcoord;
  fTree.Branch("r",&rcoord);
  fTree.Branch("z",&zcoord);
  fTree.Branch("er",&rfield);
  fTree.Branch("ez",&zfield);

  

  //loop over each node we found in the field measurements and match it to a node in the coordinates:
  for (int n = 0; n < (int) nodeVecE.size(); ++n) {
    // cout << " Node " << nodeVecE[n] << endl;
    if (nodeMap.find(nodeVecE[n]) != nodeMap.end() /*&& efSumMap.find(nodeVecE[n]) == efSumMap.end()*/) {
      rcoord=nodeMap[nodeVecE[n]].second/10.0;//convert to cm
      zcoord=nodeMap[nodeVecE[n]].first/10.0;//convert to cm
      rfield=efXVec[n]*10;//convert to V/cm
      zfield=efYVec[n]*10;//convert to V/cm
      fTree.Fill();
    }
  }

  fTree.Write();
  output.Close();
  return;
  
  
 
  }
