#include <TH3.h>


const char* GetFileExtension(const char* filename);
const char* GetFileBase(const char* filename);
TH3 * GetMaskedHist(float rmin, float rmax, float phimin, float phimax, int z, TH3* input);
  
void MaskRegion(float rmin, float rmax, float phimin, float phimax, int z=1, char *filename="nope"){
  //take in a TH3, and mask out all bins corresponding to r0<r<r1 + phi0<phi<phi1, then write it to the provided output filename?

  //char *dirname=TSystem::DirName(filename);
  //char *fullfilename=TSystem::BaseName(filename);
  const char *basename=GetFileBase(filename);
  //printf("basename immediately after return:%s\n",basename);
  const char *extname=GetFileExtension(filename);

  char *outfilename=Form("%s%s%s",basename,".masked",extname);

  printf("input file: %s\noutput file:%s\n",filename,outfilename);
  
  TFile *file= TFile::Open(filename, "READ");
  TFile *outfile=TFile::Open(outfilename,"RECREATE");
  if (!file->IsOpen()) {
    std::cerr << "Error opening file: " << filename << std::endl;
    return;
  }

  TKey* key;
  TIter nextkey(file->GetListOfKeys());
  printf("found %d entries\n",file->GetListOfKeys()->GetEntries());
  int keynum=0;
  while ((key = dynamic_cast<TKey*>(nextkey()))) {
    printf("Trying key %d (%s)",keynum, key->GetTitle());
    keynum++;
    TObject* obj = key->ReadObj();
    if (obj->IsA() == TH3::Class()) {
      TH3* th3 = dynamic_cast<TH3*>(obj);
      outfile->cd();
      TH3* outhist=GetMaskedHist(rmin,rmax,phimin,phimax,z,th3);
      outfile->Write();
    }
    delete obj;
  }

  file->Close();
  outfile->Close();
  return;
}

const char* GetFileExtension(const char* filename) {
    TString name(filename);
    Ssiz_t dotIndex = name.Last('.');
    if (dotIndex != kNPOS) {
        return name.Data() + dotIndex;
    }
    return "";
}

const char* GetFileBase(const char* filename) {
    TString name(filename);
    TString substring;
    Ssiz_t dotIndex = name.Last('.');
    printf("last period is at %d\n",dotIndex);
    if (dotIndex != kNPOS) {
      substring=name(0, dotIndex);
      char *filebase=new char[dotIndex+1];
      strncpy(filebase,substring.Data(),dotIndex);
      filebase[dotIndex]='\0';
      printf("substring=%s\n",filebase);
      return filebase;
    }
    return filename;
}

TH3 * GetMaskedHist(float rmin, float rmax, float phimin, float phimax, int z, TH3* input){
  TH3 *output=(TH3*)input.Clone(); //copy
  //go through all the bins in the range, and zero the contents

  TAxis* ax[3];
  ax[0]= output->GetXaxis();
  ax[1]= output->GetYaxis();
  ax[2]= output->GetZaxis();

  float zbins[]={0,ax[2]->GetNBins()/2,ax[2]->GetLast};
  float zminbin,zmaxbin;
  if (z<0){
    zminbin=zbins[0];
    zmaxbin=zbins[1];
  } else {
    zminbin=zbins[1];
    zmaxbin=zbins[2];
  }
  
  for (int i=ax[0]->FindBin(phimin);i<=ax[0]->FindBin(phimax);i++){
    for (int j=ax[1]->FindBin(rmin);j<=ax[1]->FindBin(rmax);j++){
      for (int k=zminbin;k<=zmaxbin;k++){
	int glob=output->GetBin(i,j,k);
	output->SetBinContent(glob,0);
      }
    }
  }
};
