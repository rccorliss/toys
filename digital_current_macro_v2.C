#include "FieldSim.h"
R__LOAD_LIBRARY(.libs/libfieldsim)

void digital_current_macro_v2(){
  printf("hello\n");
  //gSystem->Load("./libs/libfieldsim");
  //all dimensions in cm, Coulomb.
  //define a detector box.  Propagation is along z to z=0.
  const TVector3 box(20,20,100);

  //define external field strength:
  float Bmain=1.0;//Tesla.  nominal is 1.4, but now I can scale in units of tesla.
  float Emain=200;//V/cm
  double vdrift=100.0/12e-6;// 100cm/12us = our approximate drift time.
 
 
  //define a reco FieldSim
  const int rnx=8, rny=8, rnz=8;
  printf("building reco grid (%d x %d x %d)\n",rnx,rny,rnz);
  FieldSim *reco=new FieldSim(box.X(),box.Y(),box.Z(),rnx,rny,rnz,vdrift);
  reco->setFlatFields(Bmain,Emain);

  //define the MULTIPLIERS in xyz for the actual generation field grid size, to make sure it factors neatly:
  const int genmx=3,genmy=3,genmz=3;
  int gnx=rnx*genmx, gny=rny*genmy, gnz=rnz*genmz;
  printf("building gen grid (%d x %d x %d)\n",gnx,gny,gnz);
  FieldSim *gen=new FieldSim(box.X(),box.Y(),box.Z(),gnx,gny,gnz,vdrift);
  gen->setFlatFields(Bmain,Emain);

    
  //generate the lookup tables for what the field at point a,b,c is due to charge at the center of gQ[i,j,k]
  //only need to be done once per dimension/division.  Could be packed into the class but i want timing info.
  printf("populating reco lookup\n");
  reco->populate_lookup();
   printf("populating gen lookup\n");
 gen->populate_lookup();

 
  //generate a charge distribution:
  TF1 *chargedist=new TF1("chargedist","gaus(0)",0,100);
  chargedist->SetParameters(1/*norm*/,0.1/*peak in C*/,0.05/*width in C*/);
  MultiArray<float> *q=gen->q; 
  //fill our true charge grid by sampling the distribution
   printf("populating charge dist\n");
  for (int ix=0;ix<gnx;ix++){
    for (int iy=0;iy<gny;iy++){
      for (int iz=0;iz<gnz;iz++){
	q->Set(ix,iy,iz,0);//chargedist->GetRandom();
      }
    }
  }
  q->Set(0,0,0,1e-13);///that's 10^6 protons in that box.

  //fill our measured charge grid by summing the true grid:
  //maybe add noise here later.
  printf("summing for reco charge grid with q(0,0,6)=%fe-9\n",q->Get(0,0,0)*1e9);

 for (int ix=0;ix<gnx;ix++){
    for (int iy=0;iy<gny;iy++){
      for (int iz=0;iz<gnz;iz++){
	reco->q->Set(ix/genmx,iy/genmy,iz/genmz,
		     reco->q->Get(ix/genmx,iy/genmy,iz/genmz)+q->Get(ix,iy,iz));
      }
    }
  }

 //calculate our field map for the given charge configurations
 printf("populating gen fieldmap with given charge\n");

 gen->populate_fieldmap();
  printf("populating reco fieldmap with given charge\n");

 reco->populate_fieldmap();
 
 //create a particle
 TVector3 part(12.005,12.005,99.99);

 //swim it forward
 TVector3 forward=gen->swimTo(0,part,false);
 TVector3 back=gen->swimTo(part.Z(),forward,false);

 printf("part: (%2.4f,%2.4f,%2.4f)\n",part.X(),part.Y(),part.Z());
 // printf("forward: (%2.4f,%2.4f,%2.4f)\n",forward.X(),forward.Y(),forward.Z());
 TVector3 diff=part-forward;
 printf("orig-forward: (%2.4f,%2.4f,%2.4f)\n",diff.X(),diff.Y(),diff.Z());

 // printf("back: (%2.4f,%2.4f,%2.4f)\n",back.X(),back.Y(),back.Z());
 diff=forward-back;
 printf("forward-back: (%2.4f,%2.4f,%2.4f)\n",diff.X(),diff.Y(),diff.Z());
  diff=part-back;
 printf("orig-back: (%2.4fe-9,%2.4fe-9,%2.4f)\n",diff.X()*1e9,diff.Y()*1e9,diff.Z());

 //swim it backward
   

 float Bratiomin=0.0001;
 float Bratiomax=1.4;
 float Bratio;
 int nsteps=100;
 int zstepsforward=100;
 int zstepsback=20;
 float exp=1;
 float bval[nsteps];
 float xoff[nsteps];
  float yoff[nsteps];
  TVector3 midstep;
 
  for (int i=0;i<nsteps;i++){
    zstepsback=i+1;
    Bratio=1.4;
    //Bratio=Bratiomin+(Bratiomax-Bratiomin)*TMath::Power((i)/(1.0*nsteps),exp);//TMath::Exp(-100+i);
   gen->setScaleFactorB(Bratio);   
   reco->setScaleFactorB(Bratio);
   midstep=part;
   for (int j=0;j<zstepsforward;j++){
     midstep=gen->swimTo(part.Z()*(zstepsforward-1-j)/zstepsforward,midstep,true);
   }
   forward=midstep;
   for (int j=0;j<zstepsback;j++){
         midstep=reco->swimTo(part.Z()*(j+1)/zstepsback,midstep,false);
   }
   back=midstep;
   //forward=reco->swimTo(0,part,true);
   //back=reco->swimTo(part.Z(),forward,true);
   bval[i]=Bratio;
   xoff[i]=1e4*(back.X()-part.X());
   yoff[i]=1e4*(back.Y()-part.Y());
  }
   TGraph *tgx=new TGraph(nsteps,bval,xoff);
   TGraph *tgy=new TGraph(nsteps,bval,yoff);
   TGraph *tgxy=new TGraph(nsteps,xoff,yoff);
   TGraph *tgxys=new TGraph(1,xoff,yoff);
   TGraph *tgxye=new TGraph(1,xoff+nsteps-1,yoff+nsteps-1);
   tgxy->SetTitle("corrected position residual for increasing Bfield;x(cm);y(cm)");
   tgxys->SetMarkerColor(2);
   tgxys->SetMarkerStyle(kStar);
   tgxye->SetMarkerColor(3);
   tgxye->SetMarkerStyle(kStar);
   tgy->SetTitle("corrected position residual for increasing Bfield;B(T);y(cm)");
    //tgx->Draw("AC*");
    //tgy->Draw("AC*");
    tgxy->Draw("AC*");
    tgxys->Draw("same");

    TMultiGraph *mg=new TMultiGraph();
    mg->SetTitle(Form("xy offset for varying swim steps\n grids: gen(%dx%dx%d) reco(%dx%dx%d) zsteps=%d;x(um);y(um)",gnx,gny,gnz,rnx,rny,rnz,zstepsforward));
    mg->Add(tgxy);
    mg->Add(tgxys);
    mg->Add(tgxye);
    mg->Draw("AC*");
  return;
}
