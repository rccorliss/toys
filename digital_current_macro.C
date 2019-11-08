void populate_lookup(int fx,int fy,int fz,int ox,int oy,int oz, TVector3 *field_, TVector3 dim);
void populate_fieldmap(int nx, int ny, int nz, TVector3 *field_,TVector3 *partial_,double *q_);
TVector3 sum_field_at(int nx, int ny,int nz, TVector3 *partial_,double *q_,int x,int y, int z);
TVector3 calc_unit_field(TVector3 at, TVector3 from);
TVector3 fieldIntegral(float zdest,TVector3 start,int fx, int fy, int fz, TVector3 *field_, TVector3 dim);
TVector3 swimTo(float zdest,TVector3 start, int q, int fx, int fy, int fz, TVector3 *field_, TVector3 dim);

#define in6(a,b,c,d,e,f,na,nb,nc,nd,ne,nf) (a*nb*nc*nd*ne*nf + b*nc*nd*ne*nf + c*nd*ne*nf + d*ne*nf + e*nf + f)
#define in3(a,b,c,na,nb,nc) (a*nb*nc+b*nc+c)

double Bmain;
double Emain;

void digital_current_macro(){
  printf("hello\n");


  //all dimensions in cm, Coulomb.
  
  //define a detector box.  Propagation is along z to z=0.
  const TVector3 box(20,20,100);

  //define external field strength:
  Bmain=1.4;//Tesla
  Emain=200;//V/cm

  
  //define a reconstruction field grid size and the charge holders for it:
  const int rnx=3, rny=3, rnz=6;
  printf("building reco grid (%d x %d x %d)\n",rnx,rny,rnz);
  double rQ[rnx][rny][rnz];

  //define the MULTIPLIERS in xyz for the actual generation field grid size, to make sure it factors neatly:
  const int genmx=2,genmy=2,genmz=2;
  int gnx=rnx*genmx, gny=rny*genmy, gnz=rnz*genmz;
  printf("building gen grid (%d x %d x %d)\n",gnx,gny,gnz);
  double gQ[gnx][gny][gnz];
  //generate the gen lookup table for what the field at point a,b,c is due to charge at the center of gQ[i,j,k]
  //first three are the coordinates of the test charge, last three are the coordinates of the generating charge
  //  TVector3 gEpart[gnx][gny][gnz][gnx][gny][gnz];
  //seems root can't handle 6 dim arrays:
  int ncells=gnx*gny*gnz*gnx*gny*gnz;
  //max number of cells=208842;
  printf("asking for %d TVector3 array\n",ncells);
    TVector3 gEpart[ncells];

    printf("got %d TVector3 array\n",ncells);

     //populate_lookup(gnx,gny,gnz,gnx,gny,gnz,&(gEpart[0][0][0][0][0][0]),box);
  populate_lookup(gnx,gny,gnz,gnx,gny,gnz,&(gEpart[0]),box);



  //generate the reco lookup table for what the field at point a,b,c is due to charge at the center of rQ[i,j,k]
  TVector3 rEpart[rnx*rny*rnz*rnx*rny*rnz];
  
  populate_lookup(rnx,rny,rnz,rnx,rny,rnz,&(rEpart[0]),box);

  //return;
  
  //generate a charge distribution:
  TF1 *chargedist=new TF1("chargedist","gaus(0)",0,100);
  chargedist->SetParameters(1/*norm*/,0.1/*peak in C*/,0.05/*width in C*/);
  //fill our true charge grid by sampling the distribution
  for (int ix=0;ix<gnx;ix++){
    for (int iy=0;iy<gny;iy++){
      for (int iz=0;iz<gnz;iz++){
	gQ[ix][iy][iz]=0;//chargedist->GetRandom();
      }
    }
  }
  gQ[0][0][6]=1e-13; ///that's 10^6 protons in that box.

  //fill our measured charge grid by summing the true grid:
  //maybe add noise here later.
 for (int ix=0;ix<rnx;ix++){
    for (int iy=0;iy<rny;iy++){
      for (int iz=0;iz<rnz;iz++){
	rQ[ix][iy][iz]=0;
      }
    }
  }
 for (int ix=0;ix<gnx;ix++){
    for (int iy=0;iy<gny;iy++){
      for (int iz=0;iz<gnz;iz++){
	rQ[ix/genmx][iy/genmy][iz/genmz]+=gQ[ix][iy][iz];
      }
    }
  }

 //calculate our true gen field map
 TVector3 gEtot[gnx][gny][gnz];
 populate_fieldmap(gnx,gny,gnz,&(gEtot[0][0][0]),&(gEpart[0]),&(gQ[0][0][0]));

 //calculate our reco field map
 TVector3 rEtot[rnx][rny][rnz];
  populate_fieldmap(rnx,rny,rnz,&(rEtot[0][0][0]),&(rEpart[0]),&(rQ[0][0][0]));

 
 //create a particle
 TVector3 part(12.005,12.005,100);

 //swim it forward
 TVector3 forward=swimTo(0,part, +1,gnx, gny, gnz, &(gEtot[0][0][0]), box);
 TVector3 back=swimTo(part.Z(),forward, -1,gnx, gny, gnz,  &(gEtot[0][0][0]), box);

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
   
  
 float Bmin=0.0001;
 float Bmax=1.4;
 int nsteps=100;
 float exp=1;
 float bval[nsteps];
 float xoff[nsteps];
  float yoff[nsteps];

 for (int i=0;i<nsteps;i++){
   Bmain=Bmin+(Bmax-Bmin)*TMath::Power((i)/(1.0*nsteps),exp);//TMath::Exp(-100+i);
   forward=swimTo(0,part, -1,gnx, gny, gnz, &(gEtot[0][0][0]), box);
   back=swimTo(part.Z(),forward,+1,rnx,rny,rnz, &(rEtot[0][0][0]), box);
   bval[i]=Bmain;
   xoff[i]=back.X()-part.X();
   yoff[i]=back.Y()-part.Y();
 }
   TGraph *tgx=new TGraph(nsteps,bval,xoff);
   TGraph *tgy=new TGraph(nsteps,bval,yoff);
   TGraph *tgxy=new TGraph(nsteps,xoff,yoff);
   tgxy->SetTitle("corrected position residual for increasing Bfield;x(cm);y(cm)");
   tgy->SetTitle("corrected position residual for increasing Bfield;B(T);y(cm)");
    //tgx->Draw("AC*");
    tgy->Draw("AC*");
    //tgxy->Draw("AC*");
  return;
}

//void foo(int n, int m, int bar[n][m]);
void populate_fieldmap(int nx, int ny, int nz, TVector3 *field_,TVector3 *partial_,double *q_){
  //  TVector3 (*f)[nx][ny][nz];
  //  printf("%d\n",sizeof(f));
  //f= (TVector3 (*)[nx][ny][nz]) field_;
  //float (&f)[nx][ny][nz] = *reinterpret_cast<float (*)[nx][ny][nz]>(field_);

  
//sum the E field at every point in the nx by ny by nz grid of field points.
  for (int ix=0;ix<nx;ix++){
    for (int iy=0;iy<ny;iy++){
      for (int iz=0;iz<nz;iz++){
	//*f[ix][iy][iz]=sum_field_at(nx,ny,nz,partial_,q_,ix,iy,iz);
	TVector3 localf=sum_field_at(nx,ny,nz,partial_,q_,ix,iy,iz);
	field_[in3(ix,iy,iz,nx,ny,nz)]=localf;
	//if (localf.Mag()>1e-9)
	//	  printf("fieldmap@ (%d,%d,%d) mag=%f\n",ix,iy,iz,localf.Mag());
      }
    }
  }
}

TVector3 sum_field_at(int nx, int ny,int nz, TVector3 *partial_,double *q_,int x,int y, int z){
  //  TVector3 (*partial)[nx][ny][nz][nx][ny][nz]=partial_;
  //  float (*q)[nx][ny][nz]=q_;
  // printf("sum_field_at (%d,%d,%d)\n",x,y,z);

  
  //sum the E field over all nx by ny by nz cells of sources, at the specific position x,y,z.
  TVector3 sum(0,0,0);
  for (int ix=0;ix<nx;ix++){
    for (int iy=0;iy<ny;iy++){
      for (int iz=0;iz<nz;iz++){
	//sum+=*partial[x][y][z][ix][iy][iz] * *q[ix][iy][iz];
	if (x==ix && y==iy && z==iz) continue;
	sum+=partial_[in6(x,y,z,ix,iy,iz,nx,ny,nz,nx,ny,nz)]*q_[in3(ix,iy,iz,nx,ny,nz)];
      }
    }
  }
  return sum;
}


void  populate_lookup(int fx,int fy,int fz,int ox,int oy,int oz, TVector3 *field_, TVector3 dim){
  //with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.

  //  TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  //printf("populating lookup for (%dx%dx%d)x(%dx%dx%d) grid\n",fx,fy,fz,ox,oy,oz);
  TVector3 at(0,0,0);
  TVector3 from(0,0,0);
  TVector3 step(dim.X()/fx,dim.Y()/fy,dim.Z()/fz);

  for (int ifx=0;ifx<fx;ifx++){
    at.SetX(step.X()*(ifx+0.5));
    for (int ify=0;ify<fy;ify++){
      at.SetY(step.Y()*(ify+0.5));
      for (int ifz=0;ifz<fz;ifz++){
	at.SetZ(step.Z()*(ifz+0.5));

	for (int iox=0;iox<ox;iox++){
	  from.SetX(step.X()*(iox+0.5));
	  for (int ioy=0;ioy<oy;ioy++){
	    from.SetY(step.Y()*(ioy+0.5));
	    for (int ioz=0;ioz<oz;ioz++){
	      from.SetZ(step.Z()*(ioz+0.5));
	      //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
	      //printf("calc_unit_field...\n");
	      field_[in6(ifx,ify,ifz,iox,ioy,ioz,fx,fy,fz,ox,oy,oz)]=calc_unit_field(at,from);
	    }
	  }
	}
      }
    }
  }
  
  
  return;
}



TVector3 calc_unit_field(TVector3 at, TVector3 from){
  //note this is the field due to a fixed point charge in free space.
  //if doing cylindrical calcs with different boundary conditions, this needs to change.
 
  const float k=8.987*1e13;//N*cm^2/C^2 in a vacuum. N*cm^2/C units, so that we supply space charge in coulomb units.
  TVector3 delr=at-from;
  TVector3 field=delr; //to set the direction.
  if (delr.Mag()<0.00001){
    //do nothing.  the vector is already zero, which will be our approximation.
    //field.SetMag(0);//no contribution if we're in the same cell.
  } else{
    field.SetMag(k*1/(delr*delr));//scalar product on the bottom.
  }
  // printf("calc_unit_field at (%2.2f,%2.2f,%2.2f) from  (%2.2f,%2.2f,%2.2f).  Mag=%2.4fe-9\n",at.x(),at.Y(),at.Z(),from.X(),from.Y(),from.Z(),field.Mag()*1e9);
  
  return field;
}


TVector3 fieldIntegral(float zdest,TVector3 start,int fx, int fy, int fz, TVector3 *field_, TVector3 dim){
  //calculate the integrated, vectorial field along the path from start to zdest, assuming we don't cross a cell boundary
  //in the future, maybe interpolate between the values in the eight adjacent cells based on our position?  slightly less wrong, maybe?
  //TVector3 (*field)[nx][ny][nz]=field_;

  TVector3 step(dim.X()/fx,dim.Y()/fy,dim.Z()/fz);
  int dir=(start.Z()>zdest?-1:1);//+1 if going to larger z, -1 if going to smaller;

  int x=start.X()/step.X();
  int y=start.Y()/step.Y();
  int zi,zf;
  double startz,endz;

  //make sure 'zi' is always the smaller of the two numbers, for handling the partial-steps.
  if (dir>0){
    zi=start.Z()/step.Z(); //highest cell with lower bound less than position of particle
    zf=zdest/step.Z(); //highest cell with lower bound less than target position of particle
    startz=start.Z();
    endz=zdest;
  } else{
    zf=start.Z()/step.Z(); //highest cell with lower bound less than position of particle
    zi=zdest/step.Z(); //highest cell with lower bound less than target position of particle
    startz=zdest;
    endz=start.Z();
  } 

  TVector3 fieldInt(0,0,0);
    for(int i=zi;i<zf;i++){ //count the whole cell of the lower end, and skip the whole cell of the high end.
      fieldInt+=field_[in3(x,y,i,fx,fy,fz)]*step.Z();
  }
    fieldInt-=field_[in3(x,y,zi,fx,fy,fz)]*(startz-zi*step.Z());//remove the part of the low end cell we didn't travel through
    fieldInt+=field_[in3(x,y,zf,fx,fy,fz)]*(endz-zf*step.Z());//add the part of the high end cell we did travel through

    
  return dir*fieldInt;
  }


TVector3 swimTo(float zdest,TVector3 start, int q, int fx, int fy, int fz, TVector3 *field_, TVector3 dim){
  //using second order langevin expansion from http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-016.pdf
  //TVector3 (*field)[nx][ny][nz]=field_;

  //set the direction of the external fields.
  TVector3 B(0,0,1);//1.4);//static field in tesla T=Vs/m^2
  TVector3 E(0,0,1);//static field in V/cm.
  B.SetMag(Bmain);//set the magnitudes from global var
  E.SetMag(Emain);//set the magnitudes from global var

  TVector3 step(dim.X()/fx,dim.Y()/fy,dim.Z()/fz);
  int x=start.X()/step.X();
  int y=start.Y()/step.Y();
  int zi=start.Z()/step.Z();
  float zdist=zdest-start.Z();

  double vdrift=100.0/12e-6;// 100cm/12us = our approximate drift time.
  double mu=vdrift/E.Z();//cm^2/(V*s);
  double omegatau=1*mu*B.Z();//mu*Q_e*B, units cm^2/m^2
  //originally the above was q*mu*B, but 'q' is really about flipping the direction of time.  if we do this, we negate both vdrift and q, so in the end we have no charge dependence -- we 'see' the charge by noting that we're asking to drift against the overall field.
  omegatau=omegatau*1e-4;//1m/100cm * 1m/100cm to get proper unitless.
  printf("omegatau=%f\n",omegatau);
  double c0=1/(1+omegatau*omegatau);
  double c1=c0*omegatau;
  double c2=c1*omegatau;

  TVector3 fieldInt=fieldIntegral(zdest,start,fx,fy,fz,field_,dim);
  //float fieldz=field_[in3(x,y,0,fx,fy,fz)].Z()+E.Z();// *field[x][y][zi].Z();
  float fieldz=fieldInt.Z()/zdist+E.Z();// average field over the path.
  float deltaX=c0*fieldInt.X()/fieldz+c1*fieldInt.Y()/fieldz-c1*B.Y()/B.Z()*zdist+c2*B.X()/B.Z()*zdist;
  float deltaY=c0*fieldInt.Y()/fieldz-c1*fieldInt.X()/fieldz+c1*B.X()/B.Z()*zdist+c2*B.Y()/B.Z()*zdist;
  float deltaZ=0; //not correct, but small?  different E leads to different drift velocity.  see linked paper.  fix eventually.

  printf("swimTo: (%2.4f,%2.4f,%2.4f) (cell %d,%d,%d) to z=%2.4f\n",start.X(),start.Y(), start.Z(),x,y,zi,zdest);
  printf("swimTo: fieldInt=(%2.4f,%2.4f,%2.4f)\n",fieldInt.X(),fieldInt.Y(),fieldInt.Z());
  
  TVector3 dest(start.X()+deltaX,start.Y()+deltaY,zdest+deltaZ);
  
  return dest;
}
