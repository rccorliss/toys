#include "TVector3.h"
#include "CylindricalFieldSim.h"

CylindricalFieldSim::CylindricalFieldSim(float dr, float dphi, float dz,int r, int phi, int z, float vdr){
  Escale=1; Bscale=1;
  vdrift=vdr;
  dim.SetXYZ(1,0,0);
  dim.SetPerp(dr);
  dim.SetPhi(dphi);
  phispan=dphi;
  dim.SetZ(dz);
  nr=r;nphi=phi;nz=z;

  
  //create a grid with the specified dimensions
  Efield=new MultiArray<TVector3>(nr,nphi,nz);
  for (int i=0;i<Efield->Length();i++)
    Efield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid to compute the impact form each cell on each cell.
  Epartial=new MultiArray<TVector3>(nr,nphi,nz,nr,nphi,nz);
  for (int i=0;i<Epartial->Length();i++)
    Epartial->GetFlat(i)->SetXYZ(0,0,0);

  //and to hold the external fieldmap
  Eexternal=new MultiArray<TVector3>(nr,nphi,nz);
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,0);

  //ditto the external magnetic fieldmap
  Bfield=new MultiArray<TVector3>(nr,nphi,nz);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid of charges in those volumes.
  q=new MultiArray<float>(nr,nphi,nz);
  for (int i=0;i<q->Length();i++)
    *(q->GetFlat(i))=0;

  //define the step:
  step.SetXYZ(1,0,0);
  step.SetPerp(dim.Perp()/nr);
  step.SetPhi(phispan/nphi);
  step.SetZ(dim.Z()/nz);
  // printf("stepsize:  r=%f,phi=%f, wanted %f,%f\n",step.Perp(),step.Phi(),dr/r,dphi/phi);


  //define the external fields:
  setFlatFields(1.4/*tesla*/,200/*V/cm*/);
  return;
}

TVector3 CylindricalFieldSim::calc_unit_field(TVector3 at, TVector3 from){
  //note this is the field due to a fixed point charge in free space.
  //if doing cylindrical calcs with different boundary conditions, this needs to change.
 
  const float k=8.987*1e13;//N*cm^2/C^2 in a vacuum. N*cm^2/C units, so that we supply space charge in coulomb units.
  TVector3 delr=at-from;
  TVector3 field=delr; //to set the direction.
  if (delr.Mag()<0.00001){
    //do nothing.  the vector is already zero, which will be our approximation.
    //field.SetMag(0);//no contribution if we're in the same cell. -- but root warns if trying to resize something of magnitude zero.
  } else{
    field.SetMag(k*1/(delr*delr));//scalar product on the bottom.
  }
  //printf("calc_unit_field at (%2.2f,%2.2f,%2.2f) from  (%2.2f,%2.2f,%2.2f).  Mag=%2.4fe-9\n",at.x(),at.Y(),at.Z(),from.X(),from.Y(),from.Z(),field.Mag()*1e9);
  
  return field;
}

TVector3 CylindricalFieldSim::fieldIntegral(float zdest,TVector3 start){
  int dir=(start.Z()>zdest?-1:1);//+1 if going to larger z, -1 if going to smaller;

  int r=start.Perp()/step.Perp();
  int phi=start.Phi()/step.Phi();
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
      fieldInt+=Efield->Get(r,phi,i)*step.Z();
  }
    fieldInt-=Efield->Get(r,phi,zi)*(startz-zi*step.Z());//remove the part of the low end cell we didn't travel through
    fieldInt+=Efield->Get(r,phi,zf)*(endz-zf*step.Z());//add the part of the high end cell we did travel through

    
  return dir*fieldInt;
}

TVector3 CylindricalFieldSim::interpolatedFieldIntegral(float zdest,TVector3 start){
  int dir=(start.Z()>zdest?-1:1);//+1 if going to larger z, -1 if going to smaller;

  float r0=start.Perp()/step.Perp();
  float phi0=start.Phi()/step.Phi();
  int zi,zf;
  double startz,endz;
  // printf("interpolating fieldInt at  r=%f,phi=%f\n",r0,phi0);

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
  //we can always interpolate between four positions.
  //wherever x is, only one of (x+dx/2, x-dx/2) will be in a different cell, where dx is the cell size.
  //the same logic applies to y.
  //so the four cells we want to interpolate with are the ones containing:
  //(x-dx/2,y-dy/2), (x-dx/2,y-dy/2), (x+dx/2,y-dy/2), (x+dx/2,y-dy/2),
  //where x is the nominal coordinate.
  //and the coefficients ought to be (1-|x-dx/2|/dx)(1-|y-dy/2|/dy), where x and y are the nominal coordinate measured wrt the lower corner of the box.
  //which is 1 if we're in the exact center, and 0.5 if we're on the edge in one direction, or 0.25 in the corner.
  float rf[4],phif[4]; //integer=what cell we're in.  fraction=how far into the cell we are
  rf[0]=rf[1]=(r0-0.5);
  rf[2]=rf[3]=(r0+0.5);
  phif[0]=phif[2]=(phi0-0.5);
  phif[1]=phif[3]=(phi0+0.5);
  float coef[4];
  int r[4], phi[4];

  
  for (int i=0;i<4;i++){
    if (rf[i]>=nphi){//handle wrap-around:
      rf[i]-=nphi;
      phif[i]-=2*TMath::Pi();
    }
    if (rf[i]<0){//handle wrap-around:
      rf[i]+=nphi;
      phif[i]+=2*TMath::Pi();
    }
    
    r[i]=rf[i]; //get integer portion of float
    phi[i]=phif[i];
 

    //  coef[i]=(1-abs( ( (start.X()-x[i]*step.X())-step.X()/2.0 )/step.X() )  ... can simplifiy:
    float rrel=r0-r[i]; //fractional r position of point in ith cell, from -0.5 to +1.5
    float phirel=phi0-phi[i];
    coef[i]=(1-abs( rrel - 0.5))*(1-abs( phirel - 0.5));
  }
  
  TVector3 fieldInt(0,0,0);
  TVector3 partialInt(0,0,0);
  for (int c=0;c<4;c++){
    partialInt.SetXYZ(0,0,0);
    //printf("looking for element r=%d,phi=%d\n",r[c],phi[c]);
    for(int i=zi;i<zf;i++){ //count the whole cell of the lower end, and skip the whole cell of the high end.
      
      partialInt+=Efield->Get(r[c],phi[c],i)*step.Z();
    }
  
    partialInt-=Efield->Get(r[c],phi[c],zi)*(startz-zi*step.Z());//remove the part of the low end cell we didn't travel through
    partialInt+=Efield->Get(r[c],phi[c],zf)*(endz-zf*step.Z());//add the part of the high end cell we did travel through
    fieldInt+=coef[c]*partialInt;
  }
    
  return dir*fieldInt;
}

void CylindricalFieldSim::populate_fieldmap(){
//sum the E field at every point in the nr by nphi by nz grid of field points.
//could also do this 'flat', but doesn't save any time
  //printf("in pop_fieldmap, n=(%d,%d,%d)\n",nr,ny,nz);
  for (int ir=0;ir<nr;ir++){
    for (int iphi=0;iphi<nphi;iphi++){
      for (int iz=0;iz<nz;iz++){
	//*f[ix][iy][iz]=sum_field_at(nr,ny,nz,partial_,q_,ix,iy,iz);
	TVector3 localF=sum_field_at(ir,iphi,iz);
	Efield->Set(ir,iphi,iz,localF);
	//if (localF.Mag()>1e-9)
	//printf("fieldmap@ (%d,%d,%d) mag=%f\n",ir,iphi,iz,localF.Mag());
      }
    }
  }
  return;
}

void  CylindricalFieldSim::populate_lookup(){
  //with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.

  //  TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  //printf("populating lookup for (%dx%dx%d)x(%dx%dx%d) grid\n",fx,fy,fz,ox,oy,oz);
  TVector3 at(1,0,0);
  TVector3 from(1,0,0);

  for (int ifr=0;ifr<nr;ifr++){
    at.SetPerp(step.Perp()*(ifr+0.5));
    for (int ifphi=0;ifphi<nphi;ifphi++){
      at.SetPhi(step.Phi()*(ifphi+0.5));
      for (int ifz=0;ifz<nz;ifz++){
	at.SetZ(step.Z()*(ifz+0.5));
	for (int ior=0;ior<nr;ior++){
	  from.SetPerp(step.Perp()*(ior+0.5));
	  for (int iophi=0;iophi<nphi;iophi++){
	    from.SetPhi(step.Phi()*(iophi+0.5));
	    for (int ioz=0;ioz<nz;ioz++){
	      from.SetZ(step.Z()*(ioz+0.5));
	      //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
	      //printf("calc_unit_field...\n");
	      Epartial->Set(ifr,ifphi,ifz,ior,iophi,ioz,calc_unit_field(at,from));
	    }
	  }
	}
      }
    }
  }
  return;

}

void CylindricalFieldSim::setFlatFields(float B, float E){
  //printf("CylindricalFieldSim::setFlatFields(B=%f,E=%f)\n",B,E);
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,E);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,B);
  return;
}

TVector3 CylindricalFieldSim::sum_field_at(int r,int phi, int z){
 //sum the E field over all nr by ny by nz cells of sources, at the specific position x,y,z.
  TVector3 sum(0,0,0);
  for (int ir=0;ir<nr;ir++){
    for (int iphi=0;iphi<nphi;iphi++){
      for (int iz=0;iz<nz;iz++){
	//sum+=*partial[x][phi][z][ix][iphi][iz] * *q[ix][iphi][iz];
	if (r==ir && phi==iphi && z==iz) continue;
	sum+=Epartial->Get(r,phi,z,ir,iphi,iz)*q->Get(ir,iphi,iz);
      }
    }
  }
  sum+=Escale*Eexternal->Get(r,phi,z);
  //printf("summed field at (%d,%d,%d)=(%f,%f,%f)\n",x,y,z,sum.X(),sum.Y(),sum.Z());
  return sum;
}

TVector3 CylindricalFieldSim::swimTo(float zdest,TVector3 start, bool interpolate=false){

 //using second order langevin expansion from http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-016.pdf
  //TVector3 (*field)[nr][ny][nz]=field_;

  //set the direction of the external fields.
  TVector3 B=Bfield->Get(0,0,0)*Bscale;//static field in tesla T=Vs/m^2
  
  //int x=start.X()/step.X();
  //int y=start.Y()/step.Y();
  //int zi=start.Z()/step.Z();
  double zdist=zdest-start.Z();


  TVector3 fieldInt;
  if (interpolate){
    fieldInt=interpolatedFieldIntegral(zdest,start);
  }else{
    fieldInt=fieldIntegral(zdest,start);
  }
  //float fieldz=field_[in3(x,y,0,fx,fy,fz)].Z()+E.Z();// *field[x][y][zi].Z();
  double fieldz=fieldInt.Z()/zdist;// average field over the path.
  
  double mu=vdrift/fieldz;//cm^2/(V*s);
  double omegatau=1*mu*B.Z();//mu*Q_e*B, units cm^2/m^2
  //originally the above was q*mu*B, but 'q' is really about flipping the direction of time.  if we do this, we negate both vdrift and q, so in the end we have no charge dependence -- we 'see' the charge by noting that we're asking to drift against the overall field.
  omegatau=omegatau*1e-4;//1m/100cm * 1m/100cm to get proper unitless.
  //printf("omegatau=%f\n",omegatau);
  double c0=1/(1+omegatau*omegatau);
  double c1=c0*omegatau;
  double c2=c1*omegatau;

  //really this should be the integral of the ratio, not the ratio of the integrals.
  //and should be integrals over the B field, btu for now that's fixed and constant across the region, so not necessary
  //there's no reason to do this as r phi.  This is an equivalent result, since I handle everything in vectors.
  double deltaX=c0*fieldInt.X()/fieldz+c1*fieldInt.Y()/fieldz-c1*B.Y()/B.Z()*zdist+c2*B.X()/B.Z()*zdist;
  double deltaY=c0*fieldInt.Y()/fieldz-c1*fieldInt.X()/fieldz+c1*B.X()/B.Z()*zdist+c2*B.Y()/B.Z()*zdist;
  double deltaZ=0; //not correct, but small?  different E leads to different drift velocity.  see linked paper.  fix eventually.

  printf("swimTo: (%2.4f,%2.4f,%2.4f) to z=%2.4f\n",start.X(),start.Y(), start.Z(),zdest);
  printf("swimTo: fieldInt=(%2.4f,%2.4f,%2.4f)\n",fieldInt.X(),fieldInt.Y(),fieldInt.Z());
  
  TVector3 dest(start.X()+deltaX,start.Y()+deltaY,zdest+deltaZ);
  
  return dest;
  
}
