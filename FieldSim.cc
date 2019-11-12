#include "TVector3.h"
#include "FieldSim.h"

FieldSim::FieldSim(float dx, float dy, float dz,int x, int y, int z, float vdr){
  Escale=1; Bscale=1;
  vdrift=vdr;
  dim.SetXYZ(dx,dy,dz);
  nx=x;ny=y;nz=z;
  
  //create a grid with the specified dimensions
  Efield=new MultiArray<TVector3>(nx,ny,nz);
  for (int i=0;i<Efield->Length();i++)
    Efield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid to compute the impact form each cell on each cell.
  Epartial=new MultiArray<TVector3>(nx,ny,nz,nx,ny,nz);
  for (int i=0;i<Epartial->Length();i++)
    Epartial->GetFlat(i)->SetXYZ(0,0,0);

  //and to hold the external fieldmap
  Eexternal=new MultiArray<TVector3>(nx,ny,nz);
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,0);

  //ditto the external magnetic fieldmap
  Bfield=new MultiArray<TVector3>(nx,ny,nz);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid of charges in those volumes.
  q=new MultiArray<float>(nx,ny,nz);
  for (int i=0;i<q->Length();i++)
    *(q->GetFlat(i))=0;

  //define the step:
  step.SetXYZ(dim.X()/nx,dim.Y()/ny,dim.Z()/nz);

  //define the external fields:
  setFlatFields(1.4/*tesla*/,200/*V/cm*/);
  return;
}

TVector3 FieldSim::calc_unit_field(TVector3 at, TVector3 from){
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

TVector3 FieldSim::fieldIntegral(float zdest,TVector3 start){
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
      fieldInt+=Efield->Get(x,y,i)*step.Z();
  }
    fieldInt-=Efield->Get(x,y,zi)*(startz-zi*step.Z());//remove the part of the low end cell we didn't travel through
    fieldInt+=Efield->Get(x,y,zf)*(endz-zf*step.Z());//add the part of the high end cell we did travel through

    
  return dir*fieldInt;
}

void FieldSim::populate_fieldmap(){
//sum the E field at every point in the nx by ny by nz grid of field points.
//could also do this 'flat', but doesn't save any time
  //printf("in pop_fieldmap, n=(%d,%d,%d)\n",nx,ny,nz);
  for (int ix=0;ix<nx;ix++){
    for (int iy=0;iy<ny;iy++){
      for (int iz=0;iz<nz;iz++){
	//*f[ix][iy][iz]=sum_field_at(nx,ny,nz,partial_,q_,ix,iy,iz);
	TVector3 localF=sum_field_at(ix,iy,iz);
	Efield->Set(ix,iy,iz,localF);
	//if (localf.Mag()>1e-9)
	//printf("fieldmap@ (%d,%d,%d) mag=%f\n",ix,iy,iz,localF.Mag());
      }
    }
  }
  return;
}

void  FieldSim::populate_lookup(){
  //with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.

  //  TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  //printf("populating lookup for (%dx%dx%d)x(%dx%dx%d) grid\n",fx,fy,fz,ox,oy,oz);
  TVector3 at(0,0,0);
  TVector3 from(0,0,0);

  for (int ifx=0;ifx<nx;ifx++){
    at.SetX(step.X()*(ifx+0.5));
    for (int ify=0;ify<ny;ify++){
      at.SetY(step.Y()*(ify+0.5));
      for (int ifz=0;ifz<nz;ifz++){
	at.SetZ(step.Z()*(ifz+0.5));
	for (int iox=0;iox<nx;iox++){
	  from.SetX(step.X()*(iox+0.5));
	  for (int ioy=0;ioy<ny;ioy++){
	    from.SetY(step.Y()*(ioy+0.5));
	    for (int ioz=0;ioz<nz;ioz++){
	      from.SetZ(step.Z()*(ioz+0.5));
	      //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
	      //printf("calc_unit_field...\n");
	      Epartial->Set(ifx,ify,ifz,iox,ioy,ioz,calc_unit_field(at,from));
	    }
	  }
	}
      }
    }
  }
  return;

}

void FieldSim::setFlatFields(float B, float E){
  printf("FieldSim::setFlatFields(B=%f,E=%f)\n",B,E);
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,E);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,B);
  return;
}

TVector3 FieldSim::sum_field_at(int x,int y, int z){
 //sum the E field over all nx by ny by nz cells of sources, at the specific position x,y,z.
  TVector3 sum(0,0,0);
  for (int ix=0;ix<nx;ix++){
    for (int iy=0;iy<ny;iy++){
      for (int iz=0;iz<nz;iz++){
	//sum+=*partial[x][y][z][ix][iy][iz] * *q[ix][iy][iz];
	if (x==ix && y==iy && z==iz) continue;
	sum+=Epartial->Get(x,y,z,ix,iy,iz)*q->Get(ix,iy,iz);
      }
    }
  }
  sum+=Escale*Eexternal->Get(x,y,z);
  //printf("summed field at (%d,%d,%d)=(%f,%f,%f)\n",x,y,z,sum.X(),sum.Y(),sum.Z());
  return sum;
}

TVector3 FieldSim::swimTo(float zdest,TVector3 start){

 //using second order langevin expansion from http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-016.pdf
  //TVector3 (*field)[nx][ny][nz]=field_;

  //set the direction of the external fields.
  TVector3 B=Bfield->Get(0,0,0)*Bscale;//static field in tesla T=Vs/m^2
  
  //int x=start.X()/step.X();
  //int y=start.Y()/step.Y();
  //int zi=start.Z()/step.Z();
  double zdist=zdest-start.Z();


  TVector3 fieldInt=fieldIntegral(zdest,start);
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
  double deltaX=c0*fieldInt.X()/fieldz+c1*fieldInt.Y()/fieldz-c1*B.Y()/B.Z()*zdist+c2*B.X()/B.Z()*zdist;
  double deltaY=c0*fieldInt.Y()/fieldz-c1*fieldInt.X()/fieldz+c1*B.X()/B.Z()*zdist+c2*B.Y()/B.Z()*zdist;
  double deltaZ=0; //not correct, but small?  different E leads to different drift velocity.  see linked paper.  fix eventually.

  //printf("swimTo: (%2.4f,%2.4f,%2.4f) (cell %d,%d,%d) to z=%2.4f\n",start.X(),start.Y(), start.Z(),x,y,zi,zdest);
  printf("swimTo: fieldInt=(%2.4f,%2.4f,%2.4f)\n",fieldInt.X(),fieldInt.Y(),fieldInt.Z());
  
  TVector3 dest(start.X()+deltaX,start.Y()+deltaY,zdest+deltaZ);
  
  return dest;
  
}
