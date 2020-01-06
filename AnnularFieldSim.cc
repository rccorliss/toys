#include "TVector3.h"
#include "AnnularFieldSim.h"
#include "TH3F.h"

#define ALMOST_ZERO 0.00001

AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
				 int r, int roi_r0, int roi_r1, int in_rSpacing,
				 int phi, int roi_phi0, int roi_phi1, int in_phiSpacing,
				 int z, int roi_z0, int roi_z1,int in_zSpacing,
				 float vdr){
  //check well-ordering:
  if (roi_r0 >=r || roi_r1>r || roi_r0>=roi_r1){
    assert(1==2);
  }
  if (roi_phi0 >=phi || roi_phi1>phi || roi_phi0>=roi_phi1){
    printf("phi roi is out of range or spans the wrap-around.  Please spare me that math.\n");
    assert(1==2);
  }
  if (roi_z0 >=z || roi_z1>z || roi_z0>=roi_z1){
    assert(1==2);
  }

  printf("AnnularFieldSim::AnnularFieldSim with (%dx%dx%d) grid\n",r,phi,z);

  //load constants of motion, dimensions, etc:
  Escale=1; Bscale=1;
  vdrift=vdr;
  rmin=in_innerRadius; rmax=in_outerRadius;
  zmin=0;zmax=in_outerZ;

 //define the size of the volume:
  dim.SetXYZ(1,0,0);
  dim.SetPerp(rmax-rmin);
  dim.SetPhi(0);
  phispan=2*TMath::Pi();
  dim.SetZ(dz);
  

  //load parameters of the whole-volume tiling
  nr=r;nphi=phi;nz=z; //number of fundamental bins (f-bins) in each direction

  //calculate the size of an f-bin:
  //note that you have to set a non-zero value to start or perp won't update.
  step.SetXYZ(1,0,0);
  step.SetPerp(dim.Perp()/nr);
  step.SetPhi(phispan/nphi);
  step.SetZ(dim.Z()/nz);
  // printf("f-bin size:  r=%f,phi=%f, wanted %f,%f\n",step.Perp(),step.Phi(),dr/r,dphi/phi);

  
  //load parameters of our region of interest
  rmin_roi=roi_r0; phimin_roi=roi_phi0; zmin_roi=roi_z0; //lower edge of our region of interest, measured in f-bins
  rmax_roi=roi_r1; phimax_roi=roi_phi1; zmax_roi=roi_z1; //exlcuded upper edge of our region of interest, measured in f-bins

  //calculate the dimensions, in f-bins in our region of interest
  nr_roi=rmax_roi-rmin_roi;
 nphi_roi=phimax_roi-phimin_roi;
  nz_roi=zmax_roi-zmin_roi;

  //create a field grid for the roi with the specified dimensions
  Efield=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi);
  for (int i=0;i<Efield->Length();i++)
    Efield->GetFlat(i)->SetXYZ(0,0,0);
  

  //load parameters of our high-resolution behavior
 


  
  //todo:
  //1) check math to find the right spacings.  I think it's:
  //if we ask for a specific number of bins, we have to choose whether 


  /* old version where we derive the spacing from the required number of bins.
     This had a problem because if we give an nbins that doesn't have a multiple near our actual number, there would be a large leftover bin.
  r_spacing=0;
  while (r_spacing*(nr_low-1)<nr) r_spacing++;//keep incrementing until the start of the last bin doesn't fit in the region, then go back one:
  r_spacing--; //we either fit exactly, or the last bin is an odd shape.  But if we do this, we can't tell whether it's short or long.  I'd rather force it to be short.
  phi_spacing=0;
  while (phi_spacing*nphi_low<nphi) phi_spacing++;//keep incrementing until we are either exactly nr, or bigger, so our last bin is shorter than the rest.
  z_spacing=0;
  while (z_spacing*nz_low<nz) z_spacing++;//keep incrementing until we are either exactly nr, or bigger, so our last bin is shorter than the rest.
  */
  //new version:  user specifies the spacings instead.  We're now guaranteed the last bin is full size or short, which is easier to handle.
  r_spacing=in_r_spacing;
  phi_spacing=in_phi_spacing;
  z_spacing=in_z_spacing;
  nr_low=ceil(nr/(r_spacing*1.0));
  nphi_low=ceil(nphi/(phi_spacing*1.0));
  nz_low=ceil(nz/(z_spacing*1.0));
  

  //todo:
  //read the lookup on the relative highres-to-highres to get the local effects
  //then interpolate the lookup between the lowres-to-lowres on both ends to get the highest res to fill out to the next spacing
  //then interpolate the lookupbetween the lowres-to-lowres just on the destination end, summing over the whole lowres on the origin end
  //this should reduce to the current system if: highres neighborhood always encompasses the entire roi
  //or if: lowres spacing is 1 (meaning n*_low=n*)
  

  /* now divided into highres and lowres segments:
  //and the 'lookup' grid to compute the field in each cell given charge in each other cell.
  Epartial=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi,nr,nphi,nz);
  for (int i=0;i<Epartial->Length();i++)
    Epartial->GetFlat(i)->SetXYZ(0,0,0);
  */

  //create the low-res, absolute 'lookup' grid to compute the field in each cell given charge in each other cell.
  Epartial_lowres=new MultiArray<TVector3>(nr_low,nphi_low,nz_low,nr_low,nphi_low,nz_low);
  for (int i=0;i<Epartial_lowres->Length();i++)
    Epartial_lowres->GetFlat(i)->SetXYZ(0,0,0);

  //create the high-res, relative 'lookup' grid to compute the field in each cell given charge in each nearby cell.
  Epartial_highres=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi,nr_high,nphi_high,nz_high);
  for (int i=0;i<Epartial_highres->Length();i++)
    Epartial_highres->GetFlat(i)->SetXYZ(0,0,0);

  

  
  //and to hold the external fieldmap over the region of interest in high res
  Eexternal=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi);
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,0);

  //ditto the external magnetic fieldmap
  Bfield=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,0);

  //and a grid of charges in the entire source volume, regardless of roi.
  //TODO:  this could be changed to be a TH3F...
  //space charge in coulomb units.
  q=new MultiArray<float>(nr,nphi,nz);
  for (int i=0;i<q->Length();i++)
    *(q->GetFlat(i))=0;
  q_lowres=new MultiArray<float>(nr_low,nphi_low,nz_low);
  for (int i=0;i<q_lowres->Length();i++)
    *(q_lowres->GetFlat(i))=0;



  //define the external fields:
  setFlatFields(1.4/*tesla*/,200/*V/cm*/);
  return;
}

AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
				 int r, int roi_r0, int roi_r1, int in_rSpacing,
				 int phi, int roi_phi0, int roi_phi1, int in_phiSpacing,
				 int z, int roi_z0, int roi_z1,int in_zSpacing,
				 float vdr)
  :
  AnnularFieldSim(in_innerRadius, in_outerRadius, in_outerZ,
		  r, roi_r0, roi_r1, 1,
		  phi, roi_phi0, roi_phi1, 1,
		  z, roi_z0, roi_z1, 1,
		  vdr){
  {
    //just passing through for the old version again
    return;
  }
AnnularFieldSim::AnnularFieldSim(float rin, float rout, float dz,int r, int phi, int z, float vdr)
  :
  AnnularFieldSim( rin,  rout,  dz,
		   r,0, r,
		   phi, 0, phi,
		   z, 0, z,
		   vdr)
{
  //just a pass-through to go from the old style to the more detailed version.
  return;
}

TVector3 AnnularFieldSim::calc_unit_field(TVector3 at, TVector3 from){
  //note this is the field due to a fixed point charge in free space.
  //if doing cylindrical calcs with different boundary conditions, this needs to change.

  //this could check roi bounds before returning, if things start acting funny.
  
  const float k=8.987*1e13;//=1/(4*pi*eps0) in N*cm^2/C^2 in a vacuum. N*cm^2/C units, so that we supply space charge in coulomb units.
  TVector3 delr=at-from;
  TVector3 field=delr; //to set the direction.
  if (delr.Mag()<ALMOST_ZERO*ALMOST_ZERO){ //note that this has blurred units -- it should scale with all three dimensions of stepsize.  For lots of phi bins, especially, this might start to read as small before it's really small.
    //do nothing.  the vector is already zero, which will be our approximation.
    //field.SetMag(0);//no contribution if we're in the same cell. -- but root warns if trying to resize something of magnitude zero.
  } else{
    field.SetMag(k*1/(delr*delr));//scalar product on the bottom.
  }
  //printf("calc_unit_field at (%2.2f,%2.2f,%2.2f) from  (%2.2f,%2.2f,%2.2f).  Mag=%2.4fe-9\n",at.x(),at.Y(),at.Z(),from.X(),from.Y(),from.Z(),field.Mag()*1e9);
  
  return field;
}
int AnnularFieldSim::FilterPhiIndex(int phi,int range=nphi){
  //shifts phi up or down by nphi (=2pi in phi index space), and complains if this doesn't put it in range.
  int p=phi;
  if (p>=range){
    p-=range;
  }
  if (p<0){
    p+=range;
  }
  if (p>=range || p<0){
    printf("AnnularFieldSim::FilterPhiIndex asked to filter %d, which is more than range=%d out of bounds.  Check what called this.\n",phi,range);
  }
  return p;
}
  
AnnularFieldSim::BoundsCase AnnularFieldSim::GetRindexAndCheckBounds(float pos, int *r){
  float r0f=(pos-rmin)/step.Perp(); //the position in r, in units of step, starting from the low edge of the 0th bin.
  int r0=floor(r0f);
  *r=r0;

    int r0lowered_slightly=floor(r0f-ALMOST_ZERO);
  int r0raised_slightly=floor(r0f+ALMOST_ZERO); 
  if (r0lowered_slightly>=rmax_roi || r0raised_slightly<rmin_roi){
    return OutOfBounds;
  }

    //now if we are out of bounds, it is because we are on an edge, within range of ALMOST_ZERO of being in fair territory.
  if (r0>=rmax_roi){
    return OnHighEdge;
  }
  if (r0<rmin_roi){
    return OnLowEdge;
  }
  //if we're still here, we're in bounds.
  return InBounds;

}
AnnularFieldSim::BoundsCase AnnularFieldSim::GetPhiIndexAndCheckBounds(float pos, int *phi){
  float p0f=(pos)/step.Phi(); //the position in phi, in units of step, starting from the low edge of the 0th bin.
  int p0=FilterPhiIndex(floor(p0f));
  *phi=p0;

  int p0lowered_slightly=FilterPhiIndex(floor(p0f-ALMOST_ZERO));
  int p0raised_slightly=FilterPhiIndex(floor(p0f+ALMOST_ZERO));
  //annoying detail:  if we are at index 0, we might go above pmax by going down.
  // and if we are at nphi-1, we might go below pmin by going up.
  // if we are not at p0=0 or nphi-1, the slight wiggles won't move us.
  // if we are at p0=0, we are not at or above the max, and lowering slightly won't change that,
  // and is we are at p0=nphi-1, we are not below the min, and raising slightly won't change that
  if ((p0lowered_slightly>=phimax_roi && p0!=0)  || (p0raised_slightly<phimin_roi && p0!=(nphi-1))){
    return OutOfBounds;
  }
 //now if we are out of bounds, it is because we are on an edge, within range of ALMOST_ZERO of being in fair territory.
  if (p0>=phimax_roi){
    return OnHighEdge;
  }
  if (p0<phimin_roi){
    return OnLowEdge;
  }
  //if we're still here, we're in bounds.
  return InBounds;

}
AnnularFieldSim::BoundsCase AnnularFieldSim::GetZindexAndCheckBounds(float pos, int *z){
  float z0f=(pos-zmin)/step.Z(); //the position in z, in units of step, starting from the low edge of the 0th bin.
  int z0=floor(z0f);
  *z=z0;

  int z0lowered_slightly=floor(z0f-ALMOST_ZERO);
  int z0raised_slightly=floor(z0f+ALMOST_ZERO); 

  if (z0lowered_slightly>=zmax_roi || z0raised_slightly<zmin_roi){
    return OutOfBounds;
  }
  //now if we are out of bounds, it is because we are on an edge, within range of ALMOST_ZERO of being in fair territory.
  if (z0>=zmax_roi){
    return OnHighEdge;
  }
  if (z0<zmin_roi){
    return OnLowEdge;
  }
  //if we're still here, we're in bounds.
  return InBounds;
}
  


 


    

TVector3 AnnularFieldSim::fieldIntegral(float zdest,TVector3 start){
  //integrates E dz, from the starting point to the selected z position.  The path is assumed to be along z for each step, with adjustments to x and y accumulated after each step.
  // printf("AnnularFieldSim::fieldIntegral(x=%f,y=%f, z=%f) to z=%f\n",start.X(),start.Y(),start.Z(),zdest);

  int r, phi;
  bool rOkay=  (GetRindexAndCheckBounds(start.Perp(), &r) == InBounds);
  bool phiOkay=  (GetPhiIndexAndCheckBounds(start.Phi(), &phi) == InBounds);

  if (!rOkay || !phiOkay){
    printf("AnnularFieldSim::fieldIntegral asked to integrate along (r=%f,phi=%f), index=(%d,%d), which is out of bounds.  returning starting position.\n",start.Perp(),start.Phi(),r,phi);
    return start;
  }
  
  int dir=(start.Z()>zdest?-1:1);//+1 if going to larger z, -1 if going to smaller;  if they're the same, the sense doesn't matter.

  int zi, zf;
  double startz,endz;
  BoundsCase startBound,endBound;

  //make sure 'zi' is always the smaller of the two numbers, for handling the partial-steps.
  if (dir>0){
    startBound=GetZindexAndCheckBounds(start.Z(),&zi); //highest cell with lower bound less than lower bound of integral
    endBound=GetZindexAndCheckBounds(zdest,&zf); //highest cell with lower bound less than upper lower bound of integral
    startz=start.Z();
    endz=zdest;
  } else{
    endBound=GetZindexAndCheckBounds(start.Z(),&zf); //highest cell with lower bound less than lower bound of integral
    startBound=GetZindexAndCheckBounds(zdest,&zi); //highest cell with lower bound less than upper lower bound of integral
    startz=zdest;
    endz=start.Z();
  }
  bool startOkay=(startBound==InBounds); 
  bool endOkay=(endBound==InBounds || endBound==OnHighEdge); //if we just barely touch out-of-bounds on the high end, we can skip that piece of the integral
  
  if (!startOkay || !endOkay){
    printf("AnnularFieldSim::fieldIntegral asked to integrate from z=%f to %f, index=%d to %d), which is out of bounds.  returning starting position.\n",startz,endz,zi,zf);
    return start;
  }
 
  TVector3 fieldInt(0,0,0);
  // printf("AnnularFieldSim::fieldIntegral requesting (%d,%d,%d)-(%d,%d,%d) (inclusive) cells\n",r,phi,zi,r,phi,zf-1);
  for(int i=zi;i<zf;i++){ //count the whole cell of the lower end, and skip the whole cell of the high end.
    TVector3 tf=Efield->Get(r-rmin_roi,phi-phimin_roi,i-zmin_roi);
      //printf("fieldAt (%d,%d,%d)=(%f,%f,%f) step=%f\n",r,phi,i,tf.X(),tf.Y(),tf.Z(),step.Z());
      fieldInt+=tf*step.Z();
  }
  
  //since bins contain their lower bound, but not their upper, I can safely remove the unused portion of the lower cell:
    fieldInt-=Efield->Get(r-rmin_roi,phi-phimin_roi,zi-zmin_roi)*(startz-zi*step.Z());//remove the part of the low end cell we didn't travel through

//but only need to add the used portion of the upper cell if we go past the edge of it meaningfully:
    if (endz/step.Z()-zf>ALMOST_ZERO){
      //printf("endz/step.Z()=%f, zf=%f\n",endz/step.Z(),zf*1.0);
      //if our final step is actually in the next step.
      fieldInt+=Efield->Get(r-rmin_roi,phi-phimin_roi,zf-zmin_roi)*(endz-zf*step.Z());//add the part of the high end cell we did travel through
    }

    
  return dir*fieldInt;
}

TVector3 AnnularFieldSim::GetCellCenter(int r, int phi, int z){
  //returns the midpoint of the cell (halfway between each edge, not weighted center)
  
  TVector3 c(1,0,0);
  c.SetPerp((r+0.5)*step.Perp()+rmin);
  c.SetPhi((phi+0.5)*step.Phi());
  c.SetZ((z+0.5)*step.Z());

  return c;
}

TVector3 AnnularFieldSim::GetGroupCellCenter(int r0, int r1, int phi0, int phi1, int z0, int z1){
  //returns the midpoint of the cell (halfway between each edge, not weighted center)
  float ravg=(r0+r1)/2.0+0.5;
  if(phi0>phi1) phi1+=nphi;
  if(phi0>phi1) {
    printf("phi1(%d)<=phi0(%d) even after boosting phi1.  check what called this!\n",phi1,phi0);
    asser(1==2);
  }
  float phiavg=(r0+r1)/2.0+0.5;
  if (phiavg>=nphi) phiavg-=nphi;

  float zavg=(z0+z1)/2.0+0.5;

  
  TVector3 c(1,0,0);
  c.SetPerp((ravg)*step.Perp()+rmin);
  c.SetPhi((phiavg)*step.Phi());
  c.SetZ((zavg)*step.Z());

  return c;
}

TVector3 AnnularFieldSim::GetWeightedCellCenter(int r, int phi, int z){
  //returns the weighted center of the cell by volume.
  //todo:  this is vaguely hefty, and might be worth storing the result of, if speed is needed
  TVector3 c(1,0,0);

  float rin=r*step.Perp()+rmin;
  float rout=rin+step.Perp();

  float rMid=(4*TMath::Sin(step.Phi()/2)*(pow(rout,3)-pow(rin,3))
	      /(3*step.Phi()*(pow(rout,2)-pow(rin,2))));
  c.SetPerp(rMid);
  c.SetPhi((phi+0.5)*step.Phi());
  c.SetZ((z+0.5)*step.Z());

  return c;
}



TVector3 AnnularFieldSim::interpolatedFieldIntegral(float zdest,TVector3 start){
  printf("AnnularFieldSim::interpolatedFieldIntegral(x=%f,y=%f, z=%f)\n",start.X(),start.Y(),start.Z());


  float r0=(start.Perp()-rmin)/step.Perp()-0.5; //the position in r, in units of step, starting from the center of the 0th bin.
  int r0i=floor(r0); //the integer portion of the position. -- what center is below our position?
  float r0d=r0-r0i;//the decimal portion of the position. -- how far past center are we?
  int ri[4];//the r position of the four cell centers to consider
  ri[0]=ri[1]=r0i;
  ri[2]=ri[3]=r0i+1;
  float rw[4];//the weight of that cell, linear in distance from the center of it
  rw[0]=rw[1]=1-r0d; //1 when we're on it, 0 when we're at the other one.
  rw[2]=rw[3]=r0d; //1 when we're on it, 0 when we're at the other one

  bool skip[]={false,false,false,false};
  if (ri[0]<rmin_roi){
    skip[0]=true; //don't bother handling 0 and 1 in the coordinates.
    skip[1]=true;
    rw[2]=rw[3]=1; //and weight like we're dead-center on the outer cells.
  }
  if (ri[2]>=rmax_roi){
    skip[2]=true; //don't bother handling 2 and 3 in the coordinates.
    skip[3]=true;
    rw[0]=rw[1]=1; //and weight like we're dead-center on the inner cells.
  }
  
  //now repeat that structure for phi:
  float p0=(start.Phi())/step.Phi()-0.5; //the position in phi, in units of step, starting from the center of the 0th bin.
  int p0i=floor(p0); //the integer portion of the position. -- what center is below our position?
  float p0d=p0-p0i;//the decimal portion of the position. -- how far past center are we?
  int pi[4];//the phi position of the four cell centers to consider
  pi[0]=pi[2]=FilterPhiIndex(p0i);
  pi[1]=pi[3]=FilterPhiIndex(p0i+1);
  float pw[4];//the weight of that cell, linear in distance from the center of it
  pw[0]=pw[2]=1-p0d; //1 when we're on it, 0 when we're at the other one.
  pw[1]=pw[3]=p0d; //1 when we're on it, 0 when we're at the other one

  if (pi[0]<phimin_roi){
    skip[0]=true; //don't bother handling 0 and 1 in the coordinates.
    skip[2]=true;
    pw[1]=pw[3]=1; //and weight like we're dead-center on the outer cells.
  }
  if (pi[1]>=phimax_roi){
    skip[1]=true; //don't bother handling 2 and 3 in the coordinates.
    skip[3]=true;
    pw[0]=pw[2]=1; //and weight like we're dead-center on the inner cells.
  }
  

   //printf("interpolating fieldInt at  r=%f,phi=%f\n",r0,phi0);

  int dir=(start.Z()>zdest?-1:1);//+1 if going to larger z, -1 if going to smaller;  if they're the same, the sense doesn't matter.

  int zi, zf;
  double startz,endz;
  BoundsCase startBound,endBound;

  //make sure 'zi' is always the smaller of the two numbers, for handling the partial-steps.
  if (dir>0){
    startBound=GetZindexAndCheckBounds(start.Z(),&zi); //highest cell with lower bound less than lower bound of integral
    endBound=GetZindexAndCheckBounds(zdest,&zf); //highest cell with lower bound less than upper lower bound of integral
    startz=start.Z();
    endz=zdest;
  } else{
    endBound=GetZindexAndCheckBounds(start.Z(),&zf); //highest cell with lower bound less than lower bound of integral
    startBound=GetZindexAndCheckBounds(zdest,&zi); //highest cell with lower bound less than upper lower bound of integral
    startz=zdest;
    endz=start.Z();
  }
  bool startOkay=(startBound==InBounds); //maybe todo: add handling for being just below the low edge.
  bool endOkay=(endBound==InBounds || endBound==OnHighEdge); //if we just barely touch out-of-bounds on the high end, we can skip that piece of the integral
  
  if (!startOkay || !endOkay){
    printf("AnnularFieldSim::InterpolatedFieldIntegral asked to integrate from z=%f to %f, index=%d to %d), which is out of bounds.  returning starting position.\n",startz,endz,zi,zf);
    return start;
  }




  


  TVector3 fieldInt, partialInt;//where we'll store integrals as we generate them.
  
  for (int i=0;i<4;i++){
    if (skip[i]) continue; //we invalidated this one for some reason.
    partialInt.SetXYZ(0,0,0);
    printf("looking for element r=%d,phi=%d\n",ri[i],pi[i]);
    for(int j=zi;j<zf;j++){ //count the whole cell of the lower end, and skip the whole cell of the high end.
      
      partialInt+=Efield->Get(ri[i]-rmin_roi,pi[i]-phimin_roi,j-zmin_roi)*step.Z();
    }
  
    partialInt-=Efield->Get(ri[i]-rmin_roi,pi[i]-phimin_roi,zi-zmin_roi)*(startz-zi*step.Z());//remove the part of the low end cell we didn't travel through
    if (endz/step.Z()-zf>ALMOST_ZERO){
    partialInt+=Efield->Get(ri[i]-rmin_roi,pi[i]-phimin_roi,zf-zmin_roi)*(endz-zf*step.Z());//add the part of the high end cell we did travel through
    }
    fieldInt+=rw[i]*pw[i]*partialInt;
  }
    
  return dir*fieldInt;
}

void AnnularFieldSim::load_spacecharge(TH3F *hist, float zoffset, float scalefactor=1){
  //load spacecharge densities from a histogram, where scalefactor translates into local units
  //noting that the histogram limits may differ from the simulation size, and have different granularity
  //hist is assumed/required to be x=phi, y=r, z=z
  //z offset 'drifts' the charge by that distance toward z=0.

  //Get dimensions of input
  float hrmin=hist->GetYaxis()->GetXmin();
  float hrmax=hist->GetYaxis()->GetXmax();
  float hphimin=hist->GetXaxis()->GetXmin();
  float hphimax=hist->GetXaxis()->GetXmax();
  float hzmin=hist->GetZaxis()->GetXmin();
  float hzmax=hist->GetZaxis()->GetXmax();
  
  //Get number of bins in each dimension
  int hrn=hist->GetNbinsY();
  int hphin=hist->GetNbinsX();
  int hzn=hist->GetNbinsZ();

  //do some computation of steps:
  float hrstep=(hrmax-hrmin)/hrn;
  float hphistep=(hphimax-hphimin)/hphin;
  float hzstep=(hzmax-hzmin)/hzn;

  //clear the previous spacecharge dist:
    for (int i=0;i<q->Length();i++)
    *(q->GetFlat(i))=0;

  
  //loop over every bin and add that to the internal model:
  //note that bin 0 is the underflow, so we need the +1 internally

  //the minimum r we need is localr=0, hence:
  //hr=localr*dr+rmin
  //localr*dr+rmin-hrmin=hrstep*(i+0.5)
  //i=(localr*dr+rmin-hrmin)/hrstep

  for (int i=(rmin-hrmin)/hrstep;i<hrn;i++){
    float hr=hrmin+hrstep*(i+0.5);//bin center
    int localr=(hr-rmin)/step.Perp();
    if (localr<0){
      printf("Loading from histogram has r out of range! r=%f < rmin=%f\n",hr,rmin);      
      continue;
    }
    if (localr>=nr){
      printf("Loading from histogram has r out of range! r=%f > rmax=%f\n",hr,rmax);
      i=hrn; //no reason to keep looking at larger r.
      continue;
    }
    for (int j=0;j<hphin;j++){
      float hphi=hphimin+hphistep*(j+0.5); //bin center
      int localphi=hphi/step.Phi();
      if (localphi>=nphi){//handle wrap-around:
	localphi-=nphi;
      }
      if (localphi<0){//handle wrap-around:
	localphi+=nphi;
      }
      //todo:  should add ability to take in a phi- slice only
      for (int k=(zmin-(hzmin-zoffset))/hzstep;k<hzn;k++){
	float hz=hzmin-zoffset+hzstep*(k+0.5);//bin center
	int localz=hz/step.Z();
	if (localz<0){
	  printf("Loading from histogram has z out of range! z=%f < zmin=%f\n",hz,zmin);
	  continue;
	}
	if (localz>=nz){
	  printf("Loading from histogram has z out of range! z=%f > zmax=%f\n",hz,zmax);
	  k=hzn;//no reason to keep looking at larger z.
	  continue;
	}
	//can be simplified:  float vol=hzstep*(hphistep*(hr+hrstep)*(hr+hrstep) - hphistep*hr*hr);
	float vol=hzstep*hphistep*(2*hr+hrstep)*hrstep;
	float qbin=scalefactor*vol*hist->GetBinContent(hist->GetBin(j+1,i+1,k+1));
	//float qold=q->Get(localr,localphi,localz);
	//printf("loading Q=%f from hist(%d,%d,%d) into cell (%d,%d,%d), qold=%f\n",qbin,i,j,k,localr,localphi,localz,qold);
	q->Add(localr,localphi,localz,qbin);
      }
    }
  }

  //todo:  go through the q and build q_lowres.
  for (int i=0;i<q_lowres->Length();i++)
    *(q_lowres->GetFlat(i))=0;

  //fill our low-res
  //note that this assumes the last bin is short or normal length, not long.
  for (int ifr=0;ifr<nr;ifr++){
    int r_low=ifr/r_spacing;
    for (int ifphi=0;ifphi<nphi;ifphi++){
      int phi_low=ifphi/phi_spacing;
      if (phi_high>=nphi) phi_high=nphi-1;
      for (int ifz=0;ifz<nz;ifz++){
	int z_low=ifz/z_spacing;
  	q_lowres->Add(r_low,phi_low,z_low,q->Get(ifr,ifphi,ifz));
      }
    }
  }


  return;
}

void AnnularFieldSim::populate_fieldmap(){
  //sum the E field at every point in the region of interest
  // remember that Efield uses relative indices
  //printf("in pop_fieldmap, n=(%d,%d,%d)\n",nr,ny,nz);
  for (int ir=rmin_roi;ir<rmax_roi;ir++){
    for (int iphi=phimin_roi;iphi<phimax_roi;iphi++){
      for (int iz=zmin_roi;iz<zmax_roi;iz++){
	TVector3 localF=sum_field_at(ir,iphi,iz);
	Efield->Set(ir-rmin_roi,iphi-phimin_roi,iz-zmin_roi,localF);
	//if (localF.Mag()>1e-9)
	//printf("fieldmap@ (%d,%d,%d) mag=%f\n",ir,iphi,iz,localF.Mag());
      }
    }
  }
  return;
}
  
void  AnnularFieldSim::populate_lookup(){
  //with 'f' being the position the field is being measured at, and 'o' being the position of the charge generating the field.
  //remember the 'f' part of Epartial uses relative indices.
  //  TVector3 (*f)[fx][fy][fz][ox][oy][oz]=field_;
  //printf("populating lookup for (%dx%dx%d)x(%dx%dx%d) grid\n",fx,fy,fz,ox,oy,oz);


  populate_highres_lookup();
  populate_lowres_lookup();
  return;
}

void AnnularFieldSim::populate_highres_lookup(){

  TVector3 at(1,0,0);
  TVector3 from(1,0,0);

  //populate_highres_lookup();
  int r_highres_dist=(nr_highres-1)/2;
  int phi_highres_dist=(nphi_highres-1)/2;
  int z_highres_dist=(nz_highres-1)/2;



  int ncellsin[3][3][3];
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
	ncellsin[i][j][k]=0; //we could count total volume, but wihtout knowing the charge prior, it's not clear that'd be /better/
      }
    }
  }

  //todo: if this runs too slowly, I can do geometry instead of looping over all the cells that are possibly in range
  
  for (int ifr=rmin_roi;ifr<rmax_roi;ifr++){
    for (int ifphi=phimin_roi;ifphi<phimax_roi;ifphi++){
      for (int ifz=zmin_roi;ifz<zmax_roi;ifz++){
	at=GetCellCenter(ifr, ifphi, ifz);
	int r_parentlow=floor((ifr-r_highres_dist)/(r_spacing*1.0));//partly enclosed in our high-res region
	int r_parenthigh=floor((ifr+r_highres_dist)/(r_spacing*1.0))+1;//definitely not enclosed in our high-res region
	int phi_parentlow=floor((ifphi-phi_highres_dist)/(phi_spacing*1.0));
	int phi_parenthigh=floor((ifphi+phi_highres_dist)/(phi_spacing*1.0))+1; //note that this can be bigger than nphi!  We keep track of that.
	int z_parentlow=floor((ifz-z_highres_dist)/(z_spacing*1.0));
	int z_parenthigh=floor((ifz+z_highres_dist)/(z_spacing*1.0))+1;
	
	for (int ir=r_parentlow*r_spacing;ir<r_parenthigh*rspacing;ir++){
	  if (ir<0) ir=0;
	  if (ir>=nr) continue;
	  int rbin=(ir-r)+r_highres_dist;//zeroth bin when we're at max distance below, etc.
	  int rcell=1;
	  if (rbin<=0) {
	    rbin=0;
	    rcell=0;
	  }
	  if (rbin>=nr_highres){
	    rbin=nr_highres-1;
	    rcell=2;
	  }
	  for (int iphi=phi_parentlow*phi_spacing;iphi<phi_parenthigh*phispacing;iphi++){
	    //no phi checks since it's circular.
	    int phibin=(iphi-phi)+phi_highres_dist;
	    int phicell=1;
	    if (phibin<=0) {
	      phibin=0;
	      phicell=0;
	    }
	    if (phibin>=nphi_highres){
	      phibin=nphi_highres-1;
	      phicell=2;
	    }
	    for (int iz=z_parentlow*z_spacing;iz<z_parenthigh*zspacing;iz++){
	      if (iz<0) iz=0;
	      if (iz>=nz) continue;
	      int zbin=(iz-z)+z_highres_dist;
	      int zcell=1;
	      if (zbin<=0) {
		zbin=0;
		zcell=0;
	      }
	      if (zbin>=nz_highres){
		zbin=nz_highres-1;
		zcell=2;
	      }
	      from=GetCellCenter(rbin, phibin, zbin);
	      ncellsin[rcell][phicell][zcell]++;
	      int nc= ncellsin[rcell][phicell][zcell];
	      
	      TVector3 currentf=Epartial->Get(ifr-rmin_roi,ifphi-phimin_roi,ifz-zmin_roi,rbin,phibin,zbin);
	      //to keep this as the average, we multiply what's there back to its initial summed-but-not-divided value
	      //then add our new value, and the divide the new sum by the total number of cells
	      TVector3 newf=(currentf*(nc-1)+calc_unit_field(at,from))/nc;
	      Epartial_highres->Set(ifr-rmin_roi,ifphi-phimin_roi,ifz-zmin_roi,rbin,phibin,zbin,newf)
	    }
	  }
	}
      }
    }
  }
}

void AnnularFieldSim::populate_lowres_lookup(){

  TVector3 at(1,0,0);
  TVector3 from(1,0,0);

  int r_low,r_high,phi_low,phi_high,z_low,z_high;
  //todo: add checks to see if we're in the region of interest?  Otherwise, I can't scale this down to match the old.
  for (int ifr=0;ifr<nr_low;ifr++){
    r_low=ifr*r_spacing;
    r_high=r_low+r_spacing-1;
    if (r_high>=nr) r_high=nr-1;	
    for (int ifphi=0;ifphi<nphi_low;ifphi++){
      phi_low=ifphi*phi_spacing;
      phi_high=phi_low+phi_spacing-1;
      if (phi_high>=nphi) phi_high=nphi-1;
      for (int ifz=0;ifz<nz_low;ifz++){
	z_low=ifz*z_spacing;
	z_high=z_low+z_spacing-1;
	if (z_high>=nz) z_high=nz-1;
	at=GetGroupCellCenter(r_low,r_high,phi_low,phi_high,z_low,z_high);
	for (int ior=0;ior<nr_low;ior++){
	  r_low=ior*r_spacing;
	  r_high=r_low+r_spacing-1;
	  if (r_high>=nr) r_high=nr-1;
	  for (int iophi=0;iophi<nphi_low;iophi++){
	    phi_low=iophi*phi_spacing;
	    phi_high=phi_low+phi_spacing-1;
	    if (phi_high>=nphi) phi_high=nphi-1;
	    for (int ioz=0;ioz<nz_low;ioz++){
	      z_low=ioz*z_spacing;
	      z_high=z_low+z_spacing-1;
	      if (z_high>=nz) z_high=nz-1;
	      from=GetGroupCellCenter(r_low,r_high,phi_low,phi_high,z_low,z_high);

	      //*f[ifx][ify][ifz][iox][ioy][ioz]=cacl_unit_field(at,from);
	      //printf("calc_unit_field...\n");
	      Epartial_lowres->Set(ifr,ifphi,ifz,ior,iophi,ioz,calc_unit_field(at,from));
	    }
	  }
	}
      }
    }
  }
  return;

}

void AnnularFieldSim::setFlatFields(float B, float E){
  //these only cover the roi, but since we address them flat, we don't need to know that here.
  printf("AnnularFieldSim::setFlatFields(B=%f,E=%f)\n",B,E);
  printf("lengths:  Eext=%d, Bfie=%d\n",Eexternal->Length(),Bfield->Length());
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,E);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,B);
  return;
}

/*
Before the mixed-res revolution:
TVector3 AnnularFieldSim::sum_field_at(int r,int phi, int z){
 //sum the E field over all nr by ny by nz cells of sources, at the specific position r,phi,z.
  //note the specific position in Epartial is in relative coordinates.
  //printf("AnnularFieldSim::sum_field_at(r=%d,phi=%d, z=%d)\n",r,phi,z);
  TVector3 sum(0,0,0);
  for (int ir=0;ir<nr;ir++){
    for (int iphi=0;iphi<nphi;iphi++){
      for (int iz=0;iz<nz;iz++){
	//sum+=*partial[x][phi][z][ix][iphi][iz] * *q[ix][iphi][iz];
	if (r==ir && phi==iphi && z==iz) continue;
	sum+=Epartial->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi,ir,iphi,iz)*q->Get(ir,iphi,iz);
      }
    }
  }
  sum+=Escale*Eexternal->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi);
  //printf("summed field at (%d,%d,%d)=(%f,%f,%f)\n",x,y,z,sum.X(),sum.Y(),sum.Z());
  return sum;
}
*/

TVector3 AnnularFieldSim::sum_field_at(int r,int phi, int z){
 //sum the E field over all nr by ny by nz cells of sources, at the specific position r,phi,z.
  //note the specific position in Epartial is in relative coordinates.
  printf("AnnularFieldSim::sum_field_at(r=%d,phi=%d, z=%d)\n",r,phi,z);

  
  //reminder:
  //Epartial_lowres=new MultiArray<TVector3>(nr_low,nphi_low,nz_low,nr_low,nphi_low,nz_low);
  //Epartial_highres=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi,nr_high,nphi_high,nz_high);
  int r_highres_dist=(nr_highres-1)/2;
  int phi_highres_dist=(nphi_highres-1)/2;
  int z_highres_dist=(nz_highres-1)/2;

  //first, let's calculate our highres point's position in the lowres grid
  //this will be needed to do 
  

  //do the summation of the inner high-resolution region charges:
  //
  // bin 0  1 2 ...  n-2 n-1
  //  . . .|.|.|.|.|.|.|. . .
  //       | | | | | | |
  //  . . .|.|.|.|.|.|.|. . .
  //  -----+-+-+-+-+-+-+-----
  //  . . .|.|.|.|.|.|.|. . .
  //  -----+-+-+-+-+-+-+-----
  //  . . .|.|.|.|.|.|.|. . .
  //  -----+-+-+-+-+-+-+-----
  //  . . .|.|.|.|.|.|.|. . .
  //  -----+-+-+-+-+-+-+-----
  //  . . .|.|.|.|.|.|.|. . .
  //       | | | | | | |
  //  . . .|.|.|.|.|.|.|. . .
  //
  //

 
  int r_parentlow=floor((r-r_highres_dist)/(r_spacing*1.0));//partly enclosed in our high-res region
  int r_parenthigh=floor((r+r_highres_dist)/(r_spacing*1.0))+1;//definitely not enclosed in our high-res region
  int phi_parentlow=floor((phi-phi_highres_dist)/(phi_spacing*1.0));
  int phi_parenthigh=floor((phi+phi_highres_dist)/(phi_spacing*1.0))+1; //note that this can be bigger than nphi!  We keep track of that.
  int z_parentlow=floor((z-z_highres_dist)/(z_spacing*1.0));
  int z_parenthigh=floor((z+z_highres_dist)/(z_spacing*1.0))+1;

  //zero our current qlocal holder:
  for (int i=0;i<q_local->Length();i++)
    q_local->GetFlat(i)->SetXYZ(0,0,0);

  for (int ir=r_parentlow*r_spacing;ir<r_parenthigh*rspacing;ir++){
    if (ir<0) ir=0;
    if (ir>=nr) continue;
    int rbin=(ir-r)+r_highres_dist;//zeroth bin when we're at max distance below, etc.
    if (rbin<0) rbin=0;
    if (rbin>=nr_highres) rbin=nr_highres-1;
    for (int iphi=phi_parentlow*phi_spacing;iphi<phi_parenthigh*phispacing;iphi++){
      //no phi checks since it's circular.
      int phibin=(iphi-phi)+phi_highres_dist;
      if (phibin<0) phibin=0;
      if (phibin>=nphi_highres) phibin=nphi_highres-1;
      for (int iz=z_parentlow*z_spacing;iz<z_parenthigh*zspacing;iz++){
	if (iz<0) iz=0;
	if (iz>=nz) continue;
	int zbin=(iz-z)+z_highres_dist;
	if (zbin<0) zbin=0;
	if (zbin>=nz_highres) zbin=nz_highres-1;
	q_local->Add(rbin,phibin,zbin,q->Get(ir,FilterPhiIndex(iphi),iz));
      }
    }
  }

  //now q_highres is up to date for our current cell of interest.
  //start building our full sum by scaling the local lookup table by q.
  //note that the lookup table needs to have already accounted for cell centers.

  TVector3 sum(0,0,0);

  for (int ir=0;ir<nr_highres;ir++){
    for (int iphi=0;iphi<nphi_highres;iphi++){
      for (int iz=0;iz<nz_highres;iz++){
	sum+=Epartial_highres->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi,ir,iphi,iz)*q_local->Get(ir,iphi,iz);
      }
    }
  }

  //now we find our interpolated position between the eight nearby cells:
  //this is similar to the interpolated integral stuff, except we're at integer steps inside integer blocks
  //and we have z as well, now.
  bool skip[]={false,false,false,false,false,false,false,false};

  
  float r0=r/(1.0*r_spacing)-0.5; //the position in r, in units of r_spacing, starting from the center of the 0th bin.
  int r0i=floor(r0); //the integer portion of the position. -- what center is below our position?
  float r0d=r0-r0i;//the decimal portion of the position. -- how far past center are we?
  //instead of listing all eight, I'll address these as i%2, (i/2)%2 and (i/4)%2 to avoid typos
  int ri[2];//the r position of the eight cell centers to consider.
  ri[0]=r0i;
  ri[1]=r0i+1;
  float rw[2];//the weight of that cell, linear in distance from the center of it
  rw[0]=1-r0d; //1 when we're on it, 0 when we're at the other one.
  rw[1]=r0d; //1 when we're on it, 0 when we're at the other one


  
  if (ri[0]<rmin_roi/r_spacing){
    for (int i=0;i<8;i++)
      if ((i/4)%2==0)
	skip[i]=true; // don't bother handling ri[0].
    rw[1]=1; //and weight like we're dead-center on the outer cells.
  }
  if (ri[1]>=rmax_roi/r_spacing){
    for (int i=0;i<8;i++)
      if ((i/4)%2==1)
	skip[i]=true; // don't bother handling ri[1].
    rw[0]=1; //and weight like we're dead-center on the inner cells.
  }
  
  //now repeat that structure for phi:
  float p0=phi/(1.0*phi_spacing)-0.5; //the position in phi, in units of phi_spacing, starting from the center of the 0th bin.
  int p0i=floor(p0); //the integer portion of the position. -- what center is below our position?
  float p0d=p0-p0i;//the decimal portion of the position. -- how far past center are we?
  int pi[4];//the phi position of the four cell centers to consider
  pi[0]=FilterPhiIndex(p0i,nphi_low);
  pi[1]=FilterPhiIndex(p0i+1,nphi_low);
  float pw[2];//the weight of that cell, linear in distance from the center of it
  pw[0]=1-p0d; //1 when we're on it, 0 when we're at the other one.
  pw[1]=p0d; //1 when we're on it, 0 when we're at the other one

  if (pi[0]<phimin_roi/phi_spacing){
   for (int i=0;i<8;i++)
      if ((i/2)%2==0)
	skip[i]=true; // don't bother handling phii[1].
    pw[1]=1; //and weight like we're dead-center on the high cells. 

  }
  if (pi[1]>=phimax_roi/phi_spacing){
    for (int i=0;i<8;i++)
      if ((i/2)%2==1)
	skip[i]=true; // don't bother handling phii[1].
    pw[0]=1; //and weight like we're dead-center on the high cells. 
 }

  //and once more for z.  ooph.
  
  float z0=z/(1.0*z_spacing)-0.5; //the position in r, in units of r_spacing, starting from the center of the 0th bin.
  int z0i=floor(z0); //the integer portion of the position. -- what center is below our position?
  float z0d=z0-z0i;//the decimal portion of the position. -- how far past center are we?
  //instead of listing all eight, I'll address these as i%2, (i/2)%2 and (i/4)%2 to avoid typos
  int zi[2];//the r position of the eight cell centers to consider.
  zi[0]=z0i;
  zi[1]=z0i+1;
  float zw[2];//the weight of that cell, linear in distance from the center of it
  zw[0]=1-z0d; //1 when we're on it, 0 when we're at the other one.
  zw[1]=z0d; //1 when we're on it, 0 when we're at the other one


  
  if (zi[0]<zmin_roi/z_spacing){
    for (int i=0;i<8;i++)
      if ((i)%2==0)
	skip[i]=true; // don't bother handling zi[0].
    zw[1]=1; //and weight like we're dead-center on the higher cells.
  }
  if (zi[1]>=zmax_roi/z_spacing){
    for (int i=0;i<8;i++)
      if ((i)%2==1)
	skip[i]=true; // don't bother handling zi[1].
    zw[0]=1; //and weight like we're dead-center on the lower cells.
  }
  
  
  for (int ir=0;ir<nr_low;ir++){
    for (int iphi=0;iphi<nphi_low;iphi++){
      for (int iz=0;iz<nz_low;iz++){
	//conceptually: see if any full-res cell within the low-res cell is inside the high-res region:
	//the high-res region includes all indices from r-(nr_high-1)/2 to r+(nr_high-1)/2.
	//each low-res region includes all indices from ir*r_spacing to (ir+1)*r_spacing-1.
	int lrEdge={ir*r_spacing,(ir+1)*r_spacing-1};
	int hrEdge={r-r_highres_dist,r+r_highres_dist};
	if(lrEdge[0]<=hrEdge[1] && hrEdge[0]<=lrEdge[1]){
	  //if their bounds are interleaved, there is overlap, and we've already summed this region.
	  continue;
	} else {
	  for (int i=0;i<8;i++){
	    if (skip) continue;
	    sum+=(Epartial_lowres->Get(ri[(i/4)%2],pi[(i/2)%2],zi[(i)%2],ir,iphi,iz)
		  *q_lowres->Get(ir,iphi,iz))*zw[(i)%2]*pw[(i/2)%2]*rw[(i/4)%2];
	  }
	}
      }
    }
  }
  sum+=Escale*Eexternal->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi);
  //printf("summed field at (%d,%d,%d)=(%f,%f,%f)\n",x,y,z,sum.X(),sum.Y(),sum.Z());
  return sum;
}


TVector3 AnnularFieldSim::swimToInSteps(float zdest,TVector3 start,int steps=1, bool interpolate=false){
  //short-circuit if we're out of range:
  
  double zdist=zdest-start.Z();
  double zstep=zdist/steps;
  
  TVector3 ret=start;
  for (int i=0;i<steps;i++){
    int rt,pt,zt; //just placeholders
    BoundsCase zBound=GetZindexAndCheckBounds(ret.Z(),&zt);
    if (GetRindexAndCheckBounds(ret.Perp(),&rt)!=InBounds
	|| GetPhiIndexAndCheckBounds(ret.Phi(),&pt)!=InBounds
	|| (zBound!=InBounds && zBound!=OnHighEdge)){
      printf("AnnularFieldSim::swimToInSteps at step %d, asked to swim particle from (%f,%f,%f) which is outside the ROI.  Returning original position.\n",i,ret.X(),ret.Y(),ret.Z());
      return ret;
    }
    ret=swimTo(start.Z()+zstep*(i+1),ret,false);
  }
  
  return ret;
}

TVector3 AnnularFieldSim::swimTo(float zdest,TVector3 start, bool interpolate=false){

 //using second order langevin expansion from http://skipper.physics.sunysb.edu/~prakhar/tpc/Papers/ALICE-INT-2010-016.pdf
  //TVector3 (*field)[nr][ny][nz]=field_;
  int rt,pt,zt; //just placeholders
  BoundsCase zBound=GetZindexAndCheckBounds(start.Z(),&zt);
  if (GetRindexAndCheckBounds(start.Perp(),&rt)!=InBounds
      || GetPhiIndexAndCheckBounds(start.Phi(),&pt)!=InBounds
      || (zBound!=InBounds && zBound!=OnHighEdge)){
    printf("AnnularFieldSim::swimTo asked to swim particle from (%f,%f,%f) which is outside the ROI.  Returning original position.\n",start.X(),start.Y(),start.Z());
    return start;
  }
  
  //set the direction of the external fields.
  //todo: get this from a field map
  TVector3 B=Bfield->Get(0,0,0)*Bscale;//static field in tesla T=Vs/m^2
  
  //int x=start.X()/step.X();
  //int y=start.Y()/step.Y();
  //int zi=start.Z()/step.Z();
  double zdist=zdest-start.Z();

  //short-circuit if there's no travel length:
  if (TMath::Abs(zdist)<ALMOST_ZERO*step.Z()){
    printf("Asked to swim particle from (%f,%f,%f) to z=%f, which is a distance of %fcm.  Returning original position.\n",start.X(),start.Y(),start.Z(), zdest,zdist);
    return start;
  }

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

  //printf("swimTo: (%2.4f,%2.4f,%2.4f) to z=%2.4f\n",start.X(),start.Y(), start.Z(),zdest);
  //printf("swimTo: fieldInt=(%2.4f,%2.4f,%2.4f)\n",fieldInt.X(),fieldInt.Y(),fieldInt.Z());
  
  TVector3 dest(start.X()+deltaX,start.Y()+deltaY,zdest+deltaZ);
  
  return dest;
  
}


float AnnularFieldSim::RosseggerEterm(int m, int n, TVector3 at, TVector3 from){
  //from eqn 5.68 in Rossegger thesis, p 111.
  //cylindrical cavity from r=a to r=b, z=0 to z=L
  /*
  double L;
  double betaMN; // comes from solution to 5.12.
  double Esubz;
  int deltaM0;
  double eps0;
  double phi1;
  double phi;
  double r1;
  double r;
  double z1;
  double z;
  double dz;
  double dr;
  double dphi;
  
  
  */
  
  return 0;
}
