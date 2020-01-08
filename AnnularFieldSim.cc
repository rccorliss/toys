#include "TVector3.h"
#include "AnnularFieldSim.h"
#include "TH3F.h"

#define ALMOST_ZERO 0.00001

AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
				 int r, int roi_r0, int roi_r1, int in_rLowSpacing, int in_rHighSize,
				 int phi, int roi_phi0, int roi_phi1, int in_phiLowSpacing, int in_phiHighSize,
				 int z, int roi_z0, int roi_z1,int in_zLowSpacing, int in_zHighSize,
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

  //debug defaults:
  //
  debug_printActionEveryN=-1;
  debug_printCounter=0;
  //internal states used for the debug flag
  //goes with: if(debugFlag()) printf("%d: blah\n",__LINE__);
  
  
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
  dim.SetZ(in_outerZ);
  

  //load parameters of the whole-volume tiling
  nr=r;nphi=phi;nz=z; //number of fundamental bins (f-bins) in each direction
  printf("AnnularFieldSim::AnnularFieldSim set variables nr=%d, nphi=%d, nz=%d\n",nr,nphi,nz);

  //calculate the size of an f-bin:
  //note that you have to set a non-zero value to start or perp won't update.
  step.SetXYZ(1,0,0);
  step.SetPerp(dim.Perp()/nr);
  step.SetPhi(phispan/nphi);
  step.SetZ(dim.Z()/nz);
  // printf("f-bin size:  r=%f,phi=%f, wanted %f,%f\n",step.Perp(),step.Phi(),dr/r,dphi/phi);

  //create an array to store the charge in each f-bin
  q=new MultiArray<float>(nr,nphi,nz);
  for (int i=0;i<q->Length();i++)
    *(q->GetFlat(i))=0;
  
  

  
  //load parameters of our region of interest
  rmin_roi=roi_r0; phimin_roi=roi_phi0; zmin_roi=roi_z0; //lower edge of our region of interest, measured in f-bins
  rmax_roi=roi_r1; phimax_roi=roi_phi1; zmax_roi=roi_z1; //exlcuded upper edge of our region of interest, measured in f-bins
  printf("AnnularFieldSim::AnnularFieldSim set roi variables rmin=%d phimin=%d zmin=%d rmax=%d phimax=%d zmax=%d\n",
	 rmin_roi, phimin_roi, zmin_roi, rmax_roi, phimax_roi, zmax_roi);
  //calculate the dimensions, in f-bins in our region of interest
  nr_roi=rmax_roi-rmin_roi;
  nphi_roi=phimax_roi-phimin_roi;
  nz_roi=zmax_roi-zmin_roi;
  printf("AnnularFieldSim::AnnularFieldSim calc'd roi variables nr=%d nphi=%d nz=%d\n",nr_roi,nphi_roi,nz_roi);
  

  //create an array to hold the SC-induced electric field in the roi with the specified dimensions
  Efield=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi);
  for (int i=0;i<Efield->Length();i++)
    Efield->GetFlat(i)->SetXYZ(0,0,0);

  //and to hold the external electric fieldmap over the region of interest
  Eexternal=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi);
  for (int i=0;i<Eexternal->Length();i++)
    Eexternal->GetFlat(i)->SetXYZ(0,0,0);

  //ditto the external magnetic fieldmap
  Bfield=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi);
  for (int i=0;i<Bfield->Length();i++)
    Bfield->GetFlat(i)->SetXYZ(0,0,0);

  
  

  //load parameters of our high-resolution behavior
  nr_high=in_rHighSize;nphi_high=in_phiHighSize;nz_high=in_zHighSize; //dimensions, in f-bins of neighborhood of a f-bin in which we calculated the field in full resolution
  printf("AnnularFieldSim::AnnularFieldSim set highres variables nr=%d nphi=%d nz=%d \n",nr_high,nphi_high,nz_high);
	   

  //create the high-res, relative 'lookup' grid to compute the field in each f-bin in the roi given charge in each other nearby f-bin.
  //note that first and less elements in each direction average over multiple f-bins to match the l-bin spacing.
  printf("AnnularFieldSim::AnnularFieldSim building Epartial_highres with  nr_roi=%d nphi_roi=%d nz_roi=%d nr_high=%d nphi_high=%d nz_high=%d =~%2.2fM TVector3 objects\n",nr_roi,nphi_roi,nz_roi,nr_high,nphi_high,nz_high,
	 nr_roi*nphi_roi*nz_roi*nr_high*nphi_high*nz_high/(1.0e6));

  Epartial_highres=new MultiArray<TVector3>(nr_roi,nphi_roi,nz_roi,nr_high,nphi_high,nz_high);
  for (int i=0;i<Epartial_highres->Length();i++)
    Epartial_highres->GetFlat(i)->SetXYZ(0,0,0);

  printf("AnnularFieldSim::AnnularFieldSim building q_local nr_high=%d nphi_high=%d nz_high=%d =~%2.2fM floats\n",nr_high,nphi_high,nz_high,nr_high*nphi_high*nz_high/(1.0e6));
  //create a local charge grid that we will fill with the f-bin contents and match to Epartial:
  q_local=new MultiArray<float>(nr_high,nphi_high,nz_high);
  for (int i=0;i<q_local->Length();i++)
    *(q_local->GetFlat(i))=0;



  //load parameters of our low-resolution behavior
  r_spacing=in_rLowSpacing; phi_spacing=in_phiLowSpacing; z_spacing=in_zLowSpacing; //number of f-bins, in each direction, to gang together to make a single low-resolution bin (l-bin)

  //compute the number of l-bins in each direction.  In each dimension, the last bin may have a smaller number of f-bins than the rest, to match up with the total number of f-bins.
  //note that to properly handle rbins near the edge of the roi, our low-res roi needs to be larger by one bin in each direction.
  nr_low=ceil(nr/(r_spacing*1.0));
  nphi_low=ceil(nphi/(phi_spacing*1.0));
  nz_low=ceil(nz/(z_spacing*1.0));

  //calculate the lowest (inclusive) and highest (exclusive) l-bins that are in our region of interest
  rmin_roi_low=roi_r0/r_spacing;
  phimin_roi_low=roi_phi0/phi_spacing;
  zmin_roi_low=roi_z0/z_spacing; //lower edge of our region of interest, measured in l-bins
  rmax_roi_low=(roi_r1-1)/r_spacing+1;
  phimax_roi_low=(roi_phi1-1)/phi_spacing+1;
  zmax_roi_low=(roi_z1-1)/z_spacing+1; //exlcuded upper edge of our region of interest, measured in l-bins
  //that is, look at the last l-bin that is needed, and add 1 to that.

  //calculate the dimensions, in l-bins in our region of interest
  nr_roi_low=rmax_roi_low-rmin_roi_low;
  nphi_roi_low=phimax_roi_low-phimin_roi_low;
  nz_roi_low=zmax_roi_low-zmin_roi_low;

  printf("AnnularFieldSim::AnnularFieldSim building Epartial_lowres with  nr_roi_low=%d nphi_roi_low=%d nz_roi_low=%d nr_low=%d nphi_low=%d nz_low=%d =~%2.2fM TVector3 objects\n",nr_roi_low,nphi_roi_low,nz_roi_low,nr_low,nphi_low,nz_low,
	 nr_roi_low*nphi_roi_low*nz_roi_low*nr_low*nphi_low*nz_low/(1.0e6));
  printf("TVector is %lu bytes\n",sizeof(TVector3));
  //create the low-res, absolute 'lookup' grid to compute the field in each l-bin in the roi given charge in each other l-bin.
  Epartial_lowres=new MultiArray<TVector3>(nr_roi_low,nphi_roi_low,nz_roi_low,nr_low,nphi_low,nz_low);
  for (int i=0;i<Epartial_lowres->Length();i++)
    Epartial_lowres->GetFlat(i)->SetXYZ(0,0,0);

  //create a holder for the charge in each l-bin, to save time later.
  q_lowres=new MultiArray<float>(nr_low,nphi_low,nz_low);
  for (int i=0;i<q_lowres->Length();i++)
    *(q_lowres->GetFlat(i))=0;



  return;
}

AnnularFieldSim::AnnularFieldSim(float in_innerRadius, float in_outerRadius, float in_outerZ,
				 int r, int roi_r0, int roi_r1, 
				 int phi, int roi_phi0, int roi_phi1, 
				 int z, int roi_z0, int roi_z1,
				 float vdr)
  :
  AnnularFieldSim(in_innerRadius, in_outerRadius, in_outerZ,
		  r, roi_r0, roi_r1, 1,1,
		  phi, roi_phi0, roi_phi1, 1,1,
		  z, roi_z0, roi_z1, 1,1,
		  vdr)
  {
    //just passing through for the old version again
    //creates a region with high-res size of 1 (just a single cell) and low-res spacing of 1, which ought to match the behavior (with a little more overhead) from when there was no highres-lowres distinction
    printf("AnnularFieldSim::OldConstructor building AnnularFieldSim with local_size=1 in all dimensions, lowres_spacing=1 in all dimensions\n");
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
  printf("AnnularFieldSim::OldConstructor building AnnularFieldSim with roi=full in all dimensions\n");
  return;
}

TVector3 AnnularFieldSim::calc_unit_field(TVector3 at, TVector3 from){
  //if(debugFlag()) printf("%d: AnnularFieldSim::calc_unit_field(at=(r=%f,phi=%f,z=%f))\n",__LINE__,at.Perp(),at.Phi(),at.Z());
  int r_position;
  if (GetRindexAndCheckBounds(at.Perp(), &r_position)!=InBounds){
    printf("something's asking for 'at' with r=%f, which is index=%d\n",at.Perp(),r_position);
    assert(1==2);
  }
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
int AnnularFieldSim::FilterPhiIndex(int phi,int range=-1){
  if (range<0) range=nphi; //default input is range=-1.  in that case, use the intrinsic resolution of the q grid.
		 
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
    assert(1==2);
  }
  return p;
}
  
AnnularFieldSim::BoundsCase AnnularFieldSim::GetRindexAndCheckBounds(float pos, int *r){
  //if(debugFlag()) printf("%d: AnnularFieldSim::GetRindexAndCheckBounds(r=%f)\n",__LINE__,pos);

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
  // if(debugFlag()) printf("%d: AnnularFieldSim::GetPhiIndexAndCheckBounds(phi=%f)\n\n",__LINE__,pos);
  float p0f=(pos)/step.Phi(); //the position in phi, in units of step, starting from the low edge of the 0th bin.
  int phitemp=floor(p0f);
  int p0=FilterPhiIndex(phitemp);
  *phi=p0;

   phitemp=floor(p0f-ALMOST_ZERO);
  int p0lowered_slightly=FilterPhiIndex(phitemp);
   phitemp=floor(p0f+ALMOST_ZERO);
  int p0raised_slightly=FilterPhiIndex(phitemp);
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
  //if(debugFlag()) printf("%d: AnnularFieldSim::GetZindexAndCheckBounds(z=%f)\n\n",__LINE__,pos);
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
  //if(debugFlag()) printf("%d: AnnularFieldSim::fieldIntegral(x=%f,y=%f, z=%f) to z=%f\n\n",__LINE__,start.X(),start.Y(),start.Z(),zdest);

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
  TVector3 tf;
  for(int i=zi;i<zf;i++){ //count the whole cell of the lower end, and skip the whole cell of the high end.
    tf=Efield->Get(r-rmin_roi,phi-phimin_roi,i-zmin_roi);
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
    assert(1==2);
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

  float hr,hphi,hz;//the center of the histogram bin
  int localr,localphi,localz;//the f-bin of our local structure that contains that center.
  for (int i=(rmin-hrmin)/hrstep;i<hrn;i++){
    hr=hrmin+hrstep*(i+0.5);//bin center
    localr=(hr-rmin)/step.Perp();
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
      hphi=hphimin+hphistep*(j+0.5); //bin center
      localphi=hphi/step.Phi();
      if (localphi>=nphi){//handle wrap-around:
	localphi-=nphi;
      }
      if (localphi<0){//handle wrap-around:
	localphi+=nphi;
      }
      //todo:  should add ability to take in a phi- slice only
      for (int k=(zmin-(hzmin-zoffset))/hzstep;k<hzn;k++){
	hz=hzmin-zoffset+hzstep*(k+0.5);//bin center
	localz=hz/step.Z();
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
	//if(debugFlag()) printf("%d: AnnularFieldSim::load_spacecharge adding Q=%f from hist(%d,%d,%d) into cell (%d,%d,%d)\n",__LINE__,qbin,i,j,k,localr,localphi,localz);
	q->Add(localr,localphi,localz,qbin);
      }
    }
  }

  //go through the q array and build q_lowres.
  for (int i=0;i<q_lowres->Length();i++)
    *(q_lowres->GetFlat(i))=0;

  //fill our low-res
  //note that this assumes the last bin is short or normal length, not long.
  for (int ifr=0;ifr<nr;ifr++){
    int r_low=ifr/r_spacing; //index of our l-bin is just the integer division of the index of our f-bin
    for (int ifphi=0;ifphi<nphi;ifphi++){
      int phi_low=ifphi/phi_spacing;
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
  
  TVector3 localF;//holder for the summed field at the current position.
  for (int ir=rmin_roi;ir<rmax_roi;ir++){
    for (int iphi=phimin_roi;iphi<phimax_roi;iphi++){
      for (int iz=zmin_roi;iz<zmax_roi;iz++){
	localF=sum_field_at(ir,iphi,iz);
	Efield->Set(ir-rmin_roi,iphi-phimin_roi,iz-zmin_roi,localF);
	//if (localF.Mag()>1e-9)
	if(debugFlag()) printf("%d: AnnularFieldSim::populate_fieldmap fieldmap@ (%d,%d,%d) mag=%f\n",__LINE__,ir,iphi,iz,localF.Mag());
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
  int r_highres_dist=(nr_high-1)/2;
  int phi_highres_dist=(nphi_high-1)/2;
  int z_highres_dist=(nz_high-1)/2;



  static int ncellsin[3][3][3];
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      for (int k=0;k<3;k++){
	ncellsin[i][j][k]=0; //we could count total volume, but wihtout knowing the charge prior, it's not clear that'd be /better/
      }
    }
  }

  //todo: if this runs too slowly, I can do geometry instead of looping over all the cells that are possibly in range


  TVector3 currentf, newf; //the averaged field vector without, and then with the new contribution from the f-bin being considered.
  for (int ifr=rmin_roi;ifr<rmax_roi;ifr++){
    for (int ifphi=phimin_roi;ifphi<phimax_roi;ifphi++){
      for (int ifz=zmin_roi;ifz<zmax_roi;ifz++){
	//if(debugFlag()) printf("%d: AnnularFieldSim::populate_highres_lookup icell=(%d,%d,%d)\n",__LINE__,ifr,ifphi,ifz);

	at=GetCellCenter(ifr, ifphi, ifz);
	//define the parent l-bin cells we're dealing with here:
	int r_parentlow=floor((ifr-r_highres_dist)/(r_spacing*1.0));//partly enclosed in our high-res region
	int r_parenthigh=floor((ifr+r_highres_dist)/(r_spacing*1.0))+1;//definitely not enclosed in our high-res region
	int phi_parentlow=floor((ifphi-phi_highres_dist)/(phi_spacing*1.0));
	int phi_parenthigh=floor((ifphi+phi_highres_dist)/(phi_spacing*1.0))+1; //note that this can be bigger than nphi!  We keep track of that.
	int z_parentlow=floor((ifz-z_highres_dist)/(z_spacing*1.0));
	int z_parenthigh=floor((ifz+z_highres_dist)/(z_spacing*1.0))+1;


	//for every f-bin in the l-bins we're dealing with, figure out which relative highres bin it's in, and average the field vector into that bin's vector
	//note most of these relative bins have exactly one f-bin in them.  It's only the edges that can get more.
	//note this is a running average:  Anew=(Aold*Nold+V)/(Nold+1) and so on.
	//note also that we automatically skip f-bins that would've been out of the valid overall volume.
	for (int ir=r_parentlow*r_spacing;ir<r_parenthigh*r_spacing;ir++){
	  //skip parts that are out of range:
	  if (ir<0) ir=0;
	  if (ir>=nr) break;
	  
	  int rbin=(ir-ifr)+r_highres_dist;//zeroth bin when we're at max distance below, etc.
	  int rcell=1;
	  if (rbin<=0) {
	    rbin=0;
	    rcell=0;
	  }
	  if (rbin>=nr_high){
	    rbin=nr_high-1;
	    rcell=2;
	  }
	  for (int iphi=phi_parentlow*phi_spacing;iphi<phi_parenthigh*phi_spacing;iphi++){
	    //no phi out-of-range checks since it's circular.
	    int phibin=(iphi-ifphi)+phi_highres_dist;
	    int phicell=1;
	    if (phibin<=0) {
	      phibin=0;
	      phicell=0;
	    }
	    if (phibin>=nphi_high){
	      phibin=nphi_high-1;
	      phicell=2;
	    }
	    for (int iz=z_parentlow*z_spacing;iz<z_parenthigh*z_spacing;iz++){
	      if (iz<0) iz=0;
	      if (iz>=nz) break;
	      int zbin=(iz-ifz)+z_highres_dist;
	      int zcell=1;
	      if (zbin<=0) {
		zbin=0;
		zcell=0;
	      }
	      if (zbin>=nz_high){
		zbin=nz_high-1;
		zcell=2;
	      }
	      //'from' is in absolute coordinates
	      from=GetCellCenter(ir, FilterPhiIndex(iphi), iz);
	      ncellsin[rcell][phicell][zcell]++;
	      int nc= ncellsin[rcell][phicell][zcell];

	      //but Epartial is in coordinates relative to the roi
	      if (ifphi-phimin_roi<0)printf("%d: Getting with phi=%d\n",__LINE__,ifphi-phimin_roi);
	      currentf=Epartial_highres->Get(ifr-rmin_roi,ifphi-phimin_roi,ifz-zmin_roi,rbin,phibin,zbin);
	      //to keep this as the average, we multiply what's there back to its initial summed-but-not-divided value
	      //then add our new value, and the divide the new sum by the total number of cells
	      newf=(currentf*(nc-1)+calc_unit_field(at,from))*(1/(nc*1.0)); 
	      Epartial_highres->Set(ifr-rmin_roi,ifphi-phimin_roi,ifz-zmin_roi,rbin,phibin,zbin,newf);
	    }
	  }
	}
      }
    }
  }
  return;
}

void AnnularFieldSim::populate_lowres_lookup(){

  TVector3 at(1,0,0);
  TVector3 from(1,0,0);

  int fr_low,fr_high,fphi_low,fphi_high,fz_low,fz_high;//edges of the outer l-bin
  int r_low,r_high,phi_low,phi_high,z_low,z_high;//edges of the inner l-bin

  for (int ifr=rmin_roi_low;ifr<rmax_roi_low;ifr++){
    fr_low=ifr*r_spacing;
    fr_high=fr_low+r_spacing-1;
    if (fr_high>=nr) fr_high=nr-1;	
    for (int ifphi=phimin_roi_low;ifphi<phimax_roi_low;ifphi++){
      fphi_low=ifphi*phi_spacing;
      fphi_high=fphi_low+phi_spacing-1;
      if (fphi_high>=nphi) fphi_high=nphi-1;
      for (int ifz=zmin_roi_low;ifz<zmax_roi_low;ifz++){
	fz_low=ifz*z_spacing;
	fz_high=z_low+z_spacing-1;
	if (fz_high>=nz) fz_high=nz-1;
	at=GetGroupCellCenter(fr_low,fr_high,fphi_low,fphi_high,fz_low,fz_high);
	//printf("ifr=%d, rlow=%d,rhigh=%d,r_spacing=%d\n",ifr,r_low,r_high,r_spacing);
	//if(debugFlag())	  printf("%d: AnnularFieldSim::populate_lowres_lookup icell=(%d,%d,%d)\n",__LINE__,ifr,ifphi,ifz);
		
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
	      //this calc's okay.
	      Epartial_lowres->Set(ifr-rmin_roi_low,ifphi-phimin_roi_low,ifz-zmin_roi_low,ior,iophi,ioz,calc_unit_field(at,from));
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
  // if(debugFlag()) printf("%d: AnnularFieldSim::sum_field_at(r=%d,phi=%d, z=%d)\n",__LINE__,r,phi,z);

  /*
    for the near future:
  TVector3 sum=(sum_local_field_at(r, phi, z)
		+ sum_nonlocal_field_at(r,phi,z)
		+ Escale*Eexternal->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi));
  return sum;
  */
		


  

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

  int r_highres_dist=(nr_high-1)/2;
  int phi_highres_dist=(nphi_high-1)/2;
  int z_highres_dist=(nz_high-1)/2;

  int r_parentlow=floor((r-r_highres_dist)/(r_spacing*1.0));//partly enclosed in our high-res region
  int r_parenthigh=floor((r+r_highres_dist)/(r_spacing*1.0))+1;//definitely not enclosed in our high-res region
  int phi_parentlow=floor((phi-phi_highres_dist)/(phi_spacing*1.0));
  int phi_parenthigh=floor((phi+phi_highres_dist)/(phi_spacing*1.0))+1; //note that this can be bigger than nphi!  We keep track of that.
  int z_parentlow=floor((z-z_highres_dist)/(z_spacing*1.0));
  int z_parenthigh=floor((z+z_highres_dist)/(z_spacing*1.0))+1;
  //printf("AnnularFieldSim::sum_field_at parents: rlow=%d,philow=%d,zlow=%d,rhigh=%d,phihigh=%d,zhigh=%d\n",r_parentlow,phi_parentlow,z_parentlow,r_parenthigh,phi_parenthigh,z_parenthigh);

  //zero our current qlocal holder:
  for (int i=0;i<q_local->Length();i++)
    *(q_local->GetFlat(i))=0;

  //get the charge involved in the local highres block:
  for (int ir=r_parentlow*r_spacing;ir<r_parenthigh*r_spacing;ir++){
    if (ir<0) ir=0;
    if (ir>=nr) break;
    int rbin=(ir-r)+r_highres_dist;//zeroth bin when we're at max distance below, etc.
    if (rbin<0) rbin=0;
    if (rbin>=nr_high) rbin=nr_high-1;
    for (int iphi=phi_parentlow*phi_spacing;iphi<phi_parenthigh*phi_spacing;iphi++){
      //no phi range checks since it's circular.
      int phibin=(iphi-phi)+phi_highres_dist;
      if (phibin<0) phibin=0;
      if (phibin>=nphi_high) phibin=nphi_high-1;
      for (int iz=z_parentlow*z_spacing;iz<z_parenthigh*z_spacing;iz++){
	if (iz<0) iz=0;
	if (iz>=nz) break;
	//if(debugFlag()) printf("%d: AnnularFieldSim::sum_field_at, reading q in f-bin at(r=%d,phi=%d, z=%d)\n",__LINE__,ir,iphi,iz);

	int zbin=(iz-z)+z_highres_dist;
	if (zbin<0) zbin=0;
	if (zbin>=nz_high) zbin=nz_high-1;
	//printf("filtering in local highres block\n");
	q_local->Add(rbin,phibin,zbin,q->Get(ir,FilterPhiIndex(iphi),iz));
	//printf("done filtering in local highres block\n");

      }
    }
  }

  //now q_highres is up to date for our current cell of interest.
  //start building our full sum by scaling the local lookup table by q.
  //note that the lookup table needs to have already accounted for cell centers.

  TVector3 sum(0,0,0);

  //note that Epartial_highres returns zero if we're outside of the valid region.  q_local will also be zero there.
  //these are loops over the position in the epartial highres grid, so relative to the point in question:
  for (int ir=0;ir<nr_high;ir++){
    for (int iphi=0;iphi<nphi_high;iphi++){
      for (int iz=0;iz<nz_high;iz++){
	//first three are relative to the roi, last three are relative to the point in the first three.  ooph.
	if (phi-phimin_roi<0)	printf("%d: Getting with phi=%d\n",__LINE__,phi-phimin_roi);
	sum+=Epartial_highres->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi,ir,iphi,iz)*q_local->Get(ir,iphi,iz);
      }
    }
  }

  //now we find our interpolated position between the eight nearby cells:
  //this is similar to the interpolated integral stuff, except we're at integer steps inside integer blocks
  //and we have z as well, now.
  bool skip[]={false,false,false,false,false,false,false,false};

  
  float r0=r/(1.0*r_spacing)-0.5-rmin_roi_low; //the position in r, in units of r_spacing, starting from the center of the 0th l-bin in the roi.
  int r0i=floor(r0); //the integer portion of the position. -- what center is below our position?
  float r0d=r0-r0i;//the decimal portion of the position. -- how far past center are we?
  //instead of listing all eight, I'll address these as i%2, (i/2)%2 and (i/4)%2 to avoid typos
  int ri[2];//the r position of the eight cell centers to consider.
  ri[0]=r0i;
  ri[1]=r0i+1;
  float rw[2];//the weight of that cell, linear in distance from the center of it
  rw[0]=1-r0d; //1 when we're on it, 0 when we're at the other one.
  rw[1]=r0d; //1 when we're on it, 0 when we're at the other one


  //zero out if the requested l-bin is out of the roi (happens if we're on the edge of the roi) 
  if (ri[0]<0){
    for (int i=0;i<8;i++)
      if ((i/4)%2==0)
	skip[i]=true; // don't handle contributions from ri[0].
    rw[1]=1; //and weight like we're dead-center on the outer cells.
  }
  if (ri[1]>=nr_roi_low){
    for (int i=0;i<8;i++)
      if ((i/4)%2==1)
	skip[i]=true; // don't bother handling ri[1].
    rw[0]=1; //and weight like we're dead-center on the inner cells.
  }
  
  //now repeat that structure for phi:

  float p0=phi/(1.0*phi_spacing)-0.5-phimin_roi_low; //the position 'phi' in  units of phi_spacing, starting from the center of the 0th bin.
  //printf("prepping to filter low, p0=%f, phi=%d, phi_spacing=%d, phimin_roi_low=%d\n",p0,phi,phi_spacing,phimin_roi_low);

  int p0i=floor(p0); //the integer portion of the position. -- what center is below our position?
  float p0d=p0-p0i;//the decimal portion of the position. -- how far past center are we?
  int pi[4];//the phi position of the four cell centers to consider
  //printf("filtering low, p0i=%d, nphi_low=%d\n",p0i,nphi_low);
  pi[0]=FilterPhiIndex(p0i,nphi_low);
  pi[1]=FilterPhiIndex(p0i+1,nphi_low);
  //printf("done filtering low\n");
  float pw[2];//the weight of that cell, linear in distance from the center of it
  pw[0]=1-p0d; //1 when we're on it, 0 when we're at the other one.
  pw[1]=p0d; //1 when we're on it, 0 when we're at the other one

  //note that since phi wraps around itself, we have the possibility of being outside by being above/below the opposite end of where we thought we were
  //if we completely cover phi, then we'll still be in the roi.
  if (pi[0]<0 || pi[0]>=nphi_roi_low){
   for (int i=0;i<8;i++)
      if ((i/2)%2==0)
	skip[i]=true; // don't bother handling phii[1].
    pw[1]=1; //and weight like we're dead-center on the high cells. 

  }
  if (pi[1]>=nphi_roi_low || pi[1]<0){
    for (int i=0;i<8;i++)
      if ((i/2)%2==1)
	skip[i]=true; // don't bother handling phii[1].
    pw[0]=1; //and weight like we're dead-center on the high cells. 
 }

  //and once more for z.  ooph.
  
  float z0=z/(1.0*z_spacing)-0.5-zmin_roi_low; //the position in r, in units of r_spacing, starting from the center of the 0th bin.
  int z0i=floor(z0); //the integer portion of the position. -- what center is below our position?
  float z0d=z0-z0i;//the decimal portion of the position. -- how far past center are we?
  //instead of listing all eight, I'll address these as i%2, (i/2)%2 and (i/4)%2 to avoid typos
  int zi[2];//the r position of the eight cell centers to consider.
  zi[0]=z0i;
  zi[1]=z0i+1;
  float zw[2];//the weight of that cell, linear in distance from the center of it
  zw[0]=1-z0d; //1 when we're on it, 0 when we're at the other one.
  zw[1]=z0d; //1 when we're on it, 0 when we're at the other one


  
  if (zi[0]<0){
    for (int i=0;i<8;i++)
      if ((i)%2==0)
	skip[i]=true; // don't bother handling zi[0].
    zw[1]=1; //and weight like we're dead-center on the higher cells.
  }
  if (zi[1]>=nz_roi_low){
    for (int i=0;i<8;i++)
      if ((i)%2==1)
	skip[i]=true; // don't bother handling zi[1].
    zw[0]=1; //and weight like we're dead-center on the lower cells.
  }
  
  int lBinEdge[2];//lower and upper (included) edges of the low-res bin, measured in f-bins, reused per-dimension
  int fBinEdge[2];//lower and upper edge of the high-res region, measured in f-bins, reused per-dimension.  I shouldn't call this 'fbin', since it's not a single bin.
  bool overlapsR, overlapsPhi,overlapsZ; //whether we overlap in R, phi, and z.
  for (int ir=0;ir<nr_low;ir++){
    lBinEdge[0]=ir*r_spacing;
    lBinEdge[1]=(ir+1)*r_spacing-1;
    fBinEdge[0]=r-r_highres_dist;
    fBinEdge[1]=r+r_highres_dist;
    overlapsR= (lBinEdge[0]<=fBinEdge[1] && fBinEdge[0]<=lBinEdge[1]);
    for (int iphi=0;iphi<nphi_low;iphi++){
      lBinEdge[0]=iphi*phi_spacing;
      lBinEdge[1]=(iphi+1)*phi_spacing-1;
      fBinEdge[0]=phi-phi_highres_dist;
      fBinEdge[1]=phi+phi_highres_dist;
      overlapsPhi= (lBinEdge[0]<=fBinEdge[1] && fBinEdge[0]<=lBinEdge[1]);
      for (int iz=0;iz<nz_low;iz++){
	lBinEdge[0]=iz*z_spacing;
	lBinEdge[1]=(iz+1)*z_spacing-1;
	fBinEdge[0]=z-z_highres_dist;
	fBinEdge[1]=z+z_highres_dist;
	overlapsZ= (lBinEdge[0]<=fBinEdge[1] && fBinEdge[0]<=lBinEdge[1]);
	//conceptually: see if the l-bin overlaps with the high-res region:
	//the high-res region includes all indices from r-(nr_high-1)/2 to r+(nr_high-1)/2.
	//each low-res region includes all indices from ir*r_spacing to (ir+1)*r_spacing-1.
	if( overlapsR && overlapsPhi && overlapsZ){
	  //if their bounds are interleaved in all dimensions, there is overlap, and we've already summed this region.
	  continue;
	} else {
	  //if(debugFlag()) printf("%d: AnnularFieldSim::sum_field_at, considering l-bin at(r=%d,phi=%d, z=%d)\n",__LINE__,ir,iphi,iz);

	  for (int i=0;i<8;i++){
	    if (skip[i]) continue;
	    //the ri, pi, and zi elements are relative to the roi, as needed for Epartial.
	    //the ir, iphi, and iz are all absolute, as needed for q_lowres
	    if (pi[(i/2)%2]<0) printf("%d: Getting with phi=%d\n",__LINE__,pi[(i/2)%2]);
	    sum+=(Epartial_lowres->Get(ri[(i/4)%2],pi[(i/2)%2],zi[(i)%2],ir,iphi,iz)
		  *q_lowres->Get(ir,iphi,iz))*zw[(i)%2]*pw[(i/2)%2]*rw[(i/4)%2];
	  }
	}
      }
    }
  }
  sum+=Escale*Eexternal->Get(r-rmin_roi,phi-phimin_roi,z-zmin_roi);
  if(debugFlag()) printf("summed field at (%d,%d,%d)=(%f,%f,%f)\n",r,phi,z,sum.X(),sum.Y(),sum.Z());
  
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
