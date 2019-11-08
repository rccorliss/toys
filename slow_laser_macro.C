void slow_laser_macro(){

      const int nLasers=12;

  //define the basis vectors local to the laser:
    TVector3 laser_nominal[nLasers];
    TVector3 laser_transverse[nLasers];

  //define the position of the lasers:
    TVector3 laser_position[nLasers];//(0,40.0,0); //position laser emitter on the y axis.
    float laser_position_angle0=0;
    
    float angle_increment=2*TMath::Pi()/nLasers; // assume they're equally spaced
    float laser_tilt=5*TMath::Pi()/180; //and have a common tilt, leaning outward in the rz plane,
    
    //rotate each laser to the proper position and update the beam nominal and transverse direction.
    for (int i=0;i<nLasers;i++)
    {
        laser_nominal[i].SetXYZ(0,0,1);
        laser_transverse[i].SetXYZ(1,0,0);
        laser_position[i].SetXYZ(0,40.0,0);
    //rotate the basis vectors of the laser:
    laser_nominal[i].RotateX(-laser_tilt);
    laser_nominal[i].RotateZ(laser_position_angle0+angle_increment*i);
    laser_transverse[i].RotateZ(laser_position_angle0+angle_increment*i);
    //rotate the position of the laser:
    laser_position[i].RotateZ(laser_position_angle0+angle_increment*i);
    }
    
    
    
    
  //define the surface we are illuminating:
  TVector3 surface_point(0,0,100.0);
  TVector3 surface_normal(0,0,1);


 
  //set up the distributions for theta and phi angles relative to laser_nominal:
  //TF1 *laser_theta=new TF1("laser_theta","exp([0]*x)",0,TMath::Pi()/2.0);
  //TF1 *laser_theta=new TF1("laser_test","exp(-[0]*(x-[1])^2)",0,TMath::Pi()/2.0);
  TF1 *laser_phi=new TF1("laser_phi","1",0,2*TMath::Pi());
  //laser_theta->SetParameters(100,0.5);
  //laser_theta->Draw();
    
  //laser intensity as a function of distance from the beam center:
    float w0=3;
    float zsample=1.0;//all in cm.
    float lambda=600*1e-7;//10^-7 converts from nm in cm.
    float denom=w0*w0+lambda*lambda*zsample*zsample/(TMath::Pi()*TMath::Pi()*w0*w0);
    float thetamax=TMath::Pi()/4;
    float samplemax=100; //
    float samplemax2=TMath::Tan(thetamax)*zsample; // max radial distance = max angle at zsample distance.
    printf("samplemax2=%f\n",lambda);
    if (samplemax2 > samplemax)
        samplemax=samplemax2;
    TF1 *laser_r=new TF1("laser_r","x*exp(-2*x^2/[0])",0,samplemax); // intensity as function of r (note we have performed the integral over phi!)
        //TF1 *laser_r=new TF1("laser_r","x/[0]",0,TMath::Tan(thetamax)*zsample); // intensity as function of r
    laser_r->SetParameter(0,denom);
  //return;

    //TODO: find a realistic waist
    //TODO: define a lens?
    //TODO: add a legend to the plots so they show the variables used!
    //could weight by intensity if we wanted.
    //could invert the calculation to get the beam position as a function of point of intersection

    
 //plane in vector notation:
 //(v-p0).n=0 -- the position of point v relative to the defined point has no component out of the plane 
//line in vector notation:
 //v=t*L+L0 -- the position of point v is the starting position plus some multiple of its direction vector.
 //sub in:
 //(t*L+L0-p0).n=0
 //t*(L.n)+L0.n-p0.n=0
 //t=(p0-L0).n/(L.n) -- as long as L.n!=0;
 //the intersection point is:
 //v={(p0-L0).n/(L.n)}*L+L0
 //(and so the /starting/ point as a function of intersection point is:
 //L0=v-{(p0-L0).n/(L.n)}*L
 //L0=(v-p0.n/(L.n)*L)+L0.n/(L.n)*L
 //since in our case the laser position has z=0,  L0.n=Lz*zhat=0
 //L0=v-p0.n/(L.n)*L
 // )
 //conceptually, point backward from the intersection point in the negative laser direction until we hit the plane defined by the laser direction and the laser position. This tells us the 
 

 TH2F *hPhotonAtSurface=new TH2F("hPhotonAtSurface","hPhotonAtSurface",40,-80,80,40,-80,80);
 TH2F *hPhotonAngle=new TH2F("hPhotonAngle","hPhotonAngle;#theta (x);#phi (y)",50,0,2,50,0,6.5);
 TH2F *hPhotonDirection=new TH2F("hPhotonDirection","hPhotonDirection; (x); (y)",50,-2,2,50,-2,2);
 TH1F *hPhoton=new TH1F("hPhoton","hPhoton",100,0,TMath::Pi());
 
 int nPhotons=4000000;
 for (int i=0; i<nPhotons;i++)  {
     

     
     

   int L=i%nLasers;
     
     
     //hPhoton->Fill(laser_shape->GetRandom(),1.0/nPhotons);
     //float photon_theta=laser_theta->GetRandom();
     float photon_phi=laser_phi->GetRandom();
     float photon_r=laser_r->GetRandom();
     
     TVector3 photon_direction=laser_nominal[L];
     TVector3 photon_position=laser_position[L];
     TVector3 photon_offset=laser_transverse[L]*photon_r;
     photon_offset.Rotate(photon_phi,laser_nominal[L]);
     photon_position=photon_position+photon_offset;
     //hPhotonAtSurface->Fill(photon_direction.X(),photon_direction.Y());
     //hPhotonAtSurface->Fill(photon_position.X(),photon_position.Y());
     //continue;

     float denominator=photon_direction.Dot(surface_normal);
     if (abs(denominator)<0.001) { //parallel to surface.  no solution}
         continue; //skip to the next iteration through the loop
     }
     TVector3 intersect=(surface_normal.Dot(surface_point-photon_position)/denominator)*photon_direction+photon_position;
     hPhotonAtSurface->Fill(intersect.X(),intersect.Y());
     
     /*
     
     
     
   //hPhoton->Fill(laser_shape->GetRandom(),1.0/nPhotons);
   //float photon_theta=laser_theta->GetRandom();
   float photon_phi=laser_phi->GetRandom();
     float photon_r=laser_r->GetRandom();
     float  photon_theta=TMath::ATan2(photon_r,zsample);
     
   TVector3 photon_direction=laser_nominal[L];
   hPhotonAngle->Fill(photon_theta,photon_phi);//fill x and y.  annoying that it's backward from the arguments earlier, but what're you gonna do? ;)
   photon_direction.Rotate(photon_theta,laser_transverse[L]);
   photon_direction.Rotate(photon_phi,laser_nominal[L]);
   hPhotonDirection->Fill(photon_direction.X(),photon_direction.Y());

   float denominator=photon_direction.Dot(surface_normal);
   if (abs(denominator)<0.001) { //parallel to surface.  no solution}
     continue; //skip to the next iteration through the loop
   }
     TVector3 intersect=(surface_normal.Dot(surface_point-laser_position[L])/denominator)*photon_direction+laser_position[L];
   hPhotonAtSurface->Fill(intersect.X(),intersect.Y());
      */
 }
 //hPhotonDirection->Draw("colz");
 //hPhotonAngle->Draw("colz");
 hPhotonAtSurface->Draw("colz");

     
     
     //hPhoton->Draw("same");
return;
}

/* things you can look up/do:
rebuild in matlab
figure out rotations needed.
get some testable cases to confirm
pick an interesting option or two:
gaussian from last time at (0,0,1) check it looks the same
gaussian from last time at 30 degrees radially outward.
figure out if there is a root package for lookign at where lines cross planes -- if not, write the equations for it.
*/
