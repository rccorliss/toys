

void vdrift_comparison(){
  string titles[]={"angle[rad]","V_E","V_EXB","V_B"};
  //this data is from an email sent to Ross Corliss by Prakhar Garg on Feb 12, 2020
  //velocity units are cm/ns.  Settings are for E=400V/cm and B=1.4T
    //for NE_CF4 90:10
  float data[]={0,0,0,0,
	    0.174522,0.00782275,0.00042635,0.00121795,
	    0.349044,0.00741239,0.000846258,0.00237184,
	    0.523567,0.00674057,0.00125054,0.00339112,
	    0.698089,0.00583379,0.00163049,0.00420286,
	    0.872611,0.0047286,0.00197579,0.00471865,
	    1.04713,0.00348926,0.00227464,0.00483343,
	    1.22166,0.00221189,0.00251653,0.00439498,
	    1.39618  ,0.0010535,0.00269209,0.0031193,
	    1.5707,0.000371075,0.0027859,1.78E-06};


  //for NE_CF4 50:50
  /*
  float data[][]={{0,0.00822536,0,0},
	    {0.174522,0.00807386,0.000432942,0.00128091},
	    {0.349044,0.00762722,0.000855897,0.00248636},
	    {0.523567,0.00690121,0.00125955,0.00353749},
	    {0.698089,0.00593088,0.00163436,0.0043522},
	    {0.872611,0.00476267,0.0019724,0.00483794},
	    {1.04713,0.00346861,0.00226569,0.0048791},
	    {1.22166,0.0021542,0.00250777,0.00430751},
	    {1.39618,0.00100473,0.0026822,0.0028305},
	    {1.5707,0.000448308,0.0027534,1.55E-06}};
  */
  float B=14;//kGauss
  float E=400;//V/cm
  // float vdrift=data[1*4+1]
  int n=10;
  float magExB[n];
  float vTnorm[n];
  float vE[n];
  float vBnorm[n];
  float angledeg[n];
  float omtau[n];
  float omtauB[n];
  float omtau_analytic[n];
  
  printf("angle(deg)\tvE      \tvExBnorm\tomegatau\tomtau from B\tanalytic\n");
  for (int i=0;i<10;i++){
    magExB[i]=sin(data[i*4+0]);
    float magEdotB=cos(data[i*4+0]);
    if (magExB[i]<1e-8) {
      vTnorm[i]=0;
    }else{
      vTnorm[i]=data[i*4+2]/magExB[i];
    }
    
    if (magEdotB*data[i*4+2]<1e-8) {
      omtauB[i]=0;
    }else{
      omtauB[i]=data[i*4+3]/data[i*4+2]/magEdotB;
    }
    vE[i]=data[i*4+1];
    angledeg[i]=data[i*4+0]*180/3.14;
    if (vE[i]<1e-8) {
      omtau[i]=0;
    }else{
    omtau[i]=vTnorm[i]/vE[i];
    }
    omtau_analytic[i]=10*(B)*(vE[i]*1e3)/E;

    printf("%8E\t%8E\t%8E\t%8E\t%8E\t%8E\n",angledeg[i],vE[i],vTnorm[i],omtau[i],omtauB[i],omtau_analytic[i]);
  }

  TGraph *omtau_v_ang=new TGraph(9,angledeg+1,omtau+1);
  omtau_v_ang->SetTitle("#omega#tau from Magboltz E/ExB;#theta (deg);#omega#tau");
  omtau_v_ang->SetLineColor(kRed);
  TGraph *omtauB_v_ang=new TGraph(9,angledeg+1,omtauB+1);
  omtauB_v_ang->SetTitle("#omega#tau from Magboltz B/ExB;#theta (deg);#omega#tau");
  omtauB_v_ang->SetLineColor(kGreen);
  TGraph *omtau_ana_v_ang=new TGraph(9,angledeg+1,omtau_analytic+1);
 omtau_ana_v_ang->SetTitle("#omega#tau from formula;#theta (deg);#omega#tau");
  omtau_ana_v_ang->SetLineColor(kBlue);
  TCanvas *c=new TCanvas();
  omtau_v_ang->Draw();
  omtauB_v_ang->Draw("same");
  omtau_ana_v_ang->Draw("same");
  c->BuildLegend();
  return;
}
