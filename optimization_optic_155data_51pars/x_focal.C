/*===============================================================
 * Nguyen Ton Sep 12th 2016
 * Jan 2018: update on optic data, include 2.2 GeV
 * Purpose: optimize forward matrix for xfoc (in meter unit)
 * Input data: 155 data include 1.1, 1.5, 2.2 GeV
 * Output: coefficients for forward matrix
 *===============================================================*/
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TMath.h"
#include <iostream>
#include <fstream>
#include "TStyle.h"
#include <string>
int const dim=155;
int const numpar=14;
Float_t theta[dim],phi[dim],dp[dim],ytg[dim],xfoc[dim],thfoc[dim],errorz[dim],P[dim];

using namespace std;

Double_t xfunc(float delta, float yt, float phit, float thetat,float pcent, Double_t *par)
{
  Double_t value= par[0]
    +par[1]*delta
    +par[2]*delta*delta
    +par[3]*phit
    +par[4]*phit*phit
    +par[5]*yt
    +par[6]*thetat
    +par[7]*thetat*thetat
    +par[8]*thetat*yt
    +par[9]*thetat*phit
    +par[10]*pcent
    +par[11]*pcent*thetat
    +par[12]*pcent*phit
    +par[13]*pcent*yt
    ;
 return value;
}

//===============================
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t nbins = dim;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t ddelta;
   for (i=0;i<nbins; i++) {
     ddelta  = (xfoc[i]-xfunc(dp[i],ytg[i],phi[i],theta[i],P[i],par))/errorz[i];
     chisq+=ddelta*ddelta;
     double x_test=xfunc(dp[i],ytg[i],phi[i],theta[i],P[i],par);
     cout<<i<<"\t"<<xfoc[i]<<"\t"<<x_test<<"\t"<<xfoc[i]-x_test<<endl;
   }
   f = chisq;
}

void x_focal()
{
     ifstream data("full_kin_155_pnts.dat",ios_base::in);
     double th_tg,ph_tg,y_tg,dpp,yfp,phfp,xfp,thfp,err_yfp,err_xfp,err_phfp,err_thfp,P0;
     int count=0;

     TGraph * gr = new TGraph();
     gr->SetTitle("Combine all holes in Focal plane;y_{fp}(m);#phi_{fp}(rad)");
     while(!data.eof()){

       data>>ph_tg>>th_tg>>y_tg>>dpp>>xfp>>err_xfp>>yfp>>err_yfp>>phfp>>err_phfp>>thfp>>err_thfp>>P0;
       theta[count] = th_tg;
       phi[count]   = ph_tg;
       ytg[count]   = y_tg;
       dp[count]    = dpp;
       xfoc[count]  = xfp;
       thfoc[count] = thfp;
       errorz[count]= err_xfp;
       P[count]     = P0;
       //gr->SetPoint(count,xfp,dpp);
       count++;
     }

     //start TMinuit
     // TCanvas * cc = new TCanvas("cc","",800,600);
     // cc->Clear();

     TMinuit *gMinuit = new TMinuit(numpar);  //initialize TMinuit with a maximum of 5 params
     //check numpar
     gMinuit->SetFCN(fcn);

   Double_t arglist[10];
   Int_t ierflg = 0;

   arglist[0] = 1;
   gMinuit->mnexcm("SET ERR", arglist ,1,ierflg);


   arglist[0]=1;
   gMinuit->mnexcm("SET PRIntout",arglist,1,ierflg);
// Set starting values and step sizes for parameters

   Double_t vstart[numpar] = 
{
  -0.3314 ,//offset
  13.95   ,//dp
  -24.64  ,//dp^2
  -7.06   ,//phi
   3.84   ,//phi^2
  3.765   ,//y
  -0.357  ,//theta
  -100.   ,//theta*theta
  120.    ,//theta*y
  250.    ,//theta*phi
  0.001   ,//p0
  0.001   ,//p0*theta
  1.0     ,//p0*phi
  1.0     //,//p0*y
  //  10.      //p0*phi*yt
};
   Double_t step[numpar] = 
{
  0.204/10. ,//offset
  13.95/10. ,//dp
  24.64/10. ,//dp^2
  7.06/10.  ,//phi
  3.84/10.  ,//phi^2
  3.765/10. ,//y
  0.357/10. ,//theta
  10.       ,//theta*theta
  120/10.   ,//theta*y
  250./10.  ,//theta*phi
  0.001     ,//p0
  0.1       ,//p0*theta
  0.001     ,//p0*phi
  0.001     //,//p0*yt
  //  1.0        //p0*phi*y
};

   for(int ii=0;ii<numpar;ii++)
     {
       gMinuit->mnparm(ii, "a", vstart[ii], step[ii], 0,0,ierflg);
     }

// Now ready for minimization step
   arglist[0] = 500000;
   arglist[1] = 0.1;

   gMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);

   double ParErr[numpar];
   for (int ii=0;ii<numpar;ii++){
     gMinuit->GetParameter(ii,vstart[ii],ParErr[ii]);
     cout<<vstart[ii]<<endl;//" \t"<<ParErr[ii]<<endl;
   }


  Double_t amin,edm,errdef;
   Int_t nvpar,nparx,icstat;
     gMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
     cout<<" minimum chi2 ====="<<amin<<endl;
     cout<<" err definition ="<<errdef<<endl;
     cout<<" number of variable parameter ="<<nvpar<<endl;
     cout<<" number of parameters ="<<nparx<<endl;
     cout<<icstat<<endl;
     delete gMinuit;
}
