/*===============================================================
 * Nguyen Ton Sep 12th 2016
 * Jan 2018: update on optic data, include 2.2 GeV
 * Purpose: optimize forward matrix for phi focal (rad)
 * Input data: 155 data include 1.1, 1.5, 2.2 GeV
 * Output: the optimize coefficients for forward matrix
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
int const numpar=12;
double theta[dim],phi[dim],dp[dim],ytg[dim],yfoc[dim],phfoc[dim],errorz[dim],P[dim];

using namespace std;

Double_t func(float delta, float yt, float phit, float thetat, float pcent, Double_t *par)
{
  Double_t value= par[0]
    +par[1]*delta
    +par[2]*phit
    +par[3]*phit*delta
    +par[4]*phit*phit
    +par[5]*yt
    +par[6]*yt*delta
    +par[7]*thetat
    +par[8]*thetat*delta
    +par[9]*thetat*phit
    +par[10]*pcent
    +par[11]*pcent*thetat
    //    +par[12]*pcent*phit
   ;
 return value;
}

void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
   const Int_t nbins = dim;
   Int_t i;

//calculate chisquare
   Double_t chisq = 0;
   Double_t ddelta;
   for (i=0;i<nbins; i++) {
     ddelta  = (phfoc[i]-func(dp[i],ytg[i],phi[i],theta[i],P[i],par))/errorz[i];
     chisq += ddelta*ddelta;
     double tmp_test=func(dp[i],ytg[i],phi[i],theta[i],P[i],par);
     cout<<i<<"\t"<<phfoc[i]<<"\t"<<tmp_test<<"\t"<<phfoc[i]-tmp_test<<endl;
   }
   f = chisq;
}

void phi_focal()
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
       phfoc[count] = phfp;
       errorz[count]= err_phfp;
       P[count]     = P0;
       gr->SetPoint(count,yfp,phfp);
       count++;
     }

     //start TMinuit
     TCanvas * cc = new TCanvas("cc","",800,600);
     cc->Clear();
     gr->Draw("A*");
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
  0.042007 ,//offset
  0.249744 ,//dp
  3.25976  ,//phi
  -30.5807 ,//phi*dp
  -28.598  ,//phi^2
  1.36     ,//-2.19532,//y
  27.8299  ,//y*dp
  -1.1054  ,//theta
  3.35149  ,//theta*dp
  10.413   ,//theta*phi
  -0.001   ,//p0
  0.119    //,//p0*theta
  //  -0.004    //p0*phi
 };
   Double_t step[numpar] = 
{
  0.042007/10. ,//offset
  0.249744/10. ,//dp
  3.25976/10.  ,//phi
  30.5807/10.  ,//phi*dp
  28.598/10.   ,//phi^2
  2.19532/10.  ,//y
  27.8299/10.  ,//y*dp
  1.105482/10. ,//theta
  3.35149/10.  ,//theta*dp
  10.4130/10.  ,//theta*phi
  0.001        ,//p0
  0.001        //,//p0*theta
  //  0.001         //p0*phi
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
     cout<<vstart[ii]<<endl;
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
