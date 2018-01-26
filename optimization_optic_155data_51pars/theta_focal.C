/*===============================================================
 * Nguyen Ton Sep 12th 2016
 * Jan 2018: update on optic data, include 2.2 GeV
 * Purpose: optimize forward matrix for yfoc
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
int const numpar=10;
Float_t theta[dim],phi[dim],dp[dim],ytg[dim],xfoc[dim],thfoc[dim],errorz[dim],P[dim];

using namespace std;

Double_t thfunc(float delta, float yt, float phit, float thetat, float pcent, Double_t *par)
{
  Double_t value= par[0]
    +par[1]*delta
    +par[2]*phit
    +par[3]*yt
    +par[4]*thetat
    +par[5]*thetat*yt
    +par[6]*thetat*phit
    +par[7]*pcent
    +par[8]*pcent*thetat
    +par[9]*pcent*phit
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
     ddelta  = (thfoc[i]-thfunc(dp[i],ytg[i],phi[i],theta[i],P[i],par))/errorz[i];
     chisq+=ddelta*ddelta;
     double th_test=thfunc(dp[i],ytg[i],phi[i],theta[i],P[i],par);
     cout<<i<<"\t"<<thfoc[i]<<"\t"<<th_test<<"\t"<<thfoc[i]-th_test<<endl;
   }
   f = chisq;
}

void theta_focal()
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
       errorz[count]= err_thfp;
       P[count] = P0;
       gr->SetPoint(count,xfp,thfp);
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
   //double tmp_init = -1000 + 20*k;
   Double_t vstart[numpar] = 
{
  0.046127   ,//offset
  -0.039     , //dp
  1.25       ,//phi
  // -7.254,//phi^2
  -0.464     ,//y
  //-500.,//phi*y
  -0.039     ,//theta
  -30.       ,//theta*y
  -30.       ,//theta*phi
  -0.01      ,//p0
  0.02       ,//p0*theta
  -0.05       //p0*phi
 };
   Double_t step[numpar] = 
{
  0.0543/10.  ,//offset
  0.039/10.   ,//dp
  1.25/10.    ,//phi
  //7.254/10.,//phi^2
  0.464/10.   ,//y
  //50./10.,//phi*y
  0.039/10.   ,//theta
  50./10.     ,//theta*y
  30./10.     ,//theta*phi
  0.001       ,//p0
  0.001       ,//p0*theta
  0.001        //p0*phi
};



   for(int ii=0;ii<numpar;ii++)
     {
       gMinuit->mnparm(ii, "a", vstart[ii], step[ii], 0,0,ierflg);
     }

// Now ready for minimization step
   arglist[0] = 500000;
   arglist[1] = 1;

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
     //gr->SetPoint(k,tmp_init,amin);
     //}
     //gr->Draw("A*");
}
