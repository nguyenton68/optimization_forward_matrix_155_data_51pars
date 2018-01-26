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
int const numpar=15;
Float_t theta[dim],phi[dim],dp[dim],ytg[dim],yfoc[dim],phfoc[dim],errorz[dim],P[dim], I[dim];

using namespace std;

//====================
Double_t yfunc(float delta, float yt, float phit, float thetat, float pcent,Double_t *par)
{
  Double_t value= par[0]
    +par[1]*delta
    +par[2]*phit
    +par[3]*phit*delta
    +par[4]*phit*phit
    +par[5]*yt
    +par[6]*yt*delta
    +par[7]*phit*yt
    +par[8]*thetat
    +par[9]*thetat*thetat
    +par[10]*thetat*delta
    +par[11]*thetat*phit
    +par[12]*pcent
    +par[13]*pcent*thetat
    +par[14]*pcent*phit
    // +par[15]*pcent*pcent
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
     ddelta  = (yfoc[i]-yfunc(dp[i],ytg[i],phi[i],theta[i],P[i],par))/errorz[i];
     chisq+=ddelta*ddelta;
     double y_test=yfunc(dp[i],ytg[i],phi[i],theta[i],P[i],par);
     //cout<<i<<"\t"<<yfoc[i]<<"\t"<<y_test<<"\t"<<ddelta<<"\t"<<chisq<<endl;
     cout<<i+4<<"\t"<<yfoc[i]<<"\t"<<y_test<<"\t"<<yfoc[i]-y_test<<endl;
   }
   f = chisq;
   cout<<"Final chi square= "<<f<<endl;
}

void y_focal()
{
     ifstream data("full_kin_155_pnts.dat",ios_base::in);
     double th_tg,ph_tg,y_tg,dpp,yfp,phfp,xfp,thfp,err_yfp,err_xfp,err_phfp,err_thfp,P0,cur;
     int count=0;

     TGraph * gr0 = new TGraph();
     gr0->SetMarkerColor(kBlack);
     gr0->SetMarkerSize(1.0);
     gr0->SetTitle("Combine all holes in Focal plane;y_{fp}(m);#phi_{fp}(rad)");
     TGraph *gr1 = new TGraph();
     gr1->SetMarkerColor(kRed);
     gr1->SetMarkerSize(1.0);
     TMultiGraph *mgr = new TMultiGraph();
     mgr->Add(gr0);
     mgr->Add(gr1);
     TH1F * herr = new TH1F("herr","",200,0,0.01);
     while(!data.eof()){
       data>>ph_tg>>th_tg>>y_tg>>dpp>>xfp>>err_xfp>>yfp>>err_yfp>>phfp>>err_phfp>>thfp>>err_thfp>>P0;
       theta[count] = th_tg;
       phi[count]   = ph_tg;
       ytg[count]   = y_tg;
       dp[count]    = dpp;
       yfoc[count]  = yfp;
       phfoc[count] = phfp;
       errorz[count]=err_yfp;
       P[count] = P0;
       I[count] = cur;
       if(count<97){
       gr0->SetPoint(count,yfp,phfp);
       } else {
       gr1->SetPoint(count-97,yfp,phfp);
       }
       herr->Fill(err_yfp);
       count++;
     }

     //start TMinuit
     TCanvas * cc = new TCanvas("cc","",800,600);
     cc->Clear();
     mgr->Draw("A*");
     //herr->Draw();
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
  -0.035654,//offset
  0.06281,//dp
  0.376, //phi
  -20.7935, //phi*dp
  24.69,//phi^2
  0.2027,//y
  -12.0805, //y*dp
  -50.,//phi*y
  0.905,//theta
  0.08,//theta*theta
  -20.,//theta*dp
  -20., //thete*phi
  -0.05, //p0
  -0.1, //p0*theta
  -0.05//, // p0*phi
  //-0.001
 };
   Double_t step[numpar] = 
{
  0.033923/10.    ,//offset
  0.064/10.       ,//dp
  0.639284/10.    , //phi
  20.7935/10.     , //phi*dp
  21.8367/10.     ,//phi^2
  0.022990/10.    ,//y
  12.0805/10.     , //y*dp
  50./10.         ,//phi*y
  0.8626/10.      ,//theta
  10.             ,//theta*theta
  1.87960/10.     ,//theta*dp
  2.              ,//theta*phi
  0.001           ,//p0
  0.001           ,//p0*theta
  0.001//,           //p0*phi
  //0.0001
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
