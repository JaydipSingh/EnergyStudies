#include <iostream>
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

double UtilsPrint( double tl)
    {
      std::cout<<"Trajectory Track Length =>\t "<<tl<<"\t cm"<<"\n";
      return tl*100;
    }


bool StoppingMuon_track( bool E1, bool E2)

     {
        bool R=false;
         if((E1 || E2) && !(E1 && E2) ) 
         { R=true; std::cout<<"This is Stopping Muons Events "<<"\n";}
        //else
         // { R=false;
         // std::cout<<"This is Throughgoing Muons Events "<<"\n";
         //  }

       return R;

      }

void graphplot ( TGraph2D *gr)
      {
         TCanvas *tc1 = new TCanvas("tc1","Rec Space-Point",0,20,800,800);
            TH2D *h = new TH2D("h","Space",100,-720,720,100,-600,600);
             

             gStyle->SetOptStat(0);
              h->GetZaxis()->SetLimits(0.0,6500);
              gr->SetHistogram(h); 
              gr->SetMarkerColor(kRed);
              gr->SetMarkerStyle(20); 
              gr->SetMarkerSize(1.5);
         //  gr->Draw("pLINE same"); 
         //  gr->SetLineColor(kRed);
           gr->Draw("p0");
           gr->SetTitle("DUNE-FD ; X (cm); Y (cm); Z (cm)");
           gr->GetXaxis()->CenterTitle();gr->GetYaxis()->CenterTitle();gr->GetZaxis()->CenterTitle();
           gr->GetXaxis()->SetTitleOffset(1.5);gr->GetYaxis()->SetTitleOffset(1.5);gr->GetZaxis()->SetTitleOffset(1.3);

          tc1->SaveAs("/dune/app/users/jdsingh/DUNE_SP2020Work/DUNEWork/MUSUNWork/myplot/SpacePointPlot.pdf"); 
     } 
