/////////////////////////////////////////////////////////////////////////
////////////////////////////// plot3.C  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*                                                                     */
/*  PROGRAM:  plot3.C --- Plot the GAPP fit results for both the NP    */
/*                        and the SM best fit values. Take care of     */
/*                        of all 3 parameters.                         */
/*                                                                     */
/*  AUTHOR:   Kai Schmitz, HEPT Group, Michigan State University       */
/*                                                                     */
/*  DATE:     06/04/2009                                               */
/*                                                                     */
/*  SYNOPSIS:                                                          */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*  EXECUTION:                                                         */
/*                                                                     */
/*  === In the command line:                                           */
/*  ==                                                                 */
/*  ==  root -l -q -b 'plot3.C("filename","pltmd")'                    */
/*  ==                                                                 */
/*  === Defaults:                                                      */
/*  ==                                                                 */
/*  ==  filename ........... 'fp-d'                                    */
/*  ==  pltmd .............. 'phi'                                     */
/*  ==                                                                 */
/*  === Important examples:                                            */
/*  ==                                                                 */
/*  ==  root -l -q -b 'plotFit.C'                                      */
/*  ==  root -l -q -b 'plotFit.C("lr-t")'                              */
/*  ===                                                                */
/*                                                                     */
/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <math.h>
#include <cstdio>

//=======================================================================

// Main function.

void plot3(TString infile = "fp-d", TString pltmd = "cos") {

// CHECK FOR RIGHT INPUT ////////////////////////////////////////////////

   string strpltmd = pltmd, filename = infile, strfile = infile;

   if( (strpltmd.compare("cos") != 0 ) &&
       (strpltmd.compare("sin") != 0 ) &&
       (strpltmd.compare("tan") != 0 ) &&
       (strpltmd.compare("mmp") != 0 ) ) {error(4);};

// GLOBAL VARIABLES  ////////////////////////////////////////////////////

   Int_t file, point, color, style;

   Float_t fits2b, fittph, tphold, fitsph, fitcph, fitx, fitxmin, fitxmax = -1.0;

   Float_t xVal, yVal;

   Float_t xMin = 100000, xMax = -1.0, yMin = 100000, yMax = -1.0;

   Float_t MZ, MZmin, MW, MWmin, Mmin = 100000;

   Float_t Cz1, Cz2, Cz3, Cw1, Cw2, Cw3, Cw4, C1, C2;

   Float_t phiMin, phiMax, cphmin, cphmax, sphmin, sphmax;

// CUSTOMIZE PLOT ///////////////////////////////////////////////////////

   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetPalette(1);

   TCanvas *MyC = new TCanvas("MyC","Plot of the GAPP fit results",200,10,700,500);

   Float_t mmlegxmin, mmlegxmax, mmlegymin, mmlegymax;

   Float_t s2blegxmin, s2blegymin, s2blegxmax, s2blegymax;

   Float_t lblxmin, lblxmax, lblymin, lblymax;

   string plottitle = infile; //"Model: " + infile + "  |  Plot: ";

   sToUpper(plottitle);

   string xtitle, ytitle, NPleg, SMleg, display;

   NPleg = "#font[52]{M_{H}^{(NP)}, #bar{m}_{t}^{(NP)}}";
   SMleg = "#font[52]{M_{H}^{(SM)}, #bar{m}_{t}^{(SM)}}";

   if (strpltmd.compare("tan") == 0) { 

//    plottitle += "#font[42]{tan^{2}(#tilde{#phi}) over }#font[52]{#tilde{x}}#font[42]{.}";
      xtitle = "#font[52]{#tilde{x}}";
      ytitle = "#font[42]{tan^{2}(#tilde{#phi})}";
      display = "C";

      mmlegxmin = 0.15;
      mmlegxmax = 0.30;
      mmlegymin = 0.75;
      mmlegymax = 0.85;

      s2blegxmin = 0.15;
      s2blegxmax = 0.40;
      s2blegymin = 0.30;
      s2blegymax = 0.50;

      lblxmin = 0.88;
      lblxmax = 0.88;
      lblymin = 0.60;
      lblymax = 0.65;

   } else if (strpltmd.compare("cos") == 0) { 

//    plottitle += "#font[42]{cos(#tilde{#phi}) over }#font[52]{#tilde{x}}#font[42]{.}";
      xtitle = "#font[52]{#tilde{x}}";
      ytitle = "#font[42]{cos(#tilde{#phi})}";
      display = "C";

      mmlegxmin = 0.65;
      mmlegxmax = 0.80;
      mmlegymin = 0.65;
      mmlegymax = 0.75;

      s2blegxmin = 0.60;
      s2blegxmax = 0.85;
      s2blegymin = 0.30;
      s2blegymax = 0.50;

      lblxmin = 0.35;
      lblxmax = 0.50;
      lblymin = 0.45;
      lblymax = 0.50;

   } else if (strpltmd.compare("sin") == 0) {

//    plottitle += "#font[42]{sin(#tilde{#phi}) over }#font[52]{#tilde{x}}#font[42]{.}";
      xtitle = "#font[52]{#tilde{x}}";
      ytitle = "#font[42]{sin(#tilde{#phi})}";
      display = "C";

      mmlegxmin = 0.15;
      mmlegxmax = 0.30;
      mmlegymin = 0.75;
      mmlegymax = 0.85;

      s2blegxmin = 0.60;
      s2blegxmax = 0.85;
      s2blegymin = 0.30;
      s2blegymax = 0.50;

      lblxmin = 0.65;
      lblxmax = 0.80;
      lblymin = 0.60;
      lblymax = 0.65;

   } else if (strpltmd.compare("mmp") == 0) {

//    plottitle += "#font[42]{Masses of the new heavy gauge bosons.}";
      xtitle = "#font[52]{M_{Z'}}#font[42]{ (TeV)}";
      ytitle = "#font[52]{M_{W'}}#font[42]{ (TeV)}";
      display = "C";

      mmlegxmin = 0.20;
      mmlegxmax = 0.35;
      mmlegymin = 0.45;
      mmlegymax = 0.55;

      s2blegxmin = 0.20;
      s2blegxmax = 0.45;
      s2blegymin = 0.60;
      s2blegymax = 0.80;

      lblxmin = 0.65;
      lblxmax = 0.80;
      lblymin = 0.40;
      lblymax = 0.45;

   };
   
// PREPARE BOSON MASSES AND PHI BOUNDS //////////////////////////////////

   string mdl(filename,0,2);

   if ( (mdl.compare("lr") == 0) ||
        (mdl.compare("lp") == 0) ||
        (mdl.compare("hp") == 0) ||
        (mdl.compare("fp") == 0) ) {

      phiMin = 5.600; phiMax = 84.400;

      string Higgs(filename,3,1);

      if (Higgs.compare("d") == 0) {

         Cz1 = 11.95349795785275;
         Cz2 = 30.63269990028513;
         Cz3 = 42.58619785813789;
         Cw1 = 21.29309892906894;
         Cw2 = 9.339600971216193;
         Cw3 = 30.63269990028513;
         Cw4 = 42.58619785813789;

      }

      else if (Higgs.compare("t") == 0) {

         Cz1 = 5.976748978926375;
         Cz2 = 30.63269990028513;
         Cz3 = 85.17239571627579;
         Cw1 = 15.05649464522066;
         Cw2 = 3.302047590161717;
         Cw3 = 21.66058982554409;
         Cw4 = 60.22597858088265;

      }
   } 

   else if ( (mdl.compare("uu") == 0) ||
             (mdl.compare("nu") == 0) ) {

      phiMin = 10.179, phiMax = 79.821;

      C1 = 94.0397928463607;
      C2 = 77.1253849720165;

   } else {error(6);}

      cphmin = cos(TMath::Pi()*phiMin/180.0)*cos(TMath::Pi()*phiMin/180.0);
      cphmax = cos(TMath::Pi()*phiMax/180.0)*cos(TMath::Pi()*phiMax/180.0);
      sphmin = sin(TMath::Pi()*phiMin/180.0)*sin(TMath::Pi()*phiMin/180.0);
      sphmax = sin(TMath::Pi()*phiMax/180.0)*sin(TMath::Pi()*phiMax/180.0);

// LOOP OVER ROOT FILES  ////////////////////////////////////////////////
   
   for(file=0; file<=1; file++) {
   
      if(file==0) string epsfile  =  filename + "_" + strpltmd + ".eps";
      if(file==1) string filename =  filename + "_sm";
      string rootname = filename + ".root";
   
      TFile *rootfile = TFile::Open(rootname.c_str());      
      if(rootfile == NULL) error(1);
      
      TTree *tree = (TTree*)rootfile->Get(filename.c_str());
      if(tree == NULL) error(2);
   
      TBranch *fits2bbranch = (TBranch*)tree->GetBranch("fits2b");
      TBranch *fittphbranch = (TBranch*)tree->GetBranch("fittph");
      TBranch *fitxbranch   = (TBranch*)tree->GetBranch("fitx");

      if( (fits2bbranch == NULL) || 
          (fittphbranch == NULL) || 
          (fitxbranch   == NULL) ) error(3);

      tree->SetBranchAddress("fits2b",&fits2b);
      tree->SetBranchAddress("fittph",&fittph);
      tree->SetBranchAddress("fitx",  &fitx);
   
// GET ARRAYS ///////////////////////////////////////////////////////////

      Int_t Npoints = (Int_t)tree->GetEntries();

      Int_t tphStep = 0;

      Float_t tphMax = -1.0;

      for(point=0; point<Npoints; point++) {
      
         tree->GetEntry(point);
      
         if( fittph > tphMax ) {tphStep++; tphMax = fittph;}
      
      };

      const int tphSteps = tphStep;

      Float_t xArray[tphSteps], yArray[tphSteps], zArray[tphSteps];

      tphStep = -1, tphold = -1.0, fitxmin = 100000;

      for(point=0; point<Npoints; point++) {
      
         tree->GetEntry(point);

         if(fittph > tphold) {tphStep++; fitxmin = 100000;}

         fitsph = fittph / (1.0 + fittph);
         fitcph = 1.0 - fitsph; 
      
         if (strpltmd.compare("tan") == 0) {

            xVal = fitx;
            yVal = fittph;

         } else if (strpltmd.compare("cos") == 0) {

            xVal = fitx;
            yVal = sqrt(fitcph); 

         } else if (strpltmd.compare("sin") == 0) { 

            xVal = fitx; 
            yVal = sqrt(fitsph); 

         } else if (strpltmd.compare("mmp") == 0) { 

            if (fitsph != 0.0) {

               if ( (mdl.compare("lr") == 0) ||
                    (mdl.compare("lp") == 0) ||
                    (mdl.compare("hp") == 0) ||
                    (mdl.compare("fp") == 0) ) {

                  MMI(Cz1,Cz2,Cz3,Cw1,Cw2,Cw3,Cw4,fitx,fitsph,fitcph,fits2b,xVal,yVal);

               } else if ( (mdl.compare("uu") == 0) ||
                           (mdl.compare("nu") == 0) ) {

                  MMII(C1,C2,fitx,fitsph,fitcph,xVal,yVal);

               }
            }   
         }

         if( (strpltmd.compare("mmp") == 0) && (tphStep==1) ) { 
          
            xArray[0] = xArray[1]; 
            yArray[0] = yArray[1];
            zArray[0] = zArray[1]; 
 
         }

         if(fitx>fitxmax) fitxmax = fitx;

         if(fitx<fitxmin) {
 
            xArray[tphStep] = xVal;
            yArray[tphStep] = yVal;
            zArray[tphStep] = fits2b;
            fitxmin = fitx;

         }
  
         tphold = fittph;
 
      }

      if(file==0) TGraph *NPplot = new TGraph(tphSteps,xArray,yArray);
      if(file==1) TGraph *SMplot = new TGraph(tphSteps,xArray,yArray);

      TMarker *NPmrk[tphSteps], *SMmrk[tphSteps];

      MZmin = 100000, MWmin = 100000;

      for(tphStep=0; tphStep<tphSteps; tphStep++){

         marker(zArray[tphStep],color,style);

         if(file==0) { NPmrk[tphStep] = new TMarker(xArray[tphStep],yArray[tphStep],style);
                       NPmrk[tphStep]->SetMarkerSize(0.8);    
                       NPmrk[tphStep]->SetMarkerColor(color);}
         if(file==1) { SMmrk[tphStep] = new TMarker(xArray[tphStep],yArray[tphStep],style);
                       SMmrk[tphStep]->SetMarkerSize(0.8);    
                       SMmrk[tphStep]->SetMarkerColor(color);}

         if( (strpltmd.compare("mmp") != 0) || (tphStep !=0 )) {

            if (xArray[tphStep] < xMin) xMin = xArray[tphStep]; 
            if (xArray[tphStep] > xMax) xMax = xArray[tphStep]; 
            if (yArray[tphStep] < yMin) yMin = yArray[tphStep]; 
            if (yArray[tphStep] > yMax) yMax = yArray[tphStep]; 

         }

         if( ((strfile.compare("uu-d") != 0) || (strfile.compare("nu-d") != 0)) && (strpltmd.compare("cos") == 0) ) {

            fitx = xArray[tphStep];
            fitcph = yArray[tphStep]*yArray[tphStep];
            fitsph = 1.0 - fitcph;

            if( (cphmax < fitcph) && (fitcph < cphmin) ) {

               MMI(Cz1,Cz2,Cz3,Cw1,Cw2,Cw3,Cw4,fitx,fitsph,fitcph,fits2b,MZ,MW);
               if(MZ < MZmin) MZmin = MZ;
               if(MW < MWmin) MWmin = MW; 
            }
         }

         if( ((strfile.compare("uu-d") == 0) || (strfile.compare("nu-d") == 0)) && (strpltmd.compare("sin") == 0) ) {

            fitx = xArray[tphStep];
            fitsph = yArray[tphStep]*yArray[tphStep];
            fitcph = 1.0 - fitsph;

            if( (sphmin < fitsph) && (fitsph<sphmax) ) {

               MMII(C1,C2,fitx,fitsph,fitcph,MZ,MW);

               if(MZ < Mmin) { Mmin = MZ; cout << MZ << "\t" << sqrt(fitsph) << endl;}
            }
         }
      }

      cout << "(" << file << ") " << "MZmin: " << MZmin << endl;
      cout << "(" << file << ") " << "MWmin: " << MWmin << endl;
   }
 
// CREATE PLOTS /////////////////////////////////////////////////////////

   NPplot->SetLineStyle(2);
   NPplot->SetMarkerStyle(20);
   NPplot->SetMarkerSize(0.4);
   SMplot->SetMarkerStyle(20);
   SMplot->SetMarkerSize(0.4);
   
   if(strpltmd.compare("cos") == 0) {xMin = 0.0; yMin = 0.0; yMax = 1.0;}
   if(strpltmd.compare("sin") == 0) {yMin = 0.0; yMax = 1.0;}
   if(strpltmd.compare("mmp") == 0) {xMin = 0.0; xMax = 5.8/1.1; yMin = 0.0; yMax = 5;}

   TH1F* frame = MyC->DrawFrame(0.9*xMin,0.9*yMin,1.1*xMax,1.0*yMax);
   frame->SetTitle(plottitle.c_str());

   TAxis *xaxis = frame->GetXaxis();
   TAxis *yaxis = frame->GetYaxis();   
   xaxis->SetTitle(xtitle.c_str());
   xaxis->CenterTitle();
   xaxis->SetTitleOffset(1.);
   xaxis->SetDecimals();
   xaxis->SetLabelSize(0.03);
   xaxis->SetLabelOffset(0.01);
   yaxis->SetNdivisions(505);
   yaxis->SetTitle(ytitle.c_str());
   yaxis->CenterTitle();
   yaxis->SetTitleOffset(1.2);
   yaxis->SetDecimals();
   yaxis->SetLabelSize(0.03);
   yaxis->SetLabelOffset(0.01);

   TLegend *mmleg = new TLegend(mmlegxmin,mmlegymin,mmlegxmax,mmlegymax);
   mmleg->AddEntry(NPplot,NPleg.c_str(),"l");
   mmleg->AddEntry(SMplot,SMleg.c_str(),"l"); 
   mmleg->SetTextSize(0.025);
   mmleg->SetFillStyle(0);

   if( (strfile.compare("uu-d") != 0) && (strfile.compare("nu-d") != 0) ) {

      for(tphStep=0; tphStep<tphSteps; tphStep++){NPmrk[tphStep]->Draw(); SMmrk[tphStep]->Draw();}

   }

   Float_t xdummy[1] = {0.0}, ydummy[1] = {0.0};

   TGraph *circle = new TGraph(1,xdummy,ydummy);
   circle->SetMarkerStyle(24);
   circle->SetMarkerColor(kGreen+1);
   circle->SetMarkerSize(0.8); 
   
   TGraph *square = new TGraph(1,xdummy,ydummy);
   square->SetMarkerStyle(25);
   square->SetMarkerColor(kCyan+1);
   square->SetMarkerSize(0.8); 

   TGraph *triangle = new TGraph(1,xdummy,ydummy);
   triangle->SetMarkerStyle(26);
   triangle->SetMarkerColor(kBlue+1);
   triangle->SetMarkerSize(0.8); 

   TGraph *diamond = new TGraph(1,xdummy,ydummy);
   diamond->SetMarkerStyle(27);
   diamond->SetMarkerColor(kMagenta+1);
   diamond->SetMarkerSize(0.8); 

   TLegend *s2bleg = new TLegend(s2blegxmin,s2blegymin,s2blegxmax,s2blegymax); 

   s2bleg->AddEntry(circle,"#font[42]{0.00 < sin^{2}(2#tilde{#beta}) #leq 0.25}","p");
   s2bleg->AddEntry(square,"#font[42]{0.25 < sin^{2}(2#tilde{#beta}) #leq 0.50}","p");
   s2bleg->AddEntry(triangle,"#font[42]{0.50 < sin^{2}(2#tilde{#beta}) #leq 0.75}","p");
   s2bleg->AddEntry(diamond,"#font[42]{0.75 < sin^{2}(2#tilde{#beta}) #leq 1.00}","p");
   s2bleg->SetTextSize(0.025);
   s2bleg->SetFillStyle(0); 

   NPplot->Draw(display.c_str());
   SMplot->Draw(display.c_str());
   mmleg->Draw();
   if( (strfile.compare("uu-d") != 0) && (strfile.compare("nu-d") != 0) ) s2bleg->Draw();

// BOUNDS ON PHI //////////////////////////////////////////////////////

   Int_t i;

   const int iSteps = 100;    

   fitxmin = 1.0, fitxmax *= 1.5;

   Float_t deltax = (fitxmax-fitxmin)/iSteps;

   Float_t phixmin0[iSteps], phixmax0[iSteps], phixmin1[iSteps], phixmax1[iSteps];

   Float_t phiymin0[iSteps], phiymax0[iSteps], phiymin1[iSteps], phiymax1[iSteps];

   if ( (strpltmd.compare("tan") == 0) || 
        (strpltmd.compare("cos") == 0) ||
        (strpltmd.compare("sin") == 0) ) {

      for(i=0; i<100; i++) {

         fitx = fitxmin + i*deltax;

         phixmin0[i] = fitx;
         if (strpltmd.compare("tan") == 0) { phiymin0[i] = sphmin / cphmin; phiymax0[i] = sphmax / cphmax; }
         if (strpltmd.compare("cos") == 0) { phiymin0[i] = sqrt(cphmin);    phiymax0[i] = sqrt(cphmax); }
         if (strpltmd.compare("sin") == 0) { phiymin0[i] = sqrt(sphmin);    phiymax0[i] = sqrt(sphmax); }

      }

      TGraph *phiMin0 = new TGraph(iSteps,phixmin0,phiymin0);
      TGraph *phiMax0 = new TGraph(iSteps,phixmin0,phiymax0);

   } else if (strpltmd.compare("mmp") == 0) { 

      if ( (mdl.compare("lr") == 0) ||
           (mdl.compare("lp") == 0) ||
           (mdl.compare("hp") == 0) ||
           (mdl.compare("fp") == 0) ) {

         for(i=0; i<100; i++) {

            fitx = fitxmin + i*deltax;

            MMI(Cz1,Cz2,Cz3,Cw1,Cw2,Cw3,Cw4,fitx,sphmin,cphmin,0.0,phixmin0[i],phiymin0[i]);
            MMI(Cz1,Cz2,Cz3,Cw1,Cw2,Cw3,Cw4,fitx,sphmin,cphmin,1.0,phixmin1[i],phiymin1[i]);
            MMI(Cz1,Cz2,Cz3,Cw1,Cw2,Cw3,Cw4,fitx,sphmax,cphmax,0.0,phixmax0[i],phiymax0[i]);
            MMI(Cz1,Cz2,Cz3,Cw1,Cw2,Cw3,Cw4,fitx,sphmax,cphmax,1.0,phixmax1[i],phiymax1[i]);

         }

         TGraph *phiMin0 = new TGraph(iSteps,phixmin0,phiymin0);
         TGraph *phiMin1 = new TGraph(iSteps,phixmin1,phiymin1);
         TGraph *phiMax0 = new TGraph(iSteps,phixmax0,phiymax0);
         TGraph *phiMax1 = new TGraph(iSteps,phixmax1,phiymax1);

         phiMin1->SetLineStyle(7);
         phiMin1->SetMarkerStyle(22);
         phiMin1->SetMarkerSize(1.0);
         phiMax1->SetLineStyle(7);
         phiMax1->SetMarkerStyle(22);
         phiMax1->SetMarkerSize(1.0);

//       phiMin1->Draw("C");
//       phiMax1->Draw("C");

      } else if ( (mdl.compare("uu") == 0) ||
                  (mdl.compare("nu") == 0) ) {

         for(i=0; i<100; i++) {

            fitx = fitxmin + i*deltax;

            MMII(C1,C2,fitx,sphmin,cphmin,phixmin0[i],phiymin0[i]);
            MMII(C1,C2,fitx,sphmax,cphmax,phixmax0[i],phiymax0[i]);

         }

         TGraph *phiMin0 = new TGraph(iSteps,phixmin0,phiymin0);
         TGraph *phiMax0 = new TGraph(iSteps,phixmax0,phiymax0);

      }
   } 
      
   phiMin0->SetLineStyle(3);
   phiMin0->SetMarkerStyle(20);
   phiMin0->SetMarkerSize(0.4);
   phiMax0->SetLineStyle(3);
   phiMax0->SetMarkerStyle(20);
   phiMax0->SetMarkerSize(0.4);

   phiMin0->Draw("C");
   phiMax0->Draw("C");

// LABEL ALLOWED REGION ///////////////////////////////////////////////

   TPaveText *allowed = new TPaveText(lblxmin,lblymin,lblxmax,lblymax,"NDC");
   TText *text = allowed->AddText("#font[42]{allowed (95% CL)}");
   allowed->SetTextSize(0.04);
   if (strpltmd.compare("tan") == 0) text->SetTextAngle(270);
   allowed->SetFillStyle(0);
   allowed->SetLineColor(0);
   allowed->SetBorderSize(1);
   allowed->Draw();

// SAVE GRAPHIC ///////////////////////////////////////////////////////

   MyC->Print(epsfile.c_str());
  
}

//=======================================================================

// Calculate marker color and marker style.

void marker(Float_t fits2b, int& color, int& style) {

   if( (0.00<fits2b) && (fits2b<=0.25) ) {color = kGreen+1;   style = 24;}
   if( (0.25<fits2b) && (fits2b<=0.50) ) {color = kCyan+1;    style = 25;}
   if( (0.50<fits2b) && (fits2b<=0.75) ) {color = kBlue+1;    style = 26;}
   if( (0.75<fits2b) && (fits2b<=1.00) ) {color = kMagenta+1; style = 27;}

}

//=======================================================================

// Masses of the new heavy gauge bosons.

void MMI(Float_t Cz1, Float_t Cz2, Float_t Cz3, Float_t Cw1, Float_t Cw2, Float_t Cw3, Float_t Cw4, 
         Float_t x, Float_t sinphi, Float_t cosphi, Float_t sinbeta, Float_t& MZ, Float_t& MW) {

   MZ = (0.001/sqrt(sinphi*cosphi*x)) * (Cz1*cosphi*cosphi + Cz2*sinbeta + Cz3*x);
   MW = (0.001/sqrt(sinphi*x))  * (Cw1 - Cw2*cosphi*cosphi + Cw3*sinbeta + Cw4*x);

}

void MMII(Float_t C1, Float_t C2, Float_t x, Float_t sinphi, Float_t cosphi, Float_t& MZ, Float_t& MW) {

   MZ = (0.001/sqrt(sinphi*cosphi*x)) * (C1*sinphi*sinphi + C2*x);
   MW = (0.001/sqrt(sinphi*cosphi*x)) * (C1*sinphi*sinphi + C2*x);

}

//=======================================================================

// Convert strings to upper case.

void sToUpper(string& s)
{
   string::iterator i = s.begin();
   string::iterator end = s.end();

   while (i != end) {*i = toupper((unsigned char)*i); ++i;}
}

//=======================================================================

// Function to print out error messages.

void error(Int_t errorcode) {

   string errormessage = "\nERROR!\n";

   switch (errorcode){
      
      case 1: errormessage += "\nCannot open ROOT file(s)!"
                              "\nDefaults are 'fp-d.root' and 'fp-d_sm.root'."
                              "\nCheck FIRST argument in function call.";
      break;
     
      case 2: errormessage += "\nCannot open tree(s) in ROOT file(s)!"
                              "\nTree names are expected to correspond to the model names."
                              "\nCheck ROOT file or change code.";
      break;
      
      case 3: errormessage += "\nCannot open all required branches in the ROOT file(s)!"
                              "\nAll ROOT file must contain the branches"
                              "\n\n\t'fits2b, ''fittph' and 'fitx'.\n"
                              "\nCheck ROOT file(s) or change code.";
      break;
      
      case 4: errormessage += "\nIncorrect plot mode."
                              "\nPossible modes are:"
                              "\n\n\t'cos', 'sin', 'tan' and 'mmp'.\n"
                              "\nCheck SECOND argument in function call.";

      break;

      case 5: errormessage += "\nMasses of the new heavy gauge bosons"
                              "\nnot implemented for this model."
                              "\nAvailable classes of G(221) models are:"
                              "\n\n\t(BP-I,D), (BP-I,T) and (BP-II,D).\n"
                              "\nCheck (name of) the ROOT file.";

      break;

      case 6: errormessage += "\nHiggs representation not implemented for this model."
                              "\nAvailable classes of G(221) models are:"
                              "\n\n\t(BP-I,D), (BP-I,T) and (BP-II,D).\n"
                              "\nCheck (name of) the ROOT file.";

      break;

      default: errormessage += "\nError function called without specified errorcode."; 
      
   }

   errormessage += "\nAbort!\n";

   cerr << errormessage << endl; 
   
   abort();

}
