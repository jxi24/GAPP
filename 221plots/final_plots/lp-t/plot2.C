/////////////////////////////////////////////////////////////////////////
////////////////////////////// plot2.C  /////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*                                                                     */
/*  PROGRAM:  plot2.C --- Plot the GAPP fit results, two at time       */
/*                        (NP and SM).                                 */
/*                                                                     */
/*  AUTHOR:   Kai Schmitz, HEPT Group, Michigan State University       */
/*                                                                     */
/*  DATE:     06/02/2009                                               */
/*                                                                     */
/*  SYNOPSIS:                                                          */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*  EXECUTION:                                                         */
/*                                                                     */
/*  === In the command line:                                           */
/*  ==                                                                 */
/*  ==  root -l -q -b 'plot2.C("filename","pltmd")'                    */
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

void plot2(TString infile = "fp-d", TString pltmd = "cos") {

// CHECK FOR RIGHT INPUT ////////////////////////////////////////////////

   string strpltmd = pltmd, filename = infile;

   if( (strpltmd.compare("cos") != 0 ) &&
       (strpltmd.compare("sin") != 0 ) &&
       (strpltmd.compare("tan") != 0 ) &&
       (strpltmd.compare("mmp") != 0 ) ) {error(4);};

// GLOBAL VARIABLES  ////////////////////////////////////////////////////

   Int_t file, point;

   Float_t fits2b, fittph, fitsph, fitcph, fitx, xVal, yVal;

   Float_t xMin = 100000, xMax = -1.0, yMin = 100000, yMax = -1.0;

   Float_t Cz1, Cz2, Cz3, Cw1, Cw2, Cw3, Cw4, C1, C2;

// CUSTOMIZE PLOT ///////////////////////////////////////////////////////

   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetPalette(1);

   TCanvas *MyC = new TCanvas("MyC","Plot of the GAPP fit results",200,10,700,500);

   Float_t legxmin, legxmax, legymin, legymax;

   string plottitle = "Model: " + infile + "  |  Plot: ";

   string xtitle, ytitle, NPleg, SMleg, display;

   NPleg = "NP best fit values for M_{H} and #bar{m}_{t}.";
   SMleg = "SM best fit values for M_{H} and #bar{m}_{t}.";

   if (strpltmd.compare("tan") == 0) { 

      plottitle += "tan^{2}(#phi) over x.";
      xtitle = "x = u^{2}/v^{2}";
      ytitle = "tan^{2}(#phi)";
      display = "C";

      legxmin = 0.10;
      legxmax = 0.50;
      legymin = 0.80;
      legymax = 0.90;

   } else if (strpltmd.compare("cos") == 0) { 

      plottitle += "cos(#phi) over x.";
      xtitle = "x = u^{2}/v^{2}";
      ytitle = "cos(#phi)";
      display = "C";

      legxmin = 0.50;
      legxmax = 0.90;
      legymin = 0.80;
      legymax = 0.90;

   } else if (strpltmd.compare("sin") == 0) {

      plottitle += "sin(#phi) over x.";
      xtitle = "x = u^{2}/v^{2}";
      ytitle = "sin(#phi)";
      display = "C";

      legxmin = 0.50;
      legxmax = 0.90;
      legymin = 0.15;
      legymax = 0.25;

   } else if (strpltmd.compare("mmp") == 0) {

      plottitle += "Masses of the new heavy gauge bosons";
      xtitle = "M_{Z'} (TeV)";
      ytitle = "M_{W'} (TeV)";
      display = "C";

      legxmin = 0.10;
      legxmax = 0.60;
      legymin = 0.80;
      legymax = 0.90;

   };
   
// PREPARE BOSON MASSES /////////////////////////////////////////////////

   string mdl(filename,0,2);

   if ( (mdl.compare("lr") == 0) ||
        (mdl.compare("lp") == 0) ||
        (mdl.compare("hp") == 0) ||
        (mdl.compare("fp") == 0) ) {

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

      C1 = 94.0397928463607;
      C2 = 77.1253849720165;

   } else {error(6);}

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

      Int_t Nentries = (Int_t)tree->GetEntries();

      const int Npoints = Nentries;

      Float_t xArray[Npoints], yArray[Npoints];

      for(point=0; point<Npoints; point++) {
      
         tree->GetEntry(point);

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

                  xVal = (0.001/sqrt(fitsph*fitcph*fitx)) * (Cz1*fitcph*fitcph + Cz2*fits2b + Cz3*fitx);
                  yVal = (0.001/sqrt(fitsph*fitx))  * (Cw1 - Cw2*fitcph*fitcph + Cw3*fits2b + Cw4*fitx);

               } else if ( (mdl.compare("uu") == 0) ||
                           (mdl.compare("nu") == 0) ) {

                  xVal = (0.001/sqrt(fitsph*fitcph*fitx)) * (C1*fitsph*fitsph + C2*fitx);
                  yVal = (0.001/sqrt(fitsph*fitcph*fitx)) * (C1*fitsph*fitsph + C2*fitx);
            
               }
            } 

            if(point==1) { xArray[0] = xVal; yArray[0] = yVal; }
         }

         xArray[point] = xVal;
         yArray[point] = yVal;
  
         if(point>0) {

            if (xVal < xMin) xMin = xVal; if (xVal > xMax) xMax = xVal; 
            if (yVal < yMin) yMin = yVal; if (yVal > yMax) yMax = yVal; 
       
         }
      }

      if(file==0) TGraph *NPplot = new TGraph(Npoints,xArray,yArray);
      if(file==1) TGraph *SMplot = new TGraph(Npoints,xArray,yArray);

   }
 
// CREATE PLOTS /////////////////////////////////////////////////////////

   NPplot->SetLineStyle(2);
   NPplot->SetMarkerStyle(20);
   NPplot->SetMarkerSize(0.4);
   SMplot->SetMarkerStyle(20);
   SMplot->SetMarkerSize(0.4);

   TH1F* frame = MyC->DrawFrame(0.9*xMin,0.9*yMin,1.1*xMax,1.1*yMax);
   frame->SetTitle(plottitle.c_str());

   TAxis *xaxis = frame->GetXaxis();
   TAxis *yaxis = frame->GetYaxis();   
   xaxis->SetTitle(xtitle.c_str());
   xaxis->CenterTitle();
   xaxis->SetTitleOffset(1.);
   xaxis->SetDecimals();
   xaxis->SetLabelSize(0.03);
   xaxis->SetLabelOffset(0.01);
   yaxis->SetTitle(ytitle.c_str());
   yaxis->CenterTitle();
   yaxis->SetTitleOffset(1.);
   yaxis->SetDecimals();
   yaxis->SetLabelSize(0.03);
   yaxis->SetLabelOffset(0.01);

   TLegend *leg = new TLegend(legxmin,legymin,legxmax,legymax);
   leg->AddEntry(NPplot,NPleg.c_str(),"l");
   leg->AddEntry(SMplot,SMleg.c_str(),"l");
   leg->SetTextSize(0.03);
   leg->SetFillColor(0);

   NPplot->Draw(display.c_str());
   SMplot->Draw(display.c_str());
   leg->Draw();

   MyC->Print(epsfile.c_str());
  
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
