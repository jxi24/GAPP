/////////////////////////////////////////////////////////////////////////
///////////////////////////// plotFit.C  ////////////////////////////////
/////////////////////////////////////////////////////////////////////////
/*                                                                     */
/*  PROGRAM:  plotFit.C --- Plot the GAPP fit results                  */
/*                                                                     */
/*  AUTHOR:   Kai Schmitz, HEPT Group, Michigan State University       */
/*                                                                     */
/*  DATE:     04/30/2009                                               */
/*                                                                     */
/*  SYNOPSIS:                                                          */
/*                                                                     */
/*                                                                     */
/*                                                                     */
/*  EXECUTION:                                                         */
/*                                                                     */
/*  === In the command line:                                           */
/*  ==                                                                 */
/*  ==  root -l -q -b 'plotFit.C("filename","pltmd")'                  */
/*  ==                                                                 */
/*  === Defaults:                                                      */
/*  ==                                                                 */
/*  ==  filename ........... 'fp-d'                                    */
/*  ==  pltmd .............. 'tph'                                     */
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

void plotFit(TString filename = "fp-d", TString pltmd = "tph") {

// CHECK FOR RIGHT INPUT ////////////////////////////////////////////////

   string strpltmd = pltmd;

   if( strpltmd.compare("tph")  != 0 &&
       strpltmd.compare("s2b")  != 0 &&
       strpltmd.compare("mmp") != 0 ) {error(4);};
   
// OPEN THE ROOT FILE  //////////////////////////////////////////////////

   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetTitleBorderSize(0);
   gStyle->SetPalette(1);

   TCanvas *MyC = new TCanvas("MyC","Plot of the GAPP fit",200,10,700,500);

// Still to do: Automate the frame boundaries. 
   
   string strfile = filename, rootname = strfile + ".root";

   TFile *rootfile = TFile::Open(rootname.c_str());
      
   if(rootfile == NULL) error(1);
      
   TTree *tree = (TTree*)rootfile->Get(strfile.c_str());
   
   if(tree == NULL) error(2);
   
   TBranch *fits2bbranch = (TBranch*)tree->GetBranch("fits2b");
   
   TBranch *fittphbranch = (TBranch*)tree->GetBranch("fittph");
   
   TBranch *fitxbranch   = (TBranch*)tree->GetBranch("fitx");
   
   if( (fits2bbranch == NULL) || (fittphbranch == NULL) || (fitxbranch == NULL) ) error(3);
   
   Float_t fits2b, fittph, fitx;
      
   tree->SetBranchAddress("fits2b",&fits2b);

   tree->SetBranchAddress("fittph",&fittph);
   
   tree->SetBranchAddress("fitx",  &fitx);
   
   
// GET GRID /////////////////////////////////////////////////////////////

   Int_t Npoints = (Int_t)tree->GetEntries();

   Int_t point, zSteps = 0, ySteps = 0;

   Float_t fitxMin = 100000, fitxMax = -1.0, s2bMin = 100000, tphMin = 100000;

   Float_t s2bMax = -1.0, tphMax = -1.0;

   for(point=0; point<Npoints; point++) {
      
      tree->GetEntry(point);
      
      if( fits2b > s2bMax ) {zSteps++; s2bMax = fits2b;}

      if( fittph > tphMax ) {ySteps++; tphMax = fittph;}
      
   };

   const int s2bSteps = zSteps, tphSteps = ySteps, mmpSteps = Npoints;

   Float_t s2bValues[s2bSteps], tphValues[tphSteps];

   s2bMax = -1.0, tphMax = -1.0;

   int s2bStep = 0, tphStep = 0;

   for(point=0; point<Npoints; point++) {
      
      tree->GetEntry(point);
      
      if( fits2b > s2bMax ) {s2bValues[s2bStep] = fits2b; s2bStep++; s2bMax = fits2b;}

      if( fittph > tphMax ) {tphValues[tphStep] = fittph; tphStep++; tphMax = fittph;}
      
   };

// PREPARE PLOT /////////////////////////////////////////////////////////

   string plottitle = "Model: " + strfile + "  |  Plot: ";

   if( strpltmd.compare("tph") == 0 ) {

      plottitle += "tan^{2}(#phi) over x for fixed sin^{2}(2#beta)";

   }

   if( strpltmd.compare("s2b") == 0 ) {

      plottitle += "sin^{2}(2#beta) over x for fixed tan^{2}(#phi)";

   }

   if( strpltmd.compare("mmp") == 0 ) {

      plottitle += "Masses of the new heavy gauge bosons";

   }

//  PLOT DATA ///////////////////////////////////////////////////////////

   if (strpltmd.compare("tph") == 0) {

      tphMax = -1.0;

//    TGraph *tphplots[s2bSteps];
      TGraph *tphplots[s2bSteps-30];

//    for(s2bStep=0; s2bStep<s2bSteps; s2bStep++) {
      for(s2bStep=0; s2bStep<s2bSteps-30; s2bStep++) {

         Float_t tphArray[tphSteps], fitxArray[tphSteps];

         tphStep = 0;

         for(point=0; point<Npoints; point++) {
      
            tree->GetEntry(point);
      
            if(fits2b == s2bValues[s2bStep]) {

//             tphArray[tphStep]  = fittph;
               tphArray[tphStep]  = sqrt(1.0/(1.0+fittph));

               if (fittph < tphMin) tphMin = fittph;

               if (fittph > tphMax) tphMax = fittph;

               fitxArray[tphStep] = fitx;

               if (fitx < fitxMin) fitxMin = fitx;

               if (fitx > fitxMax) fitxMax = fitx;

               tphStep++;


              TMarker *m = new TMarker(fitxArray[tphStep],tphArray[tphStep],20);
              m->SetMarkerSize(2);
              m->SetMarkerColor(31+tphStep);
              m->Draw();



            }
         }

         if (s2bStep == 0) {

            TH1F* frame = MyC->DrawFrame(0.0,0.0,1.1*fitxMax,1.1);
//          TH1F* frame = MyC->DrawFrame(0.7*fitxMin,0.7*tphMin,1.1*fitxMax,1.1*tphMax);

            TAxis *xaxis = frame->GetXaxis();
            TAxis *yaxis = frame->GetYaxis();
   
            xaxis->SetTitle("x = u^{2}/v^{2}");
            xaxis->CenterTitle();
            xaxis->SetTitleOffset(1.);
            xaxis->SetDecimals();
            xaxis->SetLabelSize(0.03);
            xaxis->SetLabelOffset(0.01);
        
            yaxis->SetTitle("tan^{2}(#phi)");
            yaxis->CenterTitle();
            yaxis->SetTitleOffset(1.);
            yaxis->SetDecimals();
            yaxis->SetLabelSize(0.03);
            yaxis->SetLabelOffset(0.01);

            frame->SetTitle(plottitle.c_str());

         }

         tphplots[s2bStep] = new TGraph(tphSteps,fitxArray,tphArray);

         tphplots[s2bStep]->SetMarkerStyle(20);

         tphplots[s2bStep]->SetMarkerSize(0.4);
   
         tphplots[s2bStep]->Draw("CP");

      }
   }

   else if (strpltmd.compare("s2b") == 0) {

      s2bMax = -1.0;

      TGraph *s2bplots[tphSteps-100];

      for(tphStep=0; tphStep<tphSteps-100; tphStep++) {

         Float_t s2bArray[s2bSteps], fitxArray[s2bSteps];

         s2bStep = 0;

         for(point=0; point<Npoints; point++) {
      
            tree->GetEntry(point);
      
            if(fittph == tphValues[tphStep+20]) {

               s2bArray[s2bStep]  = fits2b;

               if (fits2b < s2bMin) s2bMin = fits2b;

               if (fits2b > s2bMax) s2bMax = fits2b;

               fitxArray[s2bStep] = fitx;

               if (fitx < fitxMin) fitxMin = fitx;

               if (fitx > fitxMax) fitxMax = fitx;

               s2bStep++;              

            }
         }

         if (tphStep == 0) {

            TH1F* frame = MyC->DrawFrame(0.8*fitxMin,0.95*s2bMin,1.2*fitxMax,1.05*s2bMax);

            TAxis *xaxis = frame->GetXaxis();
            TAxis *yaxis = frame->GetYaxis();
   
            xaxis->SetTitle("x = u^{2}/v^{2}");
            xaxis->CenterTitle();
            xaxis->SetTitleOffset(1.);
            xaxis->SetDecimals();
            xaxis->SetLabelSize(0.03);
            xaxis->SetLabelOffset(0.01);

            yaxis->SetTitle("sin^{2}(2#beta)");
            yaxis->CenterTitle();
            yaxis->SetTitleOffset(1.25);
            yaxis->SetDecimals();
            yaxis->SetLabelSize(0.03);
            yaxis->SetLabelOffset(0.01);

            frame->SetTitle(plottitle.c_str());

         }

         s2bplots[tphStep] = new TGraph(s2bSteps,fitxArray,s2bArray);

         s2bplots[tphStep]->SetMarkerStyle(20);

         s2bplots[tphStep]->SetMarkerSize(0.4);

         s2bplots[tphStep]->Draw("C");

      }
   }

   else if (strpltmd.compare("mmp") == 0) {

      Float_t mzpArray[mmpSteps], mwpArray[mmpSteps];

      Float_t mzpMin = 100000, mzpMax = -1.0, mwpMin = 100000, mwpMax = -1.0;

      Float_t fitsph, fitcph; 

      Float_t Cz1, Cz2, Cz3, Cw1, Cw2, Cw3, Cw4, C1, C2;

      string mdl(strfile,0,2);

      if ( (mdl.compare("lr") == 0) ||
           (mdl.compare("lp") == 0) ||
           (mdl.compare("hp") == 0) ||
           (mdl.compare("fp") == 0) ) {

         string Higgs(strfile,3,1);

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

         } else {error(6);}

         for(point=0; point<Npoints; point++) {
      
            tree->GetEntry(point);

            fitsph = fittph / (1.0 + fittph);

            fitcph = 1.0 - fitsph;

            if (fitsph != 0.0) {

               mzpArray[point] = (0.001/sqrt(fitsph*fitcph*fitx)) * (Cz1*fitcph*fitcph + Cz2*fits2b + Cz3*fitx);

               if (mzpArray[point] < mzpMin) mzpMin = mzpArray[point];

               if (mzpArray[point] > mzpMax) mzpMax = mzpArray[point];

               mwpArray[point] = (0.001/sqrt(fitsph*fitx))  * (Cw1 - Cw2*fitcph*fitcph + Cw3*fits2b + Cw4*fitx);

               if (mwpArray[point] < mwpMin) mwpMin = mwpArray[point];

               if (mwpArray[point] > mwpMax) mwpMax = mwpArray[point];

            } else {

               mzpArray[point] = 0.0;

               mwpArray[point] = 0.0;

            }
         }
      } 

      else if ( (mdl.compare("uu") == 0) ||
                (mdl.compare("nu") == 0) ) {

         C1 = 94.0397928463607
         C2 = 77.1253849720165

         for(point=0; point<Npoints; point++) {
      
            tree->GetEntry(point);

            fitsph = fittph / (1.0+fittph);

            fitcph = 1.0 - fitsph;

            if (fitsph != 0.0) {

               mzpArray[point] = (0.001/sqrt(fitsph*fitcph*fitx)) * (C1*fitsph*fitsph + C2*fitx);

               if (mzpArray[point] < mzpMin) mzpMin = mzpArray[point];

               if (mzpArray[point] > mzpMax) mzpMax = mzpArray[point];

               mwpArray[point] = (0.001/sqrt(fitsph*fitcph*fitx)) * (C1*fitsph*fitsph + C2*fitx);

               if (mwpArray[point] < mwpMin) mwpMin = mwpArray[point];

               if (mwpArray[point] > mwpMax) mwpMax = mwpArray[point];

            } else {

               mzpArray[point] = 0.0;

               mwpArray[point] = 0.0;

            }
         }
      }

      else {error(5);}

      TH1F* frame = MyC->DrawFrame(0.7*mzpMin,0.7*mwpMin,1.1*mzpMax,1.1*mwpMax);

      TAxis *xaxis = frame->GetXaxis();
      TAxis *yaxis = frame->GetYaxis();
   
      xaxis->SetTitle("M_{Z'} (TeV)");
      xaxis->CenterTitle();
      xaxis->SetTitleOffset(1.);
      xaxis->SetDecimals();
      xaxis->SetLabelSize(0.03);
      xaxis->SetLabelOffset(0.01);

      yaxis->SetTitle("M_{W'} (TeV)");
      yaxis->CenterTitle();
      yaxis->SetTitleOffset(1.25);
      yaxis->SetDecimals();
      yaxis->SetLabelSize(0.03);
      yaxis->SetLabelOffset(0.01);

      frame->SetTitle(plottitle.c_str());

      TGraph *mmpplot = new TGraph(Npoints,mzpArray,mwpArray);

      mmpplot->SetMarkerStyle(20);

      mmpplot->SetMarkerSize(0.3);

      mmpplot->Draw("P");

   }
  
   else error(4);

//  SAVE PLOT IN .EPS FILE /////////////////////////////////////////////?
  
   string epsfile = strfile + "_" + strpltmd + ".eps";

   MyC->Print(epsfile.c_str());
  
}

//=======================================================================

// Function to print out error messages.

void error(Int_t errorcode) {

   string errormessage = "\nERROR!\n";

   switch (errorcode){
      
      case 1: errormessage += "\nCannot open ROOT file!"
                              "\nDefault is 'fp-d.root'."
                              "\nCheck FIRST argument in function call.";
      break;
     
      case 2: errormessage += "\nCannot open tree in ROOT file!"
                              "\nThe tree name is expected to correspond to the model name."
                              "\nCheck ROOT file or change code.";
      break;
      
      case 3: errormessage += "\nCannot open all required branches in ROOT file!"
                              "\nROOT file must contain the branches"
                              "\n\n\t'fits2b', 'fittph' and 'fitx'.\n"
                              "\nCheck ROOT file or change code.";
      break;
      
      case 4: errormessage += "\nIncorrect plot mode."
                              "\nPossible modes are:"
                              "\n\n\t'tph', 's2b' and 'mzmw'.\n"
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
