//////////////////////////////////////////////////////////////////////
//////////////////////////   out2root.C   ////////////////////////////
//////////////////////////////////////////////////////////////////////
/*                                                                  */
/*  COMMENT:                                                        */
/*                                                                  */
/*  Conversion of .out to .root files.                              */
/*                                                                  */
/*  EXECUTION:                                                      */
/*                                                                  */
/*  === In the command line:                                        */
/*  ==                                                              */
/*  ==  root -l -q -b 'out2root.C("filename")'                      */
/*  ==                                                              */
/*  ==  Defaults:                                                   */
/*  ==                                                              */
/*  ==  filename .............. fp-d                                */
/*  ==                                                              */
/*  ==  Example:                                                    */
/*  ==                                                              */
/*  ==  root -l -q -b 'out2root.C("lr-t")'                          */
/*                                                                  */
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <cstdio>

// Main function that converts .out to .root files.

void out2root(TString filename="fp-d") {

//####################################################################
// Open the .out file:

   string strfile = filename;

   string outfile = "plot_" + strfile + ".out";

   ifstream infile;

   infile.open(outfile.c_str(), ifstream::in);

   // Check whether the file is available.

   if(infile.fail()) error(1);

//####################################################################
// Create the output .root file:

   string rootname;

   rootname = strfile + ".root";

   TFile *rootfile = new TFile(rootname.c_str(),"RECREATE");

//####################################################################
// Declaration of important variables:

   Float_t fits2b, fittph, fitx;

   string line,junk;

   Int_t headings = 0;

//####################################################################
// Generate the ROOT tree:

  TTree *tree = new TTree(strfile.c_str(),"Contours From GAPP Fit");

  tree->Branch("fits2b",&fits2b,"fits2b");
  tree->Branch("fittph",&fittph,"fittph");
  tree->Branch("fitx",  &fitx,  "fitx");

//####################################################################
// Reading the input file line by line and writing the root file:

   cout << "Converting plot_" << strfile << ".out to " << strfile << ".root!" << endl;

   while ( getline(infile,line) ) {

      if( headings != 0) {

            stringstream buffer(line);

            buffer >> fits2b >> fittph >> fitx;

            tree->Fill();

      }

      headings += 1;

   }

   rootfile->Write("",TObject::kOverwrite);

   rootfile->Close();
}

//=======================================================================

// Function to print out error messages.

void error(Int_t errorcode) {

   string errormessage = "\nERROR!\n";

   switch (errorcode){
      
      case 1: errormessage += "\nCannot open .out file."
                              "\nDefault is 'plot_fp-d.out'."
                              "\nCheck first argument in function call."
                              "\nType only the model name (e.g., 'fp-d').";
      break;
     
      default: errormessage += "\nError function called without specified errorcode."; 
      
   }

   errormessage += "\nAbort!\n";

   cerr << errormessage << endl; 
   
   abort();

}
