/**********************************************************************
 * This is a test file to call  fcn and chi2 fortran files 
 * Kirtimaan 11/10/2014
 * ********************************************************************/
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<vector>
#include<algorithm>
#include<iterator>
#include<ctype.h>

#include <boost/algorithm/string.hpp>
// Root includes
#include "TMinuit.h"
#include "TGraph.h"
using namespace std;
//Defintion of subroutine fcn function in F/core/chi2.f
extern "C"{
	void fcn_(int* npar, double* grad, double* fval, double* xval, int* iflag, double (*chi2)(double *, int * , double * ,double* ) );
	double chi2_(double* xval, int* npar, double* smval, double* pull);
}
//************************* fcn wrapper ************************************
void cpp_fcn( int &ntpar, double *grad, double &fval, double *xval, int iflag)
{
	
	fcn_( &ntpar, grad, &fval, xval ,&iflag,chi2_);
	//cout<<"ntpar= "<<ntpar<<" iflag= "<<iflag<<endl;
	//cout<<"fval= "<< fval <<endl;
	//for (int i=0; i< ntpar;i++){cout<<"xval["<<i<<"]= "<<xval[i]<<endl;}
}
//************************ Common block variables***************************
extern "C" {
               extern struct{
		       //*** Logical variables defined in common block of fortran common.f
		       /*
		        *       common /flags/   flagmr,flgblm,ffermi,fa2mt4,fa2mt2,fa2mt0,fla2im,
     .                 fasmt2,fasmt0,fas2mt,flagmf,fpolew,flgech,fobliq,
     .                 f4lqcd,falas2,fbayes,flagzp,flagmh,flagmt,flagmc,
     .                 flgrho,fkappa,flagS,flagT,fzprim,fsinth,fwrite,
     .                 flprob,fhiggs,flagal,fsplot
     *  flgblm and flagzp are of integer type in fortran rest are logical
     * */
                    int flagmr,flgblm,ffermi,fa2mt4,fa2mt2,fa2mt0,fla2im,
                      fasmt2,fasmt0,fas2mt,flagmf,fpolew,flgech,fobliq,
                      f4lqcd,falas2,fbayes,flagzp,flagmh,flagmt,flagmc,
                      flgrho,fkappa,flagS,flagT,fzprim,fsinth,fwrite,
                      flprob,fhiggs,flagal,fsplot;
               } flags_;
               /*      common /fitflg/  flgfitx, flgtph, flgs2b, flfout*/
               extern struct{
		       int flgfitx, flgtph, flgs2b, flfout;
	       }fitflg_;
             }
//************************ Initialize logical flags ************************
void Initlog(){
	//*** 1 =True *** 0 =False

        fitflg_.flgfitx = 1;
        fitflg_.flgtph  = 1;
        fitflg_.flgs2b  = 1;

        fitflg_.flfout = 1;

        flags_.fwrite = 0;
// ----------------------------------------

        flags_.fsplot = 1;
        flags_.flprob = 1;
// ----------------------------------------
        flags_.flagmh = 1;
// ----------------------------------------
        flags_.flagmt = 1;
        flags_.flagmc = 1;
        flags_.flagal = 1;
        flags_.flagS  = 1;
        flags_.flagT  = 1;
        flags_.flgrho = 1;
        flags_.fkappa = 1;
        flags_.fsinth = 1;
        flags_.fzprim = 1;
}             
//************************ main function ***********************************
int main()
{

	int npar=0, ipar=0;
	vector <double> vecpar;
	vector <int> parnum;
	vector <string> parname;
	vector <double> stepsize;
	vector <double> lowlim;
	vector <double> ulim;
	vector <int> fixpars;
//********************************************************************************************//	
//****************************************** READING INPUT FILE ******************************//
//********************************************************************************************//
	//********** N.B.: paramter names should not have a space in them in input file***************//
	// First read the file that contains the fit params
	const char* infile="../input/smfit.dat";
	ifstream fin2(infile);
	if (!fin2.is_open()) { cout <<"Could not find input file: \""<< infile<<"\". STOP!!!"<<endl; return 1;}
	string line;
	cout<<"Reading parameter values from file: "<<infile<<endl;
	while(getline(fin2, line)) 
	{
		line.erase(line.begin(), find_if(line.begin(), line.end(), not1(ptr_fun<int, int>(isspace)))); // trim white spaces
		if(line[0] == '#') continue; // remove comment lines
		istringstream iss(line);
		if(iss.str()=="return" )break; // stop reading file
		if((iss.str()).empty())continue;

		vector<string> tokens;
		copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(tokens)); // split string stream to vector of strings
		cout<<iss.str()<<endl;
		//********************* READING PARAMTER DEFINITIONS *************************************//
		if (isdigit(tokens[0][0])){ // check if line starts with a digit: These lines are parameter definitions
			ipar++;
			if (tokens.size() != 4 && tokens.size() !=6 ){cout<<"Error reading input file: "<<infile<<"\t tokens.size() = "<< tokens.size()<<". STOPPING !!!"<<endl; return 1;} // Check if file is read properly
			//cout<<"\t size = "<<tokens.size()<<"\t"<<tokens[0]<<"\t"<<ipar<<endl;
			if(npar < atoi(tokens[0].c_str()) ){npar= atoi(tokens[0].c_str());} // save largest value of parameter index
			//*** Save to vectors ***//
			parnum.push_back(atoi(tokens[0].c_str())); // store parameter number
			parname.push_back(tokens[1]); // store parameter name
			vecpar.push_back(atof(tokens[2].c_str())); // store parameter value
			stepsize.push_back(atof(tokens[3].c_str())); // store starting stepsize
			if(tokens.size()==6){
				lowlim.push_back(atof(tokens[4].c_str()));
				ulim.push_back(atof(tokens[5].c_str()));
			}// if(tokens.size()>4)
			else{
				lowlim.push_back(0.0);
				ulim.push_back(0.0);
			}//else
		}//if (isdigit(tokens[0]))
		//******************* READING PARAMATER FIX *******************//
		if (boost:: iequals(tokens[0],"fix")){
			for (int k =1; k< tokens.size(); k++){
				fixpars.push_back(atoi(tokens[k].c_str()));
			}
			
		}//if (boost:: iequals(tokens[0],"fix"))
	}//while(getline(fin2, line)) 

//*****************************************************//	
//******************* Initialize flags ****************//	
//*****************************************************//	
	Initlog();
	
//****************************************************//
//***************** For Debugging ********************//	
//****************************************************//
	//cout<<"Finished reading file"<<endl;
	//cout<<"Total parameter lines read = "<<ipar<<endl;
	//cout<<"Total parameters = "<<npar<<endl;
	
	//double grad[npar]={0.},xval[npar]={0.};
	//int ntpar=30; // number of paramaters to minimize over
	//double fval=0.; // return of chi2 function in fcn
	//int iflag=1;
	
	//// fill in values of array xval
	//for(int i = 0; i<parnum.size() ; i++){xval[parnum[i] -1] = vecpar[i];}
	////for(int i = 0; i<npar ; i++){cout<<setprecision(5)<<"xval["<<i<<"] = "<<xval[i]<<endl;}
	
////	fcn_( &ntpar, grad, &fval, xval ,&iflag,chi2_);
////	cout<<"Finished call: chi2 =  "<<fval<<endl;
	
	//cpp_fcn( ntpar, grad, fval, xval, iflag);
	//cout<<"Finished call: chi2 =  "<<fval<<endl;
	
	
//***************************************************************************************//	
//************************ SET UP MINUIT ************************************************//
//***************************************************************************************//	
	cout<<"=============================="<<endl;
	cout<< "Start Minuit " << endl;
	TMinuit *gMinuit = new TMinuit(30);
	gMinuit -> SetFCN(cpp_fcn);
	double arglist[10];
	int ierflg = 0;
	arglist[0] = 1;
	//******************** SET PARAMETERS *************************//
	for (int k =0 ;k < parname.size();k++)
	{
		gMinuit->mnparm(parnum[k] -1,parname[k].c_str(),vecpar[k],stepsize[k],lowlim[k],ulim[k],ierflg);
	}//for (int k =0 ;k < parname.size();k++)
	//*********************FIX PARAMETERS *************************//
	for (int k =0 ;k < fixpars.size();k++)
	{
		gMinuit->FixParameter(fixpars[k] -1);
	}
	
	gMinuit->mnexcm("CALL FCN", arglist,1 ,ierflg);
	
	gMinuit->mnexcm("MINIMIZE", arglist,1 ,ierflg);
	gMinuit->mnexcm("MINIMIZE", arglist,1 ,ierflg);
	gMinuit->mnexcm("IMPROVE", arglist,1 ,ierflg);
	gMinuit->mnexcm("IMPROVE", arglist,1 ,ierflg);
	gMinuit->mnexcm("IMPROVE", arglist,1 ,ierflg);
	gMinuit->mnexcm("SEEK", arglist,1 ,ierflg);
	gMinuit->mnexcm("MINOS", arglist,1 ,ierflg);
	gMinuit->mnexcm("HESSE", arglist,1 ,ierflg);
	
	TGraph * g1 = (TGraph*) gMinuit->Contour(100,2-1,7-1);
	
	
	
	//gMinuit->Minimize();
	//gMinuit->Minimize();
	
	
	
	
	
	return 0;
}
