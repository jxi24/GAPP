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
// Root includes
#include "TMinuit.h"
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
	cout<<"ntpar= "<<ntpar<<" iflag= "<<iflag<<endl;
	cout<<"fval= "<< fval <<endl;
	for (int i=0; i< ntpar;i++){cout<<"xval["<<i<<"]= "<<xval[i]<<endl;}
}

int main()
{

	int npar=0, ipar=0;
	vector <double> vecpar;
	vector <int> parnum;
	//****************************************** READING INPUT FILE ******************************//
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
		if(iss.str()=="return" || (iss.str()).empty())break; // stop reading file

		vector<string> tokens;
		copy(istream_iterator<string>(iss),istream_iterator<string>(),back_inserter(tokens)); // split string stream to vector of strings

		if (isdigit(tokens[0][0])){ // check if line starts with a digit
			cout<<iss.str();
			ipar++;
			cout<<"\t size = "<<tokens.size()<<"\t"<<tokens[0]<<"\t"<<ipar<<endl;
			if(npar < atoi(tokens[0].c_str()) ){npar= atoi(tokens[0].c_str());} // save largest value of parameter index
			vecpar.push_back(atof(tokens[2].c_str())); // store parameter value
			parnum.push_back(atoi(tokens[0].c_str())); // store parameter number
			}//if (isdigit(tokens[0]))
	}//while(getline(fin2, line)) 
	
	cout<<"Finished reading file"<<endl;
	cout<<"Total parameter lines read = "<<ipar<<endl;
	cout<<"Total parameters = "<<npar<<endl;
	
	double grad[npar]={0.},xval[npar]={0.};
	int ntpar=28; // number of paramaters to minimize over
	double fval=0.; // return of chi2 function in fcn
	int iflag=1;
	
	// fill in values of array xval
	for(int i = 0; i<parnum.size() ; i++){xval[parnum[i] -1] = vecpar[i];}
	for(int i = 0; i<npar ; i++){cout<<setprecision(5)<<"xval["<<i<<"] = "<<xval[i]<<endl;}
	
//	fcn_( &ntpar, grad, &fval, xval ,&iflag,chi2_);
//	cout<<"Finished call: chi2 =  "<<fval<<endl;
	
	cpp_fcn( ntpar, grad, fval, xval, iflag);
	cout<<"Finished call: chi2 =  "<<fval<<endl;
	
	return 0;
}
