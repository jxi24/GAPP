#ifndef DATA_H
#define DATA_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <vector>

class Data {
    public:
        Data(std::string dataFile);
        double GetXData(int i){return xData.at(i);};
        // GetYData returns the y prediction for a given x value
        double GetYData(double x);
        // Overloaded function to return the ith entry in the data
        double GetYData(int i){return yData.at(i);};
        // The following two getters returns the upper and lower errors respectively. The
        // error that is returned is given by: sqrt(Sys^2+Stat^2)
        double GetYUpError(double x);
        double GetYDnError(double x);
        // Overloaded functions for the ith error entry
        double GetYUpError(int i){return sqrt(yErrSysUp.at(i)*yErrSysUp.at(i)+yErrStatUp.at(i)+yErrStatUp.at(i));};
        double GetYDnError(int i){return sqrt(yErrSysDn.at(i)*yErrSysDn.at(i)+yErrStatDn.at(i)+yErrStatDn.at(i));};

    protected:

    private:
        void ParseData(std::ifstream& file);
        std::string beam;
        std::string process;
        double ecm;
        std::vector<double> xData;
        std::vector<double> yData;
        std::vector<double> yErrSysUp;
        std::vector<double> yErrSysDn;
        std::vector<double> yErrStatUp;
        std::vector<double> yErrStatDn;
};

#endif
