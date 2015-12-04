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
