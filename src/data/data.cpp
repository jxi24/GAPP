/**************************************************************
 * Class to read in data files in a given format taken from   *
 * the hepdata.cedar.ac.uk website. All files must be in the  *
 * format given by the plain text file format.                *
 *************************************************************/

#include "data.h"

Data::Data(std::string dataFile) {
    std::ifstream data;
    try{
        data.open(dataFile.c_str());

        if(data.is_open()) {
            ParseData(data);
        } else {
            throw -1;
        }
    } catch(int e) {
        std::cerr << "File cannont be opened, check to make sure " + dataFile + " exists." << std::endl;
        exit(-1);
    }
}

void Data::ParseData(std::ifstream& file) {
    std::string line;
    bool bData = false;
    int cnt = 0;
    while(!file.eof() && !bData) {
        std::getline(file,line);
        if(file.good()) {
            if(line.compare(0,2,"RE")==0) {
                std::size_t found = line.find("-");
                if(found!=std::string::npos) {
                    beam = line.substr(5,found-6);
                    process = line.substr(found+4,line.length()-found);
                    std::cout << beam << "-->" << process << std::endl;
                }
            } else if(line.compare(0,7,"SQRT(S)")==0) {
                std::size_t found = line.find("GeV");
                if(found!=std::string::npos) {
                    ecm = atof(line.substr(10,found-11).c_str());
                    std::cout << line.substr(10,found-11) << std::endl;
                    std::cout << ecm << std::endl;
                }
            } else if(line.compare(0,5,"xdesc")==0) {
                bData = true;
                std::stringstream ss(line);
                while(ss >> line) {
                    if(line == "y") cnt++;
                }
            }
        }
    }
    while(!file.eof()) {
        double x, y, yErr1Up, yErr2Up, yErr1Dn, yErr2Dn,dummy;
        file >> x >> dummy >> dummy >> y >> yErr1Up >> yErr1Dn >> yErr2Up >> yErr2Dn;
        for(int i = 0; i < cnt-1; i++) file >> dummy;
        xData.push_back(x);
        yData.push_back(y);
        yErrStatUp.push_back(yErr1Up);
        yErrStatDn.push_back(yErr1Dn);
        yErrSysUp.push_back(yErr2Up);
        yErrSysDn.push_back(yErr2Dn);
        std::cout << x << "\t" << y << "\t" << yErr1Up << "\t" << yErr1Dn << "\t" << yErr2Up << "\t" << yErr2Dn << std::endl;
    }
}
