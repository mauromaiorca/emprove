/*
    (C) Mauro Maiorca, Mauro Maiorca, 2020, Birkbeck College

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <iostream>
#include <complex>
#include <fstream>
#include <list>
#include <iterator>     // std::back_inserter
#include <cstdlib>
#include <ctime>
#include <cstring>
#include <sstream>
#include <string>
#include <regex>

		
//ITK

#define PI 3.14159265358979323846
//#include "CsvStarReadWriteAnalyseLibs.h"
#include "mrcIO.h"

#define PBSTR "************************************************************"
#define PBWIDTH 60


  typedef float WorkingPixelType;
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;


// *************************************
// inputParametersType
// **************************************
void printProgress(double percentage) {
    int val = (int) (percentage * 100);
    int lpad = (int) (percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    printf("\r%3d%% [%.*s%*s]", val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

// *************************************
// inputParametersType
// **************************************
typedef struct inputParametersType {
    char * starFileIn;
    char * stackFileOut;
    float multiplicationFactor;
}inputParametersType;


/* ******************************************
 *  USAGE
 ***************************************** */
void usage(  char ** argv ){
    std::cerr<<"invertStack (c) Mauro Maiorca\n";
    std::cerr<<"\n";
    std::cerr<<"Usage: " << argv[0] << "\n";
    std::cerr<<"      stackInput.mrcs stackOutput.mrcs [factor=-1]\n";
    exit(1);
}




// **********************************
//
//    INT MAIN
//
// **********************************
int main( int argc, char **argv ){
    
    
    if (argc<3){
        usage(argv);
    }
    
    const char * inputStack=argv[1];
    const char * outputStack=argv[2];
    double multiplyFactor=-1;
    if (argc>3){
        multiplyFactor=atof(argv[3]);
    }
    std::cerr<<"inputStack="<<inputStack<<"\n";
    std::cerr<<"outputStack="<<outputStack<<"\n";
    std::cerr<<"multiplyFactor="<<multiplyFactor<<"\n";
        
    //create empty stack
    MRCHeader HeaderTmpFile;
    readHeaderMrc(inputStack, HeaderTmpFile);
    writeEmptyMrcImage(outputStack, HeaderTmpFile);
    
    
    //readStar(rlnImageNameVector, "_rlnImageName", parameters.starFileIn);
    //std::cerr<<"processing....\n";
    //float * mapTmp = new float [nxyz];
    //float * mapCtf = new float [nxyz];
    unsigned long int nx=HeaderTmpFile.nx;
    unsigned long int ny=HeaderTmpFile.ny;
    unsigned long int nz=HeaderTmpFile.nz;
    unsigned long int nxy=nx*ny;
    unsigned long int nxyz=nx*ny*nz;

    float * sliceTmp = new float [nxy];
    for (unsigned long int ii=0; ii<nz; ii++){
        printProgress((double)ii/nz);
        readMrcSlice(inputStack, sliceTmp, HeaderTmpFile, ii);
        for (unsigned long int ij=0; ij<nxy; ij++){
            sliceTmp[ij]*=multiplyFactor;
        }
        replaceMrcSlice(outputStack, sliceTmp, HeaderTmpFile, ii);
    }
    delete [] sliceTmp;

    
    
    



}
