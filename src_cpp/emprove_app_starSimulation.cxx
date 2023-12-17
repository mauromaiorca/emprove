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
#include "CsvStarReadWriteAnalyseLibs.h"
//#include "reconstructionLib.h"
//
//#include "fourierFunctions.h"
//#include "genericFilters.h"
//#include "msaEstimate.h"
//#include "msaLibs.h"
//#include "msaCsvReadWriteAnalyse.h"
//#include "plot2d.h"

#include "mrcIO.h"
//#include "msaLibs.h"
//#include "scoringFunctionsLib.h"
//#include "preprocessingLib.h"
//#include "projectionsLib.h"
//#include "fiboLibs.h"
//#include "otsuLib.h"
//#include "riaLibs.h"
//#include "analysisLib.h"
#include "euler_libs.h"
//#include "visualizationLibs.h"

  typedef float WorkingPixelType;
  typedef float OutputPixelType;
  const unsigned int Dimension = 3;

/*
    cartesianCoordsPhiTheta.push_back(cos(phi)*sin(theta));
    cartesianCoordsPhiTheta.push_back(sin(phi)*sin(theta));
    cartesianCoordsPhiTheta.push_back(sin(phi)*sin(theta));



    double XX = cartesianCoordsPhiTheta[0];
    double YY = cartesianCoordsPhiTheta[1];
    double ZZ = cartesianCoordsPhiTheta[2];
    phi=atan2(YY,XX)/PI_180;
    theta=atan2( pow(XX*XX+YY*YY,0.5), ZZ )/PI_180;
    refinedAnglesAndOrigin[0]=phi;
    refinedAnglesAndOrigin[1]=theta;
    XX = cartesianCoordsPsi[0];
    YY = cartesianCoordsPsi[1];
    psi=atan2(YY,XX)/PI_180;
*/




double wrapAngles2(double angle, double range=360.0){
   if (angle<0){
     double tmp=range+range-angle;
     long int tmp1=ceil(tmp/range);
     angle=angle+range*tmp1;
   }
   return fmod(angle,range);
}


double  induceAngleError(double AngleDegInput, double sigma, double normalization=1.0){
  double AngleIncrement = rnd_gaus(0, sigma)/normalization;
  double newAngle=AngleDegInput+AngleIncrement; //fmod(AngleDegInput+AngleIncrement,  360.0);
  return newAngle;
}

// **************************
// **************************
//    RANDOM ANGLES PICKING
// **************************
std::vector<double>  induceEulerError(double phiAngleDegInput, double thetaAngleDegInput, double sigmaPhiDeg, double sigmaThetaDeg, double normalization=1.0){

    double PI_over_180 = M_PI/180.0;

    double phi = phiAngleDegInput*PI_over_180;
    double theta = thetaAngleDegInput*PI_over_180;
    double sigmaPhi = sigmaPhiDeg*PI_over_180;
    double sigmaTheta = sigmaThetaDeg*PI_over_180;

    double XX = cos(phi)*sin(theta);
    double YY = sin(phi)*sin(theta);
    double ZZ = cos(theta);

    double sigma = pow((sigmaPhi+sigmaTheta)/2.0,2);
    
    bool goodGosition = false;
    double newPhi=phi;//fmod(phiAngleDegInput + PhiIncrement/SinTheta,360.0);
    double newTheta=theta;//fmod(thetaAngleDegInput + ThetaIncrement,180.0);
//    int hit =0;

    do {
      double XX_1 = XX + rnd_gaus(0, sigma );
      double YY_1 = YY + rnd_gaus(0, sigma );
      double ZZ_1 = ZZ + rnd_gaus(0, sigma );
      double distance = XX_1*XX_1 + YY_1*YY_1 +ZZ_1*ZZ_1;
      //double distance1 = pow(XX-XX_1,2.0) + pow(YY-YY_1, 2.0) + pow(ZZ-ZZ_1,2.0);
//      hit++;
//      std::cerr<<distance<<" ";
      if ( distance <= 1 && distance > 0.95 /*&& distance1 < sigma*/){
        newPhi=atan2(YY_1,XX_1)/PI_over_180;
        newTheta=atan2( pow(XX_1*XX_1+YY_1*YY_1,0.5), ZZ_1 )/PI_over_180;
        goodGosition = true;
      }
    } while ( !goodGosition );
//    std::cerr<<hit<<" ";
    
    std::vector<double> outputAngles;
    outputAngles.push_back(newPhi);
    outputAngles.push_back(newTheta);
    return outputAngles;
    
}


std::vector<double>  induceEulerErrorGlobal( ){

    double PI_over_180 = M_PI/180.0;

    double newPhi=0;
    double newTheta=0;
    double precision = 1000000;
    bool goodPosition = false;
    int hit = 0;

    do {
      double XX_1 = (fmod(rnd_double(), 2*precision)-precision)/precision;
      double YY_1 = (fmod(rnd_double(), 2*precision)-precision)/precision;
      double ZZ_1 = (fmod(rnd_double(), 2*precision)-precision)/precision;
      double distance = XX_1*XX_1 + YY_1*YY_1 +ZZ_1*ZZ_1;
      hit++;
      //std::cerr<<distance<<" ";
      if ( distance <= 1 && distance > 0.9 /*&& distance1 < sigma*/){
        newPhi=atan2(YY_1,XX_1)/PI_over_180;
        newTheta=atan2( pow(XX_1*XX_1+YY_1*YY_1,0.5), ZZ_1 )/PI_over_180;
        goodPosition = true;
      }
    } while ( !goodPosition );
    //std::cerr<<hit<<" ";

    
    std::vector<double> outputAngles;
    outputAngles.push_back(newPhi);
    outputAngles.push_back(newTheta);
    return outputAngles;
    
}






unsigned long int replaceAddValueStar2(std::string itemType, std::vector<std::string> replacingValuesString, const char * filenameIn, const char * filenameOut){
    std::string dataLabelStr(relionDataStartLabel(filenameIn));
    int dataLabelStrSize=dataLabelStr.size();
    

    std::cerr<<"replaceAddValueStar filenameIn="<<filenameIn <<"   filenameOut="<<filenameOut<<"\n"; 
    //check the strings
    std::string fileInStr(filenameIn);
    std::string fileOutStr(filenameOut);
    std::string tmpFile="";
    if (fileInStr.compare(filenameOut)==0){
      tmpFile=generateTmpFilename(".star");
      copyCvsFile(filenameIn, tmpFile.c_str() );
      fileInStr=tmpFile;
    }
    std::ifstream file(fileInStr);
    std::ofstream fileOut;
    fileOut.open(filenameOut);
    fileOut.close();
    fileOut.open (filenameOut, std::ofstream::out | std::ofstream::app);    
    std::string strLine;
    unsigned long int cc = 0;
    bool got_data_images = false;
    bool got_loop = false;
    std::vector<std::string> fields;
    std::vector<int> fieldsIdx;
    int targetID = -1;
    unsigned long int headerLines = 0;

    bool newLabel = false;
    if ( getStarHeaderItemIdx(itemType, fileInStr.c_str()) < 0 ){
      newLabel = true;
    }
    

    while (std::getline(file, strLine) && ++cc < __MAX_STAR_HEADER_SIZE__){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::regex rNoSpaces("\\s+");
	std::string strNoSpaces = std::regex_replace(strLine, rNoSpaces, ",");    
       if (strNoSpaces.length()<1){
         //do nothing
         fileOut<<strLine<<"\n";
       }else if (str0.substr(0,1).compare("#")==0){
         //it is a comment, do nothing
         fileOut<<strLine<<"\n";
       }else if (!got_data_images){
          fileOut<<strLine<<"\n";
          if( strstr(strNoSpaces.substr(0,dataLabelStrSize).c_str(),dataLabelStr.c_str()) ) {
            got_data_images=true;
            //std::cerr<<"star file: data_images\n";
          }
       }else if (!got_loop){
          fileOut<<strLine<<"\n";
          if( strstr(str0.substr(0,5).c_str(),"loop_") ){
            got_loop=true;
            //std::cerr<<"got_loop\n";
          }
       }else if( strstr(str0.substr(0,1).c_str(),"_") ){ //was _rln
        fileOut<<strLine<<"\n";
        headerLines=cc+1;
        std::size_t found = str0.find(",");
        std::string headerItem=str0.substr(0,found);
        std::size_t start1=str0.find("#");
        //std::size_t start1 = headerItem.find(",");
        std::string columnItem=str0.substr(start1+1, str0.length()-1);
        std::size_t start2=columnItem.find(",");
        columnItem=str0.substr(start1+1, start2);
        int columnItemValue=atoi(columnItem.c_str())-1;
        fields.push_back(headerItem);
        fieldsIdx.push_back(columnItemValue);
        if (headerItem.compare(itemType)==0){
          targetID=columnItemValue;
        }
       }
    }

    if ( newLabel ){
      int nextIdx = StarMaxFieldIdx(fileInStr.c_str())+1;
      fileOut<<itemType<< " #"<< std::to_string(nextIdx)  <<" \n";
    }   
    unsigned long int ccCounter=0;
    file.clear();
    file.seekg(0,std::ios::beg);
    unsigned long int counter = 0;
    while (std::getline(file, strLine) ){
     if (++ccCounter >= headerLines){
        std::regex r("\\s+");
	std::string str0 = std::regex_replace(strLine, r, ",");
        std::string result=findColumnItem(str0, targetID);

        std::regex r1("[,]+");
	std::string str1 = std::regex_replace(str0,r1," ");

        std::regex r2("\\s+");
	std::string str2 = std::regex_replace(str1, r2, "");        

        if ( newLabel ){
           fileOut<<str1<<" ";
                if (counter<replacingValuesString.size()){
                  fileOut<<replacingValuesString[counter];
                }else{
                  if (str2.size()>0){
                   fileOut<<"0";
                  }
                }
            fileOut<<" ";
            counter++;
        }else{
                int beforePos=findColumnItemPosition(str0, targetID);
	        if (beforePos>=0){
          	  std::string beforeStr = str1.substr(0,beforePos);
	          fileOut<<beforeStr;
	        }
	        if (beforePos>0){
	          fileOut<<" ";
	        }
                //OK fileOut<<replacingString;
                if (counter<replacingValuesString.size()){
                  fileOut<<replacingValuesString[counter];
                }else{
                  if (str2.size()>0){
                   fileOut<<"0";
                  }
                }
	        fileOut<<" ";
	        counter++;
	        std::string afterStr = str1.substr( beforePos+result.length() );
	        fileOut<<afterStr;
        }
        fileOut<<"\n";
      }
    }
    fileOut.close();

    if (tmpFile.length()>0){
      removeCvsFile(tmpFile.c_str());
    }

    return counter;
}



unsigned long int replaceAddValueStar2(std::string itemType, std::vector<double> replacingValuesDouble, const char * filenameIn, const char * filenameOut){
    std::vector<std::string> replacingValuesString;
    for (unsigned long int ii=0; ii<replacingValuesDouble.size();ii++){
        replacingValuesString.push_back(std::to_string(replacingValuesDouble[ii]));
    }
    return replaceAddValueStar2(itemType,replacingValuesString,filenameIn,filenameOut);
}

unsigned long int replaceAddValueStar2(std::string itemType, std::vector<int> replacingValuesInt, const char * filenameIn, const char * filenameOut){
    std::vector<std::string> replacingValuesString;
    for (unsigned long int ii=0; ii<replacingValuesInt.size();ii++){
        replacingValuesString.push_back(std::to_string(replacingValuesInt[ii]));
    }
    return replaceAddValueStar2(itemType,replacingValuesString,filenameIn,filenameOut);
}



std::vector<double>  induceEulerError32(double phiAngleDegInput, double thetaAngleDegInput, double sigmaPhiDeg, double sigmaThetaDeg, double normalization=1.0){

    double PI_over_180 = M_PI/180.0;

    double phi = phiAngleDegInput*PI_over_180;
    double theta = thetaAngleDegInput*PI_over_180;
    double sigmaPhi = sigmaPhiDeg*PI_over_180;
    double sigmaTheta = sigmaThetaDeg*PI_over_180;

    double XX = cos(phi)*sin(theta);
    double YY = sin(phi)*sin(theta);
    double ZZ = cos(theta);


//    double XX_sigma = cos(sigmaPhi)*sin(sigmaTheta);
//    double YY_sigma = sin(sigmaPhi)*sin(sigmaTheta);
//    double ZZ_sigma = sin(sigmaPhi)*sin(sigmaTheta);

    double sigma = pow((sigmaPhi+sigmaTheta)/2.0,2);
    //std::cerr<<"\n\n"<<sigma<<" "<<sigmaPhi<<" "<<sigmaTheta<<"    ";
    //std::cerr<<"\n"<<XX<<" "<<YY<<" "<<ZZ<<"    ";
    
    bool goodGosition = false;
    double newPhi=phi;//fmod(phiAngleDegInput + PhiIncrement/SinTheta,360.0);
    double newTheta=theta;//fmod(thetaAngleDegInput + ThetaIncrement,180.0);
    int hit =0;

    do {
      double XX_1 = XX + rnd_gaus(0, sigma );
      double YY_1 = YY + rnd_gaus(0, sigma );
      double ZZ_1 = ZZ + rnd_gaus(0, sigma );
      double distance = XX_1*XX_1 + YY_1*YY_1 +ZZ_1*ZZ_1;
      //double distance1 = pow(XX-XX_1,2.0) + pow(YY-YY_1, 2.0) + pow(ZZ-ZZ_1,2.0);
      hit++;
      //std::cerr<<distance<<" ";
      if ( distance <= 1 && distance > 0.9 /*&& distance1 < sigma*/){
        newPhi=atan2(YY_1,XX_1)/PI_over_180;
        newTheta=atan2( pow(XX_1*XX_1+YY_1*YY_1,0.5), ZZ_1 )/PI_over_180;
        goodGosition = true;
      }
    } while ( !goodGosition );
    //std::cerr<<hit<<" ";
    //sigmaTheta


/*
    double PhiIncrement = rnd_gaus(0, sigmaPhi)/normalization;
    double ThetaIncrement = rnd_gaus(0, sigmaTheta)/normalization;
    double phiAngleRad=phiAngleDegInput*PI_on_180;
    double thetaAngleRad=thetaAngleDegInput*PI_on_180;

    double SinTheta=sin(thetaAngleRad);
    if ( SinTheta < 0.00001 && thetaAngleRad >= 0.0 ){
        SinTheta=0.00001;
    }else if ( SinTheta > -0.00001 && thetaAngleRad <= 0.0 ){
        SinTheta=-0.00001;
    }
*/
    
    std::vector<double> outputAngles;
    outputAngles.push_back(newPhi);
    outputAngles.push_back(newTheta);
    return outputAngles;
    
}







// *************************************
// inputParametersType
// **************************************
typedef struct inputParametersType {
    long int numberOfUniformViews;
    char * uniformViewStackFilename;
    char * templateStarFileUniformViews;
    char * inputFileToShake;
    double sigmaShake;
    char * parametersToShake;
    double percentageToShake;
    char * shakingString;

    char * inputFileToPreferredOrientation;
    double eulerPreferredOrientationPhi;
    double eulerPreferredOrientationTheta;
    double sigmaPreferredOrientation;
    double percentagetoPreferredOrientation;
    char * starFileOut;
    double angpix;
    
    bool addCTF;
    char * inputFileToCtf;
    double voltage;
    double defocusU;
    double defocusV;
    double defocusAngle;
    double sphAberration;
    double AmplitudeContrast;
    double voltageVariance;
    double defocusUVariance;
    double defocusVVariance;
    double defocusAngleVariance;
    double sphAberrationVariance;
    double AmplitudeContrastVariance;


    char * inputFileToUniform;
    double existingToUniformNumberViews;
}inputParametersType;


/* ******************************************
 *  USAGE
 ***************************************** */
void usage(  char ** argv ){
    std::cerr<<"starSimulation (c) Mauro Maiorca\n";
    std::cerr<<"\n";
    std::cerr<<"Usage: " << argv[0] << "\n";
    std::cerr<<"      --NewUniform numberOfUniformViews [uniformViewStackFilename.mrcs templateStarFileUniformViews.star]\n";
    std::cerr<<"      --existingToUniform inputFileToUniform.star existingToUniformNumberViews\n";
    std::cerr<<"      --shakeParameters inputFileToShake.star shakingString=parametersToShake:sigmaToShake:PercentageToShake\n";
    std::cerr<<"                                parametersToShake=origin,psi,global,defocusU,defocusV,defocusAngle  sigmaToShake=Angle_in_degrees  percentageToShake=[0..1]\n";
//    std::cerr<<"      --shakeParameters inputFileToShake.star sigmaShake parametersToShake=euler,origin,psi percentageToShake=[0..1]\n";
//    std::cerr<<"                      (NOTE: it does just psi at the moment, sorry)\n";
    std::cerr<<"      --preferredOrientation inputFileToPreferredOrientation.star eulerPreferredOrientationPhi eulerPreferredOrientationTheta sigmaPreferredOrientation percentagetoPreferredOrientation=[0..1]\n";
    std::cerr<<"      --ctf inputFileToCtf.star voltage defocusU defocusV defocusAngle sphAberration AmplitudeContrast\n";
    std::cerr<<"      --ctfVarianceError  voltageVariance defocusUVariance defocusVVariance defocusAngleVariance sphAberrationVariance AmplitudeContrastVariance\n";
    std::cerr<<"      --angpix angpix\n";
    std::cerr<<"      --o starFileOut.star\n";
    std::cerr<<"      --h (help)\n";
    std::cerr<<"\n";
    exit(1);
}

// *************************************
//
// retrieveInputParameters
// *************************************
void retrieveInputParameters(inputParametersType * parameters, int argc, char** argv){
    if ( argc < 2)
	usage(argv);

    parameters->numberOfUniformViews= -1;
    parameters->uniformViewStackFilename= NULL;
    parameters->templateStarFileUniformViews= NULL;
    parameters->inputFileToShake= NULL;
    parameters->sigmaShake= -1;
    parameters->parametersToShake= NULL;
    parameters->shakingString=NULL;
    parameters->percentageToShake= 1;
    parameters->inputFileToPreferredOrientation= NULL;
    parameters->eulerPreferredOrientationPhi= 0;
    parameters->eulerPreferredOrientationTheta= 0;
    parameters->sigmaPreferredOrientation= -1;
    parameters->percentagetoPreferredOrientation= -1;
    parameters->starFileOut= NULL;
    parameters->angpix= 1.0;
    

    parameters->addCTF=false;
    parameters->inputFileToCtf=NULL;
    parameters->voltage= 300;
    parameters->defocusU= 14558.570312;
    parameters->defocusV= 14860.139648;
    parameters->defocusAngle= 18.559999;
    parameters->sphAberration= 2.7;
    parameters->AmplitudeContrast= 0.1;


    parameters->voltageVariance= 0;
    parameters->defocusUVariance= 0;
    parameters->defocusVVariance= 0;
    parameters->defocusAngleVariance= 0;
    parameters->sphAberrationVariance= 0;
    parameters->AmplitudeContrastVariance= 0;


    parameters->inputFileToUniform= NULL;
    parameters->existingToUniformNumberViews= 1.0;
    
    for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if (!optionStr.compare("--h") ){
            usage(argv);
        }
    }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--NewUniform") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->numberOfUniformViews=atoi(argv[jj]);
            if (idx==1) parameters->uniformViewStackFilename=argv[jj];
            if (idx==2) parameters->templateStarFileUniformViews=argv[jj];  
          }
          ii+=idx;
        }
  }

    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--existingToUniform") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->inputFileToUniform=argv[jj];
              if (idx==1) parameters->existingToUniformNumberViews=atoi(argv[jj]);
            }
            ii+=idx;
          }
    }

    for (unsigned int ii=1;ii<argc;ii++){
          std::string optionStr(argv[ii]);
          if ( !optionStr.compare("--angpix") ){
            int idx=0;
            for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
              std::string subparamStr(argv[jj]);
              if (!subparamStr.substr(0,2).compare("--")) break;
              if (idx==0) parameters->angpix=atof(argv[jj]);
            }
            ii+=idx;
          }
    }
  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--shakeParameters") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->inputFileToShake=argv[jj];
            if (idx==1) parameters->shakingString=argv[jj];  
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--preferredOrientation") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->inputFileToPreferredOrientation=argv[jj];
            if (idx==1) parameters->eulerPreferredOrientationPhi=atof(argv[jj]);  
            if (idx==2) parameters->eulerPreferredOrientationTheta=atof(argv[jj]);  
            if (idx==3) parameters->sigmaPreferredOrientation=atof(argv[jj]);
            if (idx==4) parameters->percentagetoPreferredOrientation=atof(argv[jj]);
          }
          ii+=idx;
        }
  }


  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--ctf") ){
          parameters->addCTF=true;
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->inputFileToCtf=argv[jj];
            if (idx==1) parameters->voltage=atof(argv[jj]);
            if (idx==2) parameters->defocusU=atof(argv[jj]);  
            if (idx==3) parameters->defocusV=atof(argv[jj]);  
            if (idx==4) parameters->defocusAngle=atof(argv[jj]);
            if (idx==5) parameters->sphAberration=atof(argv[jj]);
            if (idx==6) parameters->AmplitudeContrast=atof(argv[jj]);
          }
          ii+=idx;
        }
  }

  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--ctfVarianceError") ){
          parameters->addCTF=true;
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->voltageVariance=atof(argv[jj]);
            if (idx==1) parameters->defocusUVariance=atof(argv[jj]);  
            if (idx==2) parameters->defocusVVariance=atof(argv[jj]);  
            if (idx==3) parameters->defocusAngleVariance=atof(argv[jj]);
            if (idx==4) parameters->sphAberrationVariance=atof(argv[jj]);
            if (idx==5) parameters->AmplitudeContrastVariance=atof(argv[jj]);
          }
          ii+=idx;
        }
  }





  for (unsigned int ii=1;ii<argc;ii++){
        std::string optionStr(argv[ii]);
        if ( !optionStr.compare("--o") ){
          int idx=0;
          for (unsigned int jj=ii+1; jj<argc ; jj++, idx++){
            std::string subparamStr(argv[jj]);
            if (!subparamStr.substr(0,2).compare("--")) break;
            if (idx==0) parameters->starFileOut=argv[jj];
          }
          ii+=idx;
        }
  }
}




// **********************************
//
//    INT MAIN
//
// **********************************
int main( int argc, char **argv ){

 inputParametersType parameters;
 retrieveInputParameters(&parameters, argc, argv);

 // *********************************
 // *****************************
 //   Uniform views for a new file
 //
 // *****************************
 if (parameters.numberOfUniformViews > 0 && parameters.uniformViewStackFilename ){
  std::cerr<<"to reimplement\n";
  /*
   if (!parameters.starFileOut){
     std::cerr<<"ERROR: define output star file\n";
     exit(1);
   }
   //std::cerr<<"      --NewUniform numberOfUniformViews [uniformViewStackFilename.mrcs templateStarFileUniformViews.star]\n";
   std::vector<double> phiV;
   std::vector<double> thetaV;
   std::vector<double> subsetV;
   long int views = generateEquallyDistributedViews (phiV, thetaV, parameters.numberOfUniformViews);
   
        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);  
        
        std::ifstream fileParticles(parameters.templateStarFileUniformViews);
        std::string strLine;
        if (parameters.templateStarFileUniformViews){
                long int startMicrograph=getStarStart(parameters.templateStarFileUniformViews);
                for (int counter=0;counter<startMicrograph;counter++){
                     std::getline(fileParticles, strLine);           
                     fileOutput << strLine <<"\n";
                   }
        }else{
	        fileOutput << "data_\n";
	        fileOutput << "loop_\n";
	        fileOutput << "_rlnImageName #1\n";
	        fileOutput << "_rlnAngleRot #2\n";
	        fileOutput << "_rlnAngleTilt #3\n";
	        fileOutput << "_rlnAnglePsi #4\n";
	        fileOutput << "_rlnOriginX #5\n";
	        fileOutput << "_rlnOriginY #6\n";
	        fileOutput << "_rlnCoordinateX #7\n";
	        fileOutput << "_rlnCoordinateY #8\n";
            fileOutput << "_rlnRandomSubset #9\n";

            fileOutput << "_rlnDetectorPixelSize #10\n";
            //fileOutput << "_rlnVoltage #11\n";
            //fileOutput << "_rlnSphericalAberration #12\n"; //0.001000

        }
        std::string templateString;
        int imageNameIdx;
        int phiIdx;
        int thetaIdx; 
        int psiIdx;
        int originXIdx;
        int originYIdx;
        int coordXIdx;
        int coordYIdx;
        int subsetIdx=-1;
     int pixelSizeIdx;
     //int voltageIdx;
     //int sphericalAberrationIdx;

        bool weHaveSubsets=false;        
        if (parameters.templateStarFileUniformViews){
             std::getline(fileParticles, strLine);
             templateString=strLine;
             //remove leading string 
             while(!templateString.empty() && std::isspace(*templateString.begin()))
               templateString.erase(templateString.begin());
             phiIdx=getStarHeaderItemIdx("_rlnAngleRot", parameters.templateStarFileUniformViews);
             thetaIdx=getStarHeaderItemIdx("_rlnAngleTilt", parameters.templateStarFileUniformViews);
             psiIdx=getStarHeaderItemIdx("_rlnAnglePsi", parameters.templateStarFileUniformViews);
             originXIdx=getStarHeaderItemIdx("_rlnOriginX", parameters.templateStarFileUniformViews);
             originYIdx=getStarHeaderItemIdx("_rlnOriginY", parameters.templateStarFileUniformViews);
             coordXIdx=getStarHeaderItemIdx("_rlnCoordinateX", parameters.templateStarFileUniformViews);
             coordYIdx=getStarHeaderItemIdx("_rlnCoordinateY", parameters.templateStarFileUniformViews);
             subsetIdx=getStarHeaderItemIdx("_rlnRandomSubset", parameters.templateStarFileUniformViews);
             imageNameIdx=getStarHeaderItemIdx("_rlnImageName", parameters.templateStarFileUniformViews);
             pixelSizeIdx=getStarHeaderItemIdx("_rlnDetectorPixelSize", parameters.templateStarFileUniformViews);
             //voltageIdx=getStarHeaderItemIdx("_rlnVoltage", parameters.templateStarFileUniformViews);
             //sphericalAberrationIdx=getStarHeaderItemIdx("_rlnSphericalAberration", parameters.templateStarFileUniformViews);
            
            if(subsetIdx<0 ){
               subsetIdx=StarMaxFieldIdx(parameters.templateStarFileUniformViews);
               fileOutput << "_rlnRandomSubset #" <<std::to_string(subsetIdx+1) <<"\n";
             }else {
                weHaveSubsets=true;
                readStar(subsetV, "_rlnRandomSubset", parameters.templateStarFileUniformViews);
             }
        }else{
             templateString=std::string("0001@")+std::string(parameters.uniformViewStackFilename)+std::string(" 0.0000  0.0000  0.0000  0.0000  0.0000   0.0000  0.0000  0 0.0000");
             imageNameIdx=0;
             phiIdx=1;
             thetaIdx=2;
             psiIdx=3;
             originXIdx=4;
             originYIdx=5;
             coordXIdx=6;
             coordYIdx=7;
	         subsetIdx=8;
             pixelSizeIdx=9;
             //voltageIdx=10;
             //sphericalAberrationIdx=11;
        }


	for (unsigned long int ii=0, counterOutput=1; ii<phiV.size(); ii++, counterOutput++ ){
	     std::string tmpTemplateString=templateString;
             std::string newMicrographNameValue=zeros_lead_to_string(counterOutput, 6)+std::string("@")+parameters.uniformViewStackFilename;
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, imageNameIdx, newMicrographNameValue);
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, phiIdx, std::to_string(phiV[ii]));
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, thetaIdx, std::to_string(thetaV[ii]));
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, psiIdx, "0.00000");
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, originXIdx, "0.00000");
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, originYIdx, "0.00000");
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, coordXIdx, "0.00000");
             tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, coordYIdx, "0.00000");
        
        tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, pixelSizeIdx, std::to_string(parameters.angpix) );
        //tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, voltageIdx, "300");
        //tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, sphericalAberrationIdx, "0.001000");

	     if(weHaveSubsets){
		 //int valueToReplace=subsetV[ii];
	         int valueToReplace=ii%2+1; 
		 //std::cerr<<valueToReplace<< " ";
                 tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, subsetIdx, std::to_string(valueToReplace));
             }else{
                 int valueToReplace=ii%2+1; 
                 tmpTemplateString=replaceValueStrlineStarFile(tmpTemplateString, subsetIdx, std::to_string(valueToReplace));
	     }


             while(!tmpTemplateString.empty() && std::isspace(*tmpTemplateString.begin()))
               tmpTemplateString.erase(tmpTemplateString.begin());             
             fileOutput << tmpTemplateString << "\n";
	}

        fileOutput.close();
        */
 }
 // *********************************
 // *****************************
 //   existingToUniform
 //
 // *****************************
 else if (parameters.inputFileToUniform && parameters.existingToUniformNumberViews > 0 && parameters.starFileOut){
  std::cerr<<"to reimplement\n";
/*
    std::vector<double> phiListIn;
    std::vector<double> thetaListIn;
    std::vector<double> phiListUniform;
    std::vector<double> thetaListUniform;
    readStar(phiListIn, "_rlnAngleRot", parameters.inputFileToUniform);
    readStar(thetaListIn, "_rlnAngleTilt", parameters.inputFileToUniform);
    //int phiIdx=getStarHeaderItemIdx("_rlnAngleRot", parameters.inputFileToUniform);
    //int thetaIdx=getStarHeaderItemIdx("_rlnAngleTilt", parameters.inputFileToUniform);
    assignClosestUniformEulerAngle ( phiListIn, thetaListIn, phiListUniform, thetaListUniform, parameters.existingToUniformNumberViews);
    //--existingToUniform inputFileToUniform.star existingToUniformNumberViews
    replaceAddValueStar("_rlnAngleRot", phiListUniform, parameters.inputFileToUniform, "____tmp_file___.star");
    replaceAddValueStar("_rlnAngleTilt", thetaListUniform, "____tmp_file___.star", parameters.starFileOut);
    removeCvsFile("____tmp_file___.star");
    */

 }
 // *********************************
 // *****************************
 //   shakeParameters 
 //
 // *****************************
 else if (parameters.inputFileToShake && parameters.shakingString){
      // --shakeParameters inputFileToShake.star sigmaShake parametersToShake=euler,origin,psi percentageToShake=[0..1]
      //std::cerr<<"      --shakeParameters inputFileToShake.star shakingString=parametersToShake:sigmaToShake:PercentageToShake\n";

        std::cerr<<"parameters.inputFileToShake="<<parameters.inputFileToShake<<"\n";
        if (!parameters.starFileOut){
          std::cerr<<"ERROR: define output star file\n";
          exit(1);
        }

        std::vector<std::vector<std::string> > optionShaking = getCsvOptionParameters(parameters.shakingString, 3);

        for (int kk=0; kk<optionShaking.size();kk++ ){
          std::cerr <<"debug optionShaking="<<optionShaking[kk][0] <<"\n";
          if (std::string("psi") == optionShaking[kk][0] ){
            std::cerr<<"shaking psi, sigma(degrees)="<< optionShaking[kk][1];
            std::vector <double> psi; 
            readStar(psi, "_rlnAnglePsi", parameters.inputFileToShake);
            for ( long int ii=0; ii<psi.size(); ii++ ){
                  double randomValue = static_cast<double>(rand()) / RAND_MAX;
                  if (randomValue<std::stod(optionShaking[kk][2])){
                    psi[ii] = induceAngleError(psi[ii], std::stod(optionShaking[kk][1]) );
                  }

            }
            replaceAddValueStar("_rlnAnglePsi", psi, parameters.inputFileToShake, parameters.starFileOut);
          }
          if (std::string("global") == optionShaking[kk][0] ){
            std::cerr<<"global shaking sigma="<<optionShaking[kk][1]<<"\n";
            double percentage = 1.0;
            if ( atof(optionShaking[kk][2].c_str()) < 1 ){
              percentage=atof(optionShaking[kk][2].c_str());
            }
            std::cerr<<"percentage shaking="<<optionShaking[kk][2]<<"\n";
            

            //std::vector <double> psi; 
            std::vector <double> rlnTheta;
            std::vector <double> rlnPhi;
            readStar(rlnPhi, "_rlnAngleRot", parameters.inputFileToShake);
            readStar(rlnTheta, "_rlnAngleTilt", parameters.inputFileToShake);	
            //readStar(psi, "_rlnAnglePsi", parameters.inputFileToShake);
            for ( long int ii=0; ii<rlnPhi.size(); ii++ ){
              double chance = fmod(rnd_double(),1000.0)/1000.0;
              if ( chance < percentage + 0.0001 ){
                //psi[ii]=fmod(rnd_double(),360000.0)/1000.0;
                //std::vector<double>  phiTheta = induceEulerError(0, 0 , 90, 90);
                std::vector<double>  phiTheta = induceEulerError32(rlnPhi[ii], rlnTheta[ii], atof(optionShaking[kk][1].c_str()), atof(optionShaking[kk][1].c_str()) );
                //std::vector<double>  phiTheta = induceEulerError32(rlnPhi[ii], rlnTheta[ii], 1, 1 );
                //std::vector<double>  induceEulerError32(double phiAngleDegInput, double thetaAngleDegInput, double sigmaPhiDeg, double sigmaThetaDeg, double normalization=1.0){

                
                rlnPhi[ii]=phiTheta[0];
                rlnTheta[ii]=phiTheta[1];
              }
            }
            replaceAddValueStar("_rlnAngleRot", rlnPhi, parameters.inputFileToShake, parameters.starFileOut);
            replaceAddValueStar("_rlnAngleTilt", rlnTheta, parameters.starFileOut, parameters.starFileOut);
            //replaceAddValueStar("_rlnAnglePsi", psi, parameters.starFileOut, parameters.starFileOut);
          }
          if (std::string("origin") == optionShaking[kk][0] ){
            std::cerr<<"shaking origin, sigma(degrees)="<< optionShaking[kk][1]<<"\n";

            std::vector <double> originX; 
            std::vector <double> originY; 
            readStar(originX, "_rlnOriginXAngst", parameters.inputFileToShake );
            readStar(originY, "_rlnOriginYAngst", parameters.inputFileToShake );
            for ( long int bbb=0; bbb<originX.size(); bbb++ ){
                  double randomValue = static_cast<double>(rand()) / RAND_MAX;
                  if (randomValue<std::stod(optionShaking[kk][2])){

//                  std::cerr << "originX=" << originX[bbb] << "   ";
//                  std::cerr << "shaking=" << std::stod(optionShaking[kk][1]) << "   ";
                    originX[bbb] = originX[bbb] + rnd_gaus(0, std::stod(optionShaking[kk][1]) );
                    originY[bbb] = originY[bbb] + rnd_gaus(0, std::stod(optionShaking[kk][1]) );
//                  std::cerr<<"   final="<<originX[bbb] <<"\n";
                  }
            }

//            std::cerr<<"parameters.starFileOut="<<parameters.starFileOut<<"\n";
            replaceAddValueStar("_rlnOriginXAngst", originX, parameters.inputFileToShake, parameters.starFileOut);
            replaceAddValueStar("_rlnOriginYAngst", originY, parameters.starFileOut, parameters.starFileOut);
          }     


          if (std::string("defocusU") == optionShaking[kk][0] ){
            std::vector<double> defocusU;
            readStar(defocusU, "_rlnDefocusU", parameters.inputFileToShake );
            for ( long int bbb=0; bbb<defocusU.size(); bbb++ ){
                  double randomValue = static_cast<double>(rand()) / RAND_MAX;
                  if (randomValue<std::stod(optionShaking[kk][2])){
                    defocusU[bbb] = defocusU[bbb] + rnd_gaus(0, std::stod(optionShaking[kk][1]) );
                  }
            }
            replaceAddValueStar("_rlnDefocusU", defocusU, parameters.starFileOut, parameters.starFileOut);
          }

          if (std::string("defocusV") == optionShaking[kk][0] ){
            std::vector<double> defocusV;
            readStar(defocusV, "_rlnDefocusV", parameters.inputFileToShake );
            for ( long int bbb=0; bbb<defocusV.size(); bbb++ ){
                  double randomValue = static_cast<double>(rand()) / RAND_MAX;
                  if (randomValue<std::stod(optionShaking[kk][2])){
                    defocusV[bbb] = defocusV[bbb] + rnd_gaus(0, std::stod(optionShaking[kk][1]) );
                  }
            }
            replaceAddValueStar("_rlnDefocusV", defocusV, parameters.starFileOut, parameters.starFileOut);
          }

          if (std::string("defocusAngle") == optionShaking[kk][0] ){
            std::vector<double> defocusAngle;
            readStar(defocusAngle, "_rlnDefocusAngle", parameters.inputFileToShake );
            for ( long int bbb=0; bbb<defocusAngle.size(); bbb++ ){
                  double randomValue = static_cast<double>(rand()) / RAND_MAX;
                  if (randomValue<std::stod(optionShaking[kk][2])){
                    defocusAngle[bbb] = defocusAngle[bbb] + rnd_gaus(0, std::stod(optionShaking[kk][1]) );
                  }
            }
            replaceAddValueStar("_rlnDefocusAngle", defocusAngle, parameters.starFileOut, parameters.starFileOut);
          }

        }

 }

 // *********************************
 // *****************************
 //   add ctf parameter 
 //
 // *****************************
 else if (parameters.inputFileToCtf && parameters.addCTF ){
  long int numItemsStar=getNumItemsStar(parameters.inputFileToCtf);
//  std::cerr<<"      --ctf inputFileToCtf.star voltage defocusU defocusV defocusAngle sphAberration AmplitudeContrast\n";
  
  
  
  std::vector<double> voltage;
  std::vector<double> defocusU;
  std::vector<double> defocusV;
  std::vector<double> defocusAngle;
  std::vector<double> sphAberration;
  std::vector<double> AmplitudeContrast;



  for (unsigned long int ii =0; ii<numItemsStar; ii++){
    voltage.push_back( rnd_gaus(parameters.voltage, parameters.voltageVariance) );
    defocusU.push_back( rnd_gaus(parameters.defocusU, parameters.defocusUVariance) );
    defocusV.push_back( rnd_gaus(parameters.defocusV, parameters.defocusVVariance) );
    defocusAngle.push_back( rnd_gaus(parameters.defocusAngle, parameters.defocusAngleVariance) );
    sphAberration.push_back( rnd_gaus(parameters.sphAberration, parameters.sphAberrationVariance) );
    AmplitudeContrast.push_back( rnd_gaus(parameters.AmplitudeContrast, parameters.AmplitudeContrastVariance) );
  }

  replaceAddValueStar( "_rlnVoltage", voltage, parameters.inputFileToCtf, parameters.inputFileToCtf);
  replaceAddValueStar("_rlnDefocusU", defocusU, parameters.inputFileToCtf, parameters.inputFileToCtf);	
  replaceAddValueStar("_rlnDefocusV", defocusV, parameters.inputFileToCtf, parameters.inputFileToCtf);
  replaceAddValueStar("_rlnDefocusAngle", defocusAngle, parameters.inputFileToCtf, parameters.inputFileToCtf);	
  replaceAddValueStar("_rlnSphericalAberration", sphAberration, parameters.inputFileToCtf, parameters.inputFileToCtf);
  replaceAddValueStar("_rlnAmplitudeContrast", AmplitudeContrast, parameters.inputFileToCtf, parameters.inputFileToCtf);	

 }
 // *********************************
 // *****************************
 //   add preferred orientation 
 //
 // *****************************
 else if (parameters.inputFileToPreferredOrientation && parameters.sigmaPreferredOrientation >=0  && parameters.percentagetoPreferredOrientation >0 ){
 //--preferredOrientation inputFileToPreferredOrientation.star eulerPreferredOrientationPhi eulerPreferredOrientationTheta sigmaPreferredOrientation percentagetoPreferredOrientation=[0..1]
   if (!parameters.starFileOut){
     std::cerr<<"ERROR: define output star file\n";
     exit(1);
   }
   std::vector<double> rlnPhi;
   std::vector<double> rlnTheta;
   readStar(rlnPhi, "_rlnAngleRot", parameters.inputFileToPreferredOrientation);
   readStar(rlnTheta, "_rlnAngleTilt", parameters.inputFileToPreferredOrientation);	
   int phiIdx=getStarHeaderItemIdx("_rlnAngleRot", parameters.inputFileToPreferredOrientation);
   int thetaIdx=getStarHeaderItemIdx("_rlnAngleTilt", parameters.inputFileToPreferredOrientation);
   
   std::vector<long int> toProcessIdx;
  for (long int ii=0; ii<rlnPhi.size(); ii++){
      toProcessIdx.push_back(ii);
  }
  std::srand(std::time(0));
  std::random_shuffle ( toProcessIdx.begin(), toProcessIdx.end() );
  for(long int ii=0; ii < toProcessIdx.size()*parameters.percentagetoPreferredOrientation; ii++){
    

//   ThetaFactor=math.cos( math.radians( preferredAngles[0] % 181.0 )/(math.pi) )
//   Theta=np.random.normal((preferredAngles[1]%180.0)*ThetaFactor,variance)

      
 /*
  FIRST ATTEMPT, NO GOOD BUT WORKS
    double phiAngleDeg=parameters.eulerPreferredOrientationPhi+rnd_gaus(0, parameters.sigmaPreferredOrientation);
//    double phiAngleRad=    phiAngleDeg*PI/180.0;
    
//    double thetaAngleDeg = parameters.eulerPreferredOrientationTheta+rnd_gaus(0, parameters.sigmaPreferredOrientation);
      double thetaAngleDeg = 90+rnd_gaus(0, parameters.sigmaPreferredOrientation);

      //    double thetaAngleRad = thetaAngleDeg*PI/180.0;
      
    rlnPhi[toProcessIdx[ii]]=phiAngleDeg;// *(1-cos(thetaAngleRad));
    rlnTheta[toProcessIdx[ii]]=thetaAngleDeg;// *cos(thetaAngleRad);
*/
      
      //rlnPhi[toProcessIdx[ii]]=phiAngleDeg;//*(1-cos(thetaAngleRad));
      //rlnTheta[toProcessIdx[ii]]=thetaAngleDeg;//*cos(thetaAngleRad);
      
//      std::cerr<<"error Phi="<<parameters.eulerPreferredOrientationPhi<<"   theta="<<parameters.eulerPreferredOrientationTheta<<"\n";
      std::vector<double> errorAngle=induceEulerError(parameters.eulerPreferredOrientationPhi, parameters.eulerPreferredOrientationTheta, parameters.sigmaPreferredOrientation, parameters.sigmaPreferredOrientation);


//      int numGlobalOccurrences=10;
//      int numLocalOccurrences=10;
//      int startLocalAfterSuccesfulIte=3;
     // std::vector<double> errorAngle=induceGlobalEulerError(numOccurrences);


      
      rlnPhi[toProcessIdx[ii]]=errorAngle[0];
      rlnTheta[toProcessIdx[ii]]=errorAngle[1];
      
      
    //std::cerr<<" " << toProcessIdx[ii]<< "   ->   "<<rlnPhi[toProcessIdx[ii]] <<"\n";
  }

        std::ofstream fileOutput;
        fileOutput.open(parameters.starFileOut);
        fileOutput.close();
        fileOutput.open (parameters.starFileOut, std::ofstream::out | std::ofstream::app);   
        long int startMicrograph=getStarStart(parameters.inputFileToPreferredOrientation);
        std::ifstream fileParticles(parameters.inputFileToPreferredOrientation);
        std::string strLine;
        for (int counter=0;counter<startMicrograph;counter++){
           std::getline(fileParticles, strLine);
           fileOutput << strLine <<"\n";
        }
        for (unsigned long int ii=0, counter=0; ii<rlnPhi.size(); ii++){
	     std::getline(fileParticles, strLine);
	     std::string tmpString=strLine;
             while(!tmpString.empty() && std::isspace(*tmpString.begin()))
               tmpString.erase(tmpString.begin());
             tmpString=replaceValueStrlineStarFile(tmpString, phiIdx, std::to_string(rlnPhi[ii]));
             tmpString=replaceValueStrlineStarFile(tmpString, thetaIdx, std::to_string(rlnTheta[ii]));
	     fileOutput << tmpString <<"\n";

	}

	fileOutput.close();

 }



}

