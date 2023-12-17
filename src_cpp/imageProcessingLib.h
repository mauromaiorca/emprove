#ifndef __IMAGEPROCESSING_LIB__
#define __IMAGEPROCESSING_LIB__

#include <math.h>
#include <iostream>
//#include <iomanip>      // std::setprecision
//#include <complex>
#include <stdio.h>
#include <stdlib.h>
//#include <time.h>
#include <cstdlib>
#include "scores.h"



#include "mrcIO.h"


// ***********************************
//  
// 
template <typename T>
void blurMaskEdgeOnImage(T* I, T* maskI, T*imageOut, double sigma, int kernelLength , const unsigned long int nx, const unsigned long int ny){

//std::cerr<<"sigma="<< sigma << "   kernelLengthint="<< kernelLength<<"\n";

const int nxy = nx * ny;
const int length=kernelLength;
const int nxK =  (2*length+1);
const int nyK = (2*length+1);
const int kernelSize = nxK*nyK;
float * kernel = new float [ kernelSize ];
float * tmpI = new float [ nxy ];
//double angpix = 1;


double sum = 0;
double counter=0;

for (int jjT=-length; jjT<length+1; jjT++){
  for(int iiT=-length; iiT<length+1; iiT++){
          counter++;
        double distanceSquared =  pow(jjT,2.0)+pow(iiT,2.0);
        int X=iiT+length;
        int Y=jjT+length;
        kernel [  X + (Y * nxK) ] = 1 / (sigma * sigma * 2*3.1415927) 
              * exp (- (distanceSquared) / (2.0 * sigma * sigma));
        sum+=kernel [  X + (Y * nxK) ];
  }
}

//correct compensate for discretization (sum needs to be 1):
double normSum = (double((1.0-sum)))/double(kernelSize);
for (int ii=0; ii<kernelSize; ii++){
        kernel [ii] += normSum;
}


sum = 0;
for (int jjT=-length; jjT<length+1; jjT++){
  for(int iiT=-length; iiT<length+1; iiT++){
        //double distanceSquared =  pow(jjT,2.0)+pow(iiT,2.0);
        int X=iiT+length;
        int Y=jjT+length;
        sum+=kernel [ X +Y * nxK ];
 //       std::cerr << "kernel ["<< X +Y * nxK << "]=";
 //       std::cerr << kernel [ X +Y * nxK ] << " ";
  }
 // std::cerr<< " \n";
}
//std::cerr << " \n";
//std::cerr << "sum=" << sum <<" \n";



for (unsigned long int YY=length; YY<ny-length; YY++){
 for(unsigned long int XX=length; XX<nx-length; XX++){
   int numBlack=0;
   int numWhite=0;

   for (int jjT=-length; jjT<length; jjT++){
     for(int iiT=-length; iiT<length; iiT++){
        unsigned long int idx = (XX+iiT)+(YY+jjT)*nx;
        if ( maskI[idx] > 0.5f){
                numWhite++;  
        }else{
                numBlack++; 
        }
     }
   }


   //
   if (numBlack>0 && numWhite>0){
        double sum = 0;
        for (int jjT=-length; jjT<length; jjT++){
          for(int iiT=-length; iiT<length; iiT++){
                unsigned long int idx = (XX+iiT)+(YY+jjT)*nx;
                unsigned long int kernelIdx=iiT+length+(jjT+length)*nxK;
                sum+= kernel[kernelIdx]*I[idx];
          }
        }
        tmpI [ XX + YY * nx ] = sum;
   }else{
     tmpI [ XX + YY * nx ] = I [ XX + YY * nx ];
   }
 }
 

}

 for (unsigned long int ii=0;ii<(unsigned long int)nxy;ii++){
   imageOut[ii]= tmpI[ii];
 }

delete [] kernel;
delete [] tmpI;
}


// ************************
//
//      blurMaskedImageEdge2D
//
template <typename T, typename U>
void blurMaskedImageEdge2D(T* imageIn, U* IMask, T* imageOut, unsigned long int nx, unsigned long int ny,  std::vector<std::vector<double> > listImageEdgesSigmaWidthList, double maskThreshold = 0.9){
  unsigned long int nxy = nx * ny; 
  if (listImageEdgesSigmaWidthList.size()==0){
    //do nothing
    for (unsigned long int ij=0; ij<nxy; ij++){
      imageOut[ij]=imageIn[ij];
    }
  }else{
    for (unsigned int hh=0; hh<listImageEdgesSigmaWidthList.size(); hh++){
      if (hh==0 ){
        if ( listImageEdgesSigmaWidthList[0][0]==0 ){
            //do nothing
            for (unsigned long int ij=0; ij<nxy; ij++){
              imageOut[ij]=imageIn[ij];
            }
          }else{
            blurMaskEdgeOnImage(imageIn, IMask, imageOut, listImageEdgesSigmaWidthList[0][0], listImageEdgesSigmaWidthList[0][1], nx, ny, maskThreshold);
          }
        }else{
            blurMaskEdgeOnImage(imageOut, IMask, imageOut, listImageEdgesSigmaWidthList[hh][0], listImageEdgesSigmaWidthList[hh][1], nx, ny, maskThreshold);
        }
    }
  }
}



// *************************************
//
/** This function recursively implements a multiple Otsu thresholding algorithm for
 * multiple thresolds. The code is based on the paper:
 * Liao, P-S. & Chung, P-C. (September 2001), "A fast algorithm for multilevel thresholding", Journal of Information Science and Engineering 17 (5): 713-727
 */
//
// **************************************
void computeThresholdsRecursively(const double * H, //!< lookup table for the histogram
    double * t, //!< vector with tcomputed hreshold values
    const unsigned long int bins, //!< histogram bin size
    const unsigned long int NT, //!< number of threshold
    bool firstCall = true, int start=1, int end=1,  unsigned long int n=1, double * maxSig = NULL, double * Sq=NULL, unsigned int * m=NULL){

		if(firstCall){
			start=1;
			end=bins-NT;
			n=1;
			maxSig =new double (0);
			Sq=new double[NT+1];
			m=new unsigned int[NT+2];
		}

		for (int i = start; i < end; i++) {
			m[n]=i;
			if(n==1) {
				Sq[0] = H[1+bins*i];
			}
			if(n>1 && n <= NT){
				Sq[n-1] = H[start+bins*i];
			}
			if(n==NT){
				Sq[n] = H[i+1+bins*(bins-1)];
				double SqVal = 0;
				for (unsigned long int kk=0; kk <= NT; kk++){
					SqVal += Sq[kk];
				}
				//std::cerr<< "Sq="<<SqVal <<"\n";
				if (*maxSig < SqVal)	{
					for (unsigned long int cc=0; cc<NT; cc++){
							t[cc] = m[cc+1];
					}
					*maxSig = SqVal;
				}

			}else{
				computeThresholdsRecursively(H,t, bins, NT, false, i+1,end+1, n+1, maxSig, Sq, m);
			}
		}

		if(firstCall){
			delete maxSig;
			delete [] Sq;
			delete [] m;
		}

}



// *************************************
/** Compute Otsu threshold in an image */
// computeOtsuThresholds
// **************************************
template<typename T>
void computeOtsuThresholds ( T * image, //!< input image
    double * t, //!< output vector of threshold values
    const unsigned long int nt, //!< number of threshold values
    const unsigned long int bins, //!< number of bins in the histogram
    const unsigned long int nxyz, //!< size of the image
    bool verbose = false  //!< verbose output
){
		double epsilon = 0.01;

		// ----------------------------
		// COMPUTE HISTOGRAM
		// std::cerr<<" COMPUTE HISTOGRAM\n";
	  	double * histogram = new double [bins];
	  	buildHistogram(image, histogram,  bins, nxyz);
		T minImageValue=image[0], maxImageValue=image[0];
		for (unsigned long int ii = 1; ii< nxyz; ii++){
			if (minImageValue>image[ii]) minImageValue=image[ii];
			if (maxImageValue<image[ii]) maxImageValue=image[ii];
		}
		for (unsigned long int ii = 0; ii< bins; ii++){
			histogram[ii]=0.0f;
		}
		double ratioVal = 0;
		if (maxImageValue > minImageValue){
		 ratioVal=((double)((double)bins-1.0f))/((double)((double)maxImageValue-(double)minImageValue));
		}
		for (unsigned long int ii = 0; ii< nxyz; ii++){
			double val = (image[ii] - minImageValue)*ratioVal;
                        int index = floor(val);
			histogram[index]++;
		}

		//std::cerr<<"Histogram Statistics: min value="<< minImageValue <<"   max value=" <<  maxImageValue <<"\n";


		// ----------------------------
		//BUILD LOOKUP TABLES
		//std::cerr<<" BUILD LOOKUP TABLES\n";
	  	double * P = new double[bins*bins];
	  	double * S = new double[bins*bins];
	  	double * H = new double[bins*bins];
			for (unsigned long int i = 0; i < bins*bins; i++){
				P[i]=S[i]=H[i]=0.0f;
			}
			//diagonal (row 0 is all zero)
			//std::cerr<<" diagonal\n";
			for (unsigned long int i=1; i < bins; i++){
				unsigned long int diagonalIndex = i+i*bins;
				 P[diagonalIndex] = histogram[i];
				 S[diagonalIndex] = ((double) i)*histogram[i];
			}
			// calculate first row (row 0 is all zero)
			//std::cerr<<" calculate first row\n";
			//unsigned long int firstRowIndex = bins;
			for (unsigned long int i=1; i < bins-1; i++) {
			  P [1+(i+1)*bins] = P[1+i*bins] + histogram[i+1];
			  S [1+(i+1)*bins] = S[1+i*bins] + ((double) (i+1))*histogram[i+1];
			}
			// using row 1 to calculate others
			//std::cerr<<" using row 1 to calculate others\n";
			for (unsigned long int i=2; i < bins; i++){
			  for (unsigned long int j=i+1; j < bins; j++) {
				P[i+j*bins] = P[1+j*bins] - P[1+(i-1)*bins];
				S[i+j*bins] = S[1+j*bins] - S[1+(i-1)*bins];
			  }
			}
			// calculate H
			//std::cerr<<" calculate H\n";
			for (unsigned long int i=1; i < bins; i++){
			  for (unsigned long int j=i+1; j < bins; j++){
			  		unsigned long int ijIndex = i+j*bins;
					if (pow (P[ ijIndex ],2) > epsilon)
					  H[ ijIndex ] = (S[ ijIndex ]*S[ ijIndex ])/P[ ijIndex ];
					else
					  H[ ijIndex ] = 0.0f;
			  }
			}
	  	delete [] P;
	  	delete [] S;


		// ----------------------------
		// COMPUTE THRESHOLD VALUES
		if (verbose) std::cerr<<"[COMPUTE THRESHOLD VALUES...";
		computeThresholdsRecursively(H, t, bins, nt);
		//computeThresholdsIteratively(H, t, bins, nt);
		if (verbose) std::cerr<<" DONE]\n";
	  	delete [] H;
	  	delete [] histogram;

		//std::cerr<<"Binned Thresholds=[";
		//for (int i =0; i<nt; i++){
			//std::cerr<<t[i];
			//if (i<nt-1)
				//std::cerr<<",";
		//}
		//std::cerr<<"]\n";

		// ----------------------------
		// NORMALIZE HISTOGRAM VALUES
		// TO ORIGINAL IMAGE VALUES
		//std::cerr<<"[NORMALIZE HISTOGRAM VALUES...";
		double ratioValInv = 0;
		if (maxImageValue > minImageValue && bins>1){
		 ratioValInv=(double)((double)maxImageValue-(double)minImageValue)/(double)((double)bins-1.0f);
		}
		for(unsigned int tt = 0; tt < nt; tt++ ){
			t [ tt ] = minImageValue + t [ tt ]*ratioValInv;
		}
		//std::cerr<<" DONE]\n";
}






template<typename T, typename U>
void maskedThreshold ( T * image, //!< input image
    U * mask, //!< mask image
    U * maskOut, //!< mask image
    const unsigned long int nx, //!< size of the image
    const unsigned long int ny, //!< size of the image
    const unsigned long int nz, //!< size of the image
    const int numThresholds = 2
){
  unsigned long int nxyz = nx * ny * nz;
  unsigned long int counter = 0;
  for (unsigned long int ijk =0 ; ijk< nxyz; ijk++){
    maskOut[ijk]=0;
    if ( mask[ijk] > 0.001 ){
      counter++;
    }
  }
  float * values = new float [counter];

  for (unsigned long int ijk =0, tmpCounter = 0 ; ijk < nxyz; ijk++){
    if ( mask[ijk] > 0.001 ){
      values[tmpCounter]=image[ijk];
      tmpCounter++;
    }
  }

  //const int numThresholds = 2;
  double * threshold= new double [numThresholds];
  int bins = 40;
  computeOtsuThresholds ( values,  threshold, numThresholds, bins, counter);
  //std::cerr<<"threshold="<<threshold<<"\n";
  for (unsigned long int ijk =0 ; ijk < nxyz; ijk++){
    //maskOut[ijk] = image[ijk];
    if ( mask[ijk] > 0.001 ){
      if (image[ijk] > threshold[numThresholds-1]){
        maskOut[ijk] = 1;
      }
    }
  }
  delete [] values;
  delete [] threshold;
}





template<typename T, typename U>
std::vector<double> computeMaskedAnalysis ( T * image, //!< input image
    U * mask, //!< mask image
    const unsigned long int nx, //!< size of the image
    const unsigned long int ny, //!< size of the image
    const unsigned long int nz //!< size of the image
){
  std::vector<double> result;
  unsigned long int nxyz = nx * ny * nz;
  unsigned long int counter = 0;
  for (unsigned long int ijk =0 ; ijk< nxyz; ijk++){
    if ( mask[ijk] > 0.001 ){
      counter++;
    }
  }
  if (counter < 2) {
    result.push_back(-1);
    result.push_back(-1);
    return result;
  }  
  float * values = new float [counter];

  long double average=0;
  unsigned long int tmpCounter = 0;
  for (unsigned long int ijk =0; ijk < nxyz; ijk++){
    if ( mask[ijk] > 0.001 ){
      values[tmpCounter]=image[ijk];
      average+=image[ijk];
      tmpCounter++;
    }
  }
  average/=(long double)(tmpCounter);


  long double averageMin=0;
  unsigned long int tmpCounterMax = 0;
  long double averageMax=0;  
  unsigned long int tmpCounterMin = 0;

  const int numThresholds = 2;
  double * threshold= new double [numThresholds];
  int bins = 40;
  computeOtsuThresholds ( values,  threshold, numThresholds, bins, counter);
  //std::cerr<<"threshold="<<threshold<<"\n";
  for (unsigned long int ijk =0 ; ijk < nxyz; ijk++){
    //maskOut[ijk] = image[ijk];
    if ( mask[ijk] > 0.001 ){
      if (image[ijk] > threshold[numThresholds-1]){
        averageMin+=image[ijk];
        tmpCounterMax++;
      }else{
        averageMax+=image[ijk];
        tmpCounterMin++;
      }
    }
  }

  if (tmpCounterMin>1){
    averageMin/=(long double)(tmpCounterMin);
    result.push_back(averageMin);
  }else{
    result.push_back(-1);
  }
  if (tmpCounterMax>1){
    averageMax/=(long double)(tmpCounterMax);
    result.push_back(averageMax);
  }else{
    result.push_back(-1);
  }


  delete [] values;
  delete [] threshold;
  return result;
}


template<typename T>
std::vector<double> rotAverage ( T * image, //!< input image
    const unsigned long int nx, //!< size of the image
    const unsigned long int ny, //!< size of the image
    const unsigned long int nz //!< size of the image
){
  std::vector<double> returnVector;
  unsigned long int nxyz=nx*ny*nz;
  //T * idxI=new T [nxyz];

  //  get the size of the buffer
  unsigned long int sizeBuf=0;
  double nx0=nx/2.0;
  double ny0=ny/2.0;
  double nz0=ny/2.0;
  for (unsigned long int kk=0; kk<nz; kk++){
    double distanceZ=pow((double)kk-nz0,2.0);
    for (unsigned long int jj=0; jj<ny; jj++){
      double distanceY=pow((double)jj-ny0,2.0);
      for (unsigned long int ii=0; ii<nx; ii++){
        double distanceX=pow((double)ii-nz0,2.0);
        unsigned long int distance= floor(pow( distanceX+distanceY+distanceZ,0.5 ));
        if ( distance > sizeBuf ) sizeBuf = distance;

      }
    }
  }
  sizeBuf++;
  double * bufferI= new double [sizeBuf];
  double * bufferCounterI= new double [sizeBuf];
  for (unsigned long i=0; i<sizeBuf; i++){
    bufferI[i]=bufferCounterI[i]=0;
  }


  for (unsigned long int kk=0; kk<nz; kk++){
    double distanceZ=pow((double)kk-nz0,2.0);
    for (unsigned long int jj=0; jj<ny; jj++){
      double distanceY=pow((double)jj-ny0,2.0);
      for (unsigned long int ii=0; ii<nx; ii++){
        double distanceX=pow((double)ii-nz0,2.0);
        unsigned long int distance= floor(pow( distanceX+distanceY+distanceZ,0.5 ));
        unsigned long int idx =ii+jj*nx+kk*(nx*ny);
        //idxI[idx] = distance;
        if (distance<sizeBuf ){
          bufferI[distance]+=distance;
          bufferCounterI[distance]++;
        }
      }
    }
  }
  for (long int i=0; i<sizeBuf; i++){
    if ( bufferCounterI[i] > 0 ){
      returnVector.push_back(bufferI[i]/bufferCounterI[i]);
    }
  }

  //writeMrcImage("idx.mrc", idxI,  nx,  ny,  nz,  1);
  //delete [] idxI;
  delete [] bufferI;
  delete [] bufferCounterI;
  return returnVector;
}








template<typename T>
void amplitudeNormalization ( T * AmplitudesImage, //!< input image
    T * amplitudesToReplace, //!< input image
    const unsigned long int nx, //!< size of the image
    const unsigned long int ny, //!< size of the image
    const unsigned long int nz //!< size of the image
){

  std::vector<double> meanAmplitudesIn=rotAverage(AmplitudesImage, nx, ny, nz);
  std::vector<double> meanAmplitudesOut=rotAverage(amplitudesToReplace, nx, ny, nz);
  std::vector<double> normalizingVector;
  const long double thresholdLowAmplitude = 0.0000000001;
  for (unsigned long int ii=0; ii<meanAmplitudesIn.size(); ii++){
    double normValue=meanAmplitudesIn[ii];
    if ( meanAmplitudesIn[ii] > thresholdLowAmplitude ){
      normValue = meanAmplitudesOut[ii]/meanAmplitudesIn[ii];
    }
    normalizingVector.push_back(normValue);
  }

  unsigned long int nxyz=nx*ny*nz;
  unsigned long int sizeBuf=meanAmplitudesIn.size();
  double nx0=nx/2.0;
  double ny0=ny/2.0;
  double nz0=ny/2.0;


  for (unsigned long int kk=0; kk<nz; kk++){
    double distanceZ=pow((double)kk-nz0,2.0);
    for (unsigned long int jj=0; jj<ny; jj++){
      double distanceY=pow((double)jj-ny0,2.0);
      for (unsigned long int ii=0; ii<nx; ii++){
        double distanceX=pow((double)ii-nz0,2.0);
        unsigned long int distance= floor(pow( distanceX+distanceY+distanceZ,0.5 ));
        unsigned long int idx =ii+jj*nx+kk*(nx*ny);
        //idxI[idx] = distance;
        
        if ( distance < sizeBuf ){
          AmplitudesImage[idx] = AmplitudesImage[idx] * normalizingVector[distance];
        }
      }
    }
  }
  //writeMrcImage("idx.mrc", idxI,  nx,  ny,  nz,  1);
  //delete [] idxI;
}



#endif
