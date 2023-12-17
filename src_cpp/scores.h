#ifndef __SCORING_FUNCTIONS_LIB___
#define __SCORING_FUNCTIONS_LIB___

#include <math.h>
#include <iostream>
#include <vector>
#include "derivative_libs.h"

//Normalize
std::vector<double> normalize(float * image1, //!< input image
        float * image2, //!< input image
        float * mask,
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        double maskThreshold = 0.1
){

  std::vector<double> result;

  float max1=0;
  float min1=0;
  float max2=0;
  float min2=0;
  double mean1=0;
  double mean2=0;
  long int counter = 0;
  for (unsigned long int ijk = 0; ijk< nxyz; ijk++){
  if (mask[ijk]>maskThreshold){
    if (counter==0){
      max1=image1[ijk];
      max2=image2[ijk];
      min1=image1[ijk];
      min2=image2[ijk];
    }else{
      if (max1<image1[ijk]) max1=image1[ijk];
      if (max2<image2[ijk]) max2=image2[ijk];
      if (min1>image1[ijk]) min1=image1[ijk];
      if (min2>image2[ijk]) min2=image2[ijk];
    }
    mean1+=image1[ijk];
    mean2+=image2[ijk];
    counter++;
   }
  }
  if (counter>0){
    mean1/=counter;
    mean2/=counter;
  }
  double diffI1 = (double)(max1-min1);
  double diffI2 = (double)(max2-min2);
  
  result.push_back(mean1);
  result.push_back(diffI1);
  result.push_back(mean2);
  result.push_back(diffI2);
  return result;

}


// *************************************
/** build the grayscale histogram from an image */
// buildHistogram
// **************************************
template<typename T>
void buildHistogram(const T* image, //!< input image
        double * histogram, //!< output histogram
        const unsigned long int bins, //!< number of bins for the histogram
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        const T * mask = NULL,
        double minMaskValue = 0.0
){
  if (mask==NULL){ //NO MASK
	double min=image[0], max=image[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
		if (min>image[ii]) min=image[ii];
		if (max<image[ii]) max=image[ii];
	}
    for (unsigned long int ii = 0; ii< bins; ii++){
		histogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
		double val = 0;
		if (max-min>0){
		 val = (image[ii] - min)*(bins-1)/(max-min);
		}
		int index = floor(val);
		histogram[index]++;
	}
   }else{ //MASK
	double min=image[0], max=image[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
	   if(mask[ii]>minMaskValue){
		if (min>image[ii]) min=image[ii];
		if (max<image[ii]) max=image[ii];
           }
	}
        for (unsigned long int ii = 0; ii< bins; ii++){
		histogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
	   if(mask[ii]>minMaskValue){
		double val = 0;
		if (max-min>0){
		 val = (image[ii] - min)*(bins-1)/(max-min);
		}
		int index = floor(val);
		histogram[index]++;
           }
	}
   }
 }



// *************************************
//  SSIM
// **************************************
template<typename T>
double SSIM(const T* GT, const T* I2, unsigned long int nxyz, const T* mask=NULL, double minMaskValue = 0.1){

  const T* I1 = GT;
  double meanI1=0;
  double meanI2=0;
  double sdI1=0;
  double sdI2=0;
  double sigmaI12=0;
  unsigned long int counter  = 0;
  double min = 0;
  double max = 1;
  double dynamicRange=1;

  if (mask){
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
        if ( mask[ijk] > minMaskValue ){
            if (counter==0){
                meanI1=I1[ijk];
                meanI2=I2[ijk];
            }else{
                meanI1+=I1[ijk];
                meanI2+=I2[ijk];
                if (max<GT[ijk]) max=GT[ijk];
                if (min>GT[ijk]) min=GT[ijk];

            }
            counter++;
        }
    }
  }else{
    counter = nxyz;
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
        meanI1+=I1[ijk];
        meanI2+=I2[ijk];
        if (max<GT[ijk]) max=GT[ijk];
        if (min>GT[ijk]) min=GT[ijk];
    }
  }
  if (max-min>1)
    dynamicRange=max-min;


  if (counter > 0){
      meanI1/=counter;
      meanI2/=counter;
  }
  if (mask){
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
      if ( mask[ijk] > minMaskValue ){
        sdI1+=pow(I1[ijk]-meanI1,2.0);
        sdI2+=pow(I2[ijk]-meanI2,2.0);
        sigmaI12+=(I1[ijk]-meanI1)*(I2[ijk]-meanI2);
      }
    }
  }else{
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
      sdI1+=pow(I1[ijk]-meanI1,2.0);
      sdI2+=pow(I2[ijk]-meanI2,2.0);
      sigmaI12+=(I1[ijk]-meanI1)*(I2[ijk]-meanI2);
    } 
  }

  if (counter > 1.0){
      sdI1=pow((1.0/(counter-1.0))*counter,0.5);
      sdI2=pow((1.0/(counter-1.0))*counter,0.5);
      sigmaI12=pow((1.0/(counter-1.0))*sigmaI12,0.5);
  }

  //C1=(K1*L)^2, where L is the dynamic range of the pixel values, and K1 << 1
  dynamicRange=1.0;
  double C1 = pow(0.00000001*dynamicRange,2);
  double C2 = pow(0.00000001*dynamicRange,2);
  double C3 = C2/2.0;
  double luminance=(2*meanI1*meanI2+C1)/(meanI1*meanI1+meanI2*meanI2+C1);
  double contrast=(2*sdI1*sdI2+C2)/(sdI1*sdI1+sdI2*sdI2+C2);
  double structure=(sigmaI12+C3)/(sdI1*sdI2+C3);
  double ssim=luminance*luminance + contrast*contrast + structure*structure;
  return ssim;
}



// *************************************
//  PSNR
// **************************************
template<typename T>
double PSNR(const T* GT, const T* I2, unsigned long int nxyz, const T* mask=NULL, double minMaskValue = 0.1){

  double maxGT=0;
  double MSE=0;
  unsigned long int counter  = 0;

  if (mask){
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
        if ( mask[ijk] > minMaskValue ){
            if (counter==0){
                maxGT=GT[ijk];
            }else{
                if ( maxGT<GT[ijk] ){
                  maxGT=GT[ijk];
                }
                MSE+=pow(GT[ijk]-I2[ijk], 2.0);
            }
            counter++;
        }
    }
  }else{
    counter = nxyz;
    for (unsigned long int ijk=0; ijk<nxyz;ijk++){
      if ( maxGT<GT[ijk] ){
        maxGT=GT[ijk];
      }
      counter++;
      MSE+=pow(GT[ijk]-I2[ijk], 2.0);
    }
  }
//  std::cerr<<"counter="<<counter<<"\n";
//  std::cerr<<"maxGT="<<maxGT<<"\n";
//  std::cerr<<"MSE_unnormalized="<<MSE<<"\n";
  MSE=(1.0/counter)*MSE;
//  std::cerr<<"MSE="<<MSE<<"\n";
  double PSNR_val=20.0*log10(maxGT/pow(MSE,0.5));
//  std::cerr<<"PSNR_val="<<PSNR_val<<"\n";
  if (PSNR_val != PSNR_val){
    PSNR_val = 0;
  }

  return PSNR_val;
}



// *************************************
//  crossCorrelationDistance
// **************************************
template<typename T>
double crossCorrelationDistance(const T* I1, const T* I2, unsigned long int nxyz, const T* mask=NULL, double minMaskValue = 0){
 double mean1 = 0;
 double mean2 = 0;
 double totalValues = 0;
 if (mask==NULL){
  for (unsigned long int ii = 0; ii<nxyz; ii++){
   mean1+=I1[ii];
   mean2+=I2[ii];
  }
  totalValues=nxyz;
 }else{
  for (unsigned long int ii = 0; ii<nxyz; ii++){
   if(mask[ii]>minMaskValue){
     mean1+=I1[ii];
     mean2+=I2[ii];
     totalValues++;
   }
  }
 }
 if ( totalValues > 0 ){
    mean1/=totalValues;
    mean2/=totalValues;
 }
 double numer=0;
 double denom1=0;
 double denom2=0; 
 if (mask==NULL){
         for (unsigned long int ii = 0; ii<nxyz; ii++){
          numer+=(I1[ii]-mean1)*(I2[ii]-mean2);
          denom1+=pow(I1[ii]-mean1,2);
          denom2+=pow(I2[ii]-mean2,2);
         }
 }else{
         for (unsigned long int ii = 0; ii<nxyz; ii++){
          if(mask[ii]>minMaskValue){
                  numer+=(I1[ii]-mean1)*(I2[ii]-mean2);
                  denom1+=pow(I1[ii]-mean1,2);
                  denom2+=pow(I2[ii]-mean2,2);
          }
         }
 }
 
 double CCC=0;
 double denomin=pow(denom1,0.5)*pow(denom2,0.5);
 
 if (numer==0.f){
   CCC=0;
 }else if (pow(denomin,2)<0.0000000001f){
  CCC=1;
 }else{
   CCC=numer/denomin;
 }
 return CCC;
 
}




// *************************************
//  SCI: structural Cross Correlation Index
// **************************************
template<typename T>
double SCI( T* I1, 
             T* I2, 
            unsigned long int nx, unsigned long int ny, unsigned long int nz, 
            double sigma, 
            const T* mask=NULL, 
            double minMaskValue = 0.1){
 
  unsigned long int nxyz = nx * ny * nz;
  double angpix = 1;
  int padding = 2;
  //compute first and second derivatives
  T * I1_proc= new T [ nxyz ];
  T * I2_proc= new T [ nxyz ];

double CC = crossCorrelationDistance(I1, I2,  nxyz, mask, minMaskValue);
if (CC != CC){
  CC=0;
} else if ( CC < 0 ){
  CC=0;
}
//std::cerr<<"CC="<<CC<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I2,	I2_proc);
  double CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //  std::cerr<<"first derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"first derivative Y="<<CC_tmp<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"second derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"score full="<<CC<<"\n";
  //std::cerr<<"second derivative Y="<<CC_tmp<<"\n";

  delete [] I1_proc;
  delete [] I2_proc;

  return CC;
}




// *************************************
//  SCI: Structural Difference Image 
// **************************************

template<typename T>
double SD_norm ( T* I1, 
              T* I2, 
              T* OutI, unsigned long int nxyz){

  long double mean1=0;
  long double mean2=0;
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    mean1+=I1[ijk];
    mean2+=I2[ijk];
  }
  mean1/=nxyz;
  mean2/=nxyz;

  double mean = 0;
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]=exp(-0.1*pow( I1[ijk]-mean1 -I2[ijk]+mean2, 2.0));
  }
//  mean /= (double)nxyz;
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]/=5.0;
  }
}


// *************************************
//  SCI: Structural Difference Image 
// **************************************
template<typename T>
double SDIM ( T* I1, 
              T* I2, 
              T* OutI,
              unsigned long int nx, unsigned long int ny, unsigned long int nz, 
              double sigma, 
              const T* mask=NULL, 
              double minMaskValue = 0.1 ) {
 
  unsigned long int nxyz = nx * ny * nz;
  double angpix = 1;
  int padding = 2;
  //compute first and second derivatives
  T * I1_proc= new T [ nxyz ];
  T * I2_proc= new T [ nxyz ];
  T * tmpI= new T [ nxyz ];


  double CC = crossCorrelationDistance(I1, I2,  nxyz, mask, minMaskValue);
  if (CC != CC){
    CC=0;
  } else if ( CC < 0 ){
    CC=0;
  }
  SD_norm(I1,I2,tmpI,nxyz);  
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]=tmpI[ijk];
  }
  writeMrcImage("I0.mrc", tmpI, nx,ny,nz);
  //std::cerr<<"CC="<<CC<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I2,	I2_proc);
  double CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    //OutI[ijk]+=tmpI[ijk];
    OutI[ijk]+= tmpI[ijk]; //OK
  }
  writeMrcImage("I1x.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //  std::cerr<<"first derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]+=tmpI[ijk];
  }
  writeMrcImage("I1y.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"first derivative Y="<<CC_tmp<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]+=tmpI[ijk];
  }
  writeMrcImage("I1xx.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }
  //std::cerr<<"second derivative X="<<CC_tmp<<"\n";


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I2,	I2_proc);
  CC_tmp = crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  SD_norm(I1_proc,I2_proc,tmpI,nxyz);
  for (unsigned long int ijk=0; ijk<nxyz; ijk++){
    OutI[ijk]+=tmpI[ijk];
  }
  writeMrcImage("I1yy.mrc", tmpI, nx,ny,nz);
  if (CC_tmp != CC_tmp){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp;
  }


  delete [] I1_proc;
  delete [] I2_proc;
  delete [] tmpI;
}





// *************************************
//  SCI: structural Cross Correlation Index
// **************************************
template<typename T>
double SCI_sqr( T* I1, 
             T* I2, 
            unsigned long int nx, unsigned long int ny, unsigned long int nz, 
            double sigma, 
            const T* mask=NULL, 
            double minMaskValue = 0.1){
 
  double minTolerance = 0.0;
  unsigned long int nxyz = nx * ny * nz;
  double angpix = 1;
  int padding = 2;
  //compute first and second derivatives
  T * I1_proc= new T [ nxyz ];
  T * I2_proc= new T [ nxyz ];

  T * I1_procReverse= new T [ nxyz ];
  T * I2_procReverse= new T [ nxyz ];


double CC = crossCorrelationDistance(I1, I2,  nxyz, mask, minMaskValue);
if (CC != CC){
  CC=0;
} else if ( CC < 0 ){
  CC=0;
}
//std::cerr<<"CC="<<CC<<"\n";

  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	1,	 (T*) I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  //writeMrcImage("I2_proc.mrc", I2_proc, nx,ny,nz);
  //writeMrcImage("I2_procReverse.mrc", I2_procReverse, nx,ny,nz);
  double CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dx tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	1,	(T*)  I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dy tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 0,	2,	(T*)  I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dxx tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }


  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I1,	I1_proc);
  gaussRecursiveDerivatives1D (sigma, nx,	ny, nz,	angpix,	padding, 1,	2,	(T*)  I2,	I2_proc);
  for (unsigned long int ijk=0;ijk<nxyz;ijk++){
    I1_proc[ijk]=pow(I1_proc[ijk],2.0);
    I2_proc[ijk]=pow(I2_proc[ijk],2.0);
  }
  CC_tmp1 = minTolerance+crossCorrelationDistance(I1_proc, I2_proc,  nxyz, mask, minMaskValue);
  //std::cerr<<"dyy tmp1="<<CC_tmp1<<"\n";
  if (CC_tmp1 != CC_tmp1){
    CC = 0;
  } else if ( CC < 0 ){
    CC=0;
  }else{
    CC *= CC_tmp1;
  }

  delete [] I1_proc;
  delete [] I2_proc;
  delete [] I1_procReverse;
  delete [] I2_procReverse;

  return CC;
}



template<typename T>
double ssdDistance(const T* I1, const T* I2, unsigned long int nxyz, const T* mask=NULL){
 double SSD=0;
 if (mask==NULL){
    for (unsigned long int ii = 0; ii<nxyz; ii++){
       SSD+=pow(I1[ii]-I2[ii],2.0);
    }
 }else{
    for (unsigned long int ii = 0; ii<nxyz; ii++){
       if (mask[ii]>0){  
         SSD+=pow(I1[ii]-I2[ii],2.0);
       }
    }
 }
 
 SSD=pow(SSD,0.5);
 return SSD;
}

template<typename T>
double ssdDistance(std::vector<T> I1, std::vector<T> I2){
 double SSD=0;
 unsigned long int size = I1.size();
 for (unsigned long int ii = 0; ii<size; ii++){
  SSD+=pow(I1[ii]-I2[ii],2.0);
 }
 SSD=pow(SSD,0.5);
 return SSD;
}

template<typename T, typename U>
void squareDistanceImage(const T* I1, const T* I2, U* I12, unsigned long int nxyz, const T* mask=NULL){
 //double SD=0;
 if (mask==NULL){
   for (unsigned long int ii = 0; ii<nxyz; ii++){
     I12[ii]=pow(I1[ii]-I2[ii],2.0);
   }
 }else{
   for (unsigned long int ii = 0; ii<nxyz; ii++){
       if (mask[ii]>0){  
         I12[ii]=pow(I1[ii]-I2[ii],2.0);
       }else{
         I12[ii]=0;
       }
   }
 }
}

template<typename T>
double sdDistance(const T* I1, const T* I2, unsigned long int nxyz, const T* mask=NULL){
 double SD=0;
 if (mask==NULL){
   for (unsigned long int ii = 0; ii<nxyz; ii++){
     SD+=I1[ii]-I2[ii];
   }
 }else{
   for (unsigned long int ii = 0; ii<nxyz; ii++){
       if (mask[ii]>0){  
         SD+=I1[ii]-I2[ii];
       }
   }
 }
 return SD;
}

template<typename T>
double sdDistance(std::vector<T> I1, std::vector<T> I2){
 double SD=0;
 unsigned long int size = I1.size();
 for (unsigned long int ii = 0; ii<size; ii++){
  SD+=I1[ii]-I2[ii];
 }
 return SD;
}


// *************************************
//  buildJointHistogram
// **************************************
template<typename T>
void buildJointHistogram(const T* image1, //!< input image1
        const T* image2, //!< input image2
        double * jointHistogram, //!< output histogram
        const unsigned long int bins, //!< number of bins for the histogram
        const unsigned long int nxyz,  //!< size (in pixel) of the image
        const T* mask=NULL,
        double minMaskValue=0
){
  //std::cerr<<"starting HIST\n";
  if (mask==NULL){
	double min1=image1[0], max1=image1[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
		if (min1>image1[ii]) min1=image1[ii];
		if (max1<image1[ii]) max1=image1[ii];
	}
	double min2=image2[0], max2=image2[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
		if (min2>image2[ii]) min2=image2[ii];
		if (max2<image2[ii]) max2=image2[ii];
	}
        for (unsigned long int ii = 0; ii< bins*bins; ii++){
		jointHistogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
               	double val1 = 0;
               	double val2 = 0;
	        if (max1>min1)
        		val1 = (image1[ii] - min1)*(bins-1)/(max1-min1);
	        if (max2>min2)
        		val2 = (image2[ii] - min2)*(bins-1)/(max2-min2);
	        int index1 = floor(val1);
	        int index2 = floor(val2);
		jointHistogram[index1+bins*index2]++;
	}

  } else { // MASK
	double min1=image1[0], max1=image1[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
	   if (  mask[ii] > minMaskValue ){
		if (min1>image1[ii] ) min1=image1[ii];
		if (max1<image1[ii] ) max1=image1[ii];
	   }
	}
	double min2=image2[0], max2=image2[0];
	for (unsigned long int ii = 1; ii< nxyz; ii++){
	   if ( ((long int) mask[ii]) > (long int)minMaskValue ){
		if ( min2>image2[ii] ) min2=image2[ii];
		if ( max2<image2[ii] ) max2=image2[ii];
           }
	}
  for (unsigned long int ii = 0; ii< bins*bins; ii++){
		jointHistogram[ii]=0.0f;
	}
	for (unsigned long int ii = 0; ii< nxyz; ii++){
	   if ( ((long int) mask[ii]) > (long int)minMaskValue){
               	double val1 = 0;
               	double val2 = 0;
	        if (max1>min1)
        		val1 = (image1[ii] - min1)*(bins-1)/(max1-min1);
	        if (max2>min2)
        		val2 = (image2[ii] - min2)*(bins-1)/(max2-min2);
	        int index1 = floor(val1);
	        int index2 = floor(val2);
		jointHistogram[index1+bins*index2]++;
           }
	}
  }
  //std::cerr<<"ending HIST\n";
}


template<typename T, typename U>
std::vector<double> JointHistogramLinearRegression(const T* jointHistogram, //!< joint histogram input
                                    const int numOfBins,
                                    U* LineDistributionImage,
                                    int mainLineID,
                                    T* jointHistogramUp=NULL,
                                    T* jointHistogramDown=NULL,
                                    int verbosity = 1){
                                                                     
                                    
 std::vector<double> vectorOutput;

 //compute mean
 double meanX=0;
 double meanY=0;
 double counter=0;
 //fill the image
 for ( int ii=0;ii<numOfBins;ii++){
  for ( int jj=0;jj<numOfBins;jj++){
    double value = jointHistogram[ii+numOfBins*jj];
    if (value>0){
     meanX+=value*ii;
     meanY+=value*jj;
     counter+=value;
    }
  }
 }
 if (counter>0){
  meanX/=counter;
  meanY/=counter;
 }
 double covariance=0;
 double variance=0;
 for ( int ii=0;ii<numOfBins;ii++){
  for ( int jj=0;jj<numOfBins;jj++){
    double value = jointHistogram[ii+numOfBins*jj];
    if (value>0){
     covariance+=value*(ii-meanX)*(jj-meanY);
     variance+=value*pow(ii-meanX,2.0);
    }
  }
 }
 double slope=0;
 double intercept=0;
 if (variance>0){
  slope=covariance/variance;
 }
 intercept=meanY-slope*meanX;
 
 
 double slopePependicular=0;
 double interceptPerpendicular=0;
 if (slope!=0){
  slopePependicular=-1.0/slope;
 }
 interceptPerpendicular=meanY-slopePependicular*meanX;
 double angularDifference=(atan(1.0)-atan(slope))*180.0/3.14159265359;
 //slope=1;
 if (verbosity > 0 ){
   std::cerr<<"slope="<<slope<<"    intercept="<<intercept<<"\n";
   std::cout<<"angleDifference="<< angularDifference <<"\n";
 }
 vectorOutput.push_back(slope);
 vectorOutput.push_back(intercept);
 vectorOutput.push_back(angularDifference);
 
 //fill the image and the two up/down histograms
 
 for ( int ii=0;ii<numOfBins;ii++){
    int y=floor(intercept+ii*slope);
    if (y>0 && y<numOfBins){
     LineDistributionImage[ii+y*numOfBins]=mainLineID;
    }
    
    //FILL HISTOGRAMS UP AND DOWN, IF REQUESTED (NOT NULL)
    if (jointHistogramUp != NULL && jointHistogramDown!=NULL)
     for ( int jj=0;jj<numOfBins;jj++){
      if (y==jj) {
       jointHistogramUp[ii+jj*numOfBins]=0;
       jointHistogramDown[ii+jj*numOfBins]=0;
      }else{
       if (y>jj){
        jointHistogramUp[ii+jj*numOfBins]=jointHistogram[ii+jj*numOfBins];
        jointHistogramDown[ii+jj*numOfBins]=0;
       }else{
        jointHistogramUp[ii+jj*numOfBins]=0;
        jointHistogramDown[ii+jj*numOfBins]=jointHistogram[ii+jj*numOfBins];
       }
      }
     }


    //write this extra info only if you have histograms up and down to take care of
    if (jointHistogramUp != NULL && jointHistogramDown!=NULL){
     LineDistributionImage[ii+ii*numOfBins]=3;
    
     //perpendicular
     int yP=floor(interceptPerpendicular+ii*slopePependicular);    
     if (yP>0 && yP<numOfBins){
      LineDistributionImage[ii+yP*numOfBins]=2;
     }
    }
    
 }
 return vectorOutput;

}

// *************************************
//  computeEntropy
// **************************************
template<typename T>
double computeEntropy(const T* data1, unsigned long int nxyz, const T* mask=NULL){
 double entropy = 0;
 double sum=0;
 if (mask==NULL){
         for (unsigned long int ii=0; ii<nxyz; ii++){
          sum+=data1[ii];
         }
         if (sum == 0)
          return 0.0;
         for (unsigned long int ii=0; ii<nxyz; ii++){
          const double probability = data1[ii] / sum;
          if (probability > 0)
           entropy += - probability * log( probability ) / log( 2.0 );
         }
         return entropy;
 }else{ //MASK
         for (unsigned long int ii=0; ii<nxyz; ii++){
          if (mask[ii]>0){
            sum+=data1[ii];
          }
         }
         if (sum == 0)
          return 0.0;
         for (unsigned long int ii=0; ii<nxyz; ii++){
            if (mask[ii]>0){
                  const double probability = data1[ii] / sum;
                  if (probability > 0){
                   entropy += - probability * log( probability ) / log( 2.0 );
                  }              
            }
         }
         return entropy;
 }
}


// *************************************
//  NormalisedMutualInformation
// **************************************
template<typename T>
double NormalisedMutualInformation(const T* I1Target, 
                            const T* I2,
                            unsigned long int nxyz, 
                            unsigned long int numBins,
                            const T* mask=NULL,
                            double minMaskValue=0){

 double * histogram1 = new double [numBins];
 double * histogram2 = new double [numBins];
 double * histogram12 = new double [numBins*numBins];
 buildHistogram(I1Target, histogram1, numBins, nxyz,mask,minMaskValue);
 buildHistogram(I2, histogram2, numBins, nxyz,mask,minMaskValue);
 buildJointHistogram(I1Target, I2, histogram12, numBins, nxyz,mask,minMaskValue);
  
 double H1 = computeEntropy(histogram1, numBins);
 double H2 = computeEntropy(histogram2, numBins);
 double H12 = computeEntropy(histogram12, numBins*numBins);
 //std::cerr<<"H1="<<H1<<"     H2="<<H2<<"    H12="<<H12<<"\n";
 //double MutualInformation = H1+H2-H12;
 double NormalizedMutualInformation2 = 0;
 if (H12!=0) NormalizedMutualInformation2 = ( H1 + H2 ) / (2*H12);
 
 delete [] histogram1;
 delete [] histogram2;
 delete [] histogram12;
 return NormalizedMutualInformation2;

}



// *************************************
//  JointEntropyDistance
// **************************************
template<typename T>
double JointEntropyDistance(const T* I1Target, 
                            const T* I2,
                            unsigned long int nxyz, 
                            unsigned long int numBins){

 
 double * H12 = new double [numBins*numBins];
 double * H11 = new double [numBins*numBins];
 buildJointHistogram(I1Target, I2, H12, numBins, nxyz);
 buildJointHistogram(I1Target, I1Target, H11, numBins, nxyz);
 double distances = 0.0;
 double counter = 0.0001;
 double distance1 = 0.0;
 double distance2 = 0.0;
 
 for (unsigned long int ii = 0; ii<numBins; ii++){
  for (unsigned long int jj = 0; jj<numBins; jj++){
   if (H11[jj+ii*numBins] != 0){
    distance1+=H11[jj+ii*numBins];
    distance2+=H12[jj+ii*numBins];
    counter+=1;
   }
  }
 }

 distances = 0;
 if (distance1 != 0.0f ){
  distances=1.0-distance2/distance1;
 }
 delete [] H11;
 delete [] H12;
 return distances;


}



// *************************************
//  JointEntropyLocalDistance
// **************************************
template<typename T>
std::vector<float> imageEntropyDistance(const T* I1Target, const T* I2, T* resultI1, T* resultI2, unsigned long int nx, unsigned long int ny, unsigned long int nz, unsigned long int numBins){
 unsigned long int nxyz=nx*ny*nz;

 double * H12 = new double [numBins*numBins];
 double * H11 = new double [numBins*numBins];
 buildJointHistogram(I1Target, I2, H12, numBins, nxyz);
 buildJointHistogram(I1Target, I1Target, H11, numBins, nxyz);
 double distances = 0.0;
 double counter = 0.0001;
 double distance1 = 0.0;
 double distance2 = 0.0;
 
 double * probabilityVector1 = new double [numBins];
 double * probabilityVector2 = new double [numBins];
 for (unsigned long int ii = 0; ii<numBins; ii++){
  probabilityVector1[ii]=0;
  probabilityVector2[ii]=0;
 }
 
 for (unsigned long int ii = 0; ii<numBins; ii++){
  for (unsigned long int jj = 0; jj<numBins; jj++){
   if (H11[jj+ii*numBins] != 0){
    distance1+=H11[jj+ii*numBins];
    distance2+=H12[jj+ii*numBins];
    probabilityVector1[ii]+=H11[jj+ii*numBins];
    probabilityVector2[ii]+=H12[jj+ii*numBins];    
    //distances+=pow(H11[jj+ii*numBins]-H12[jj+ii*numBins],2);
    counter+=1;
   }
  }
 }
 distances=distance2/distance1;

 double * probabilityVector=&probabilityVector1[0];
 for (unsigned long int ii = 0; ii<numBins; ii++){
  if (probabilityVector1[ii]!=0){
   probabilityVector[ii]=probabilityVector2[ii]/probabilityVector1[ii];
  }else{
   probabilityVector[ii]=0;
  }
 }
 
 //
 double min1=I1Target[0], max1=I1Target[0];
 for (unsigned long int ii = 1; ii< nxyz; ii++){
  if (min1>I1Target[ii]) min1=I1Target[ii];
  if (max1<I1Target[ii]) max1=I1Target[ii];
 }
 double min2=I2[0], max2=I2[0];
 for (unsigned long int ii = 1; ii< nxyz; ii++){
  if (min2>I2[ii]) min2=I2[ii];
  if (max2<I2[ii]) max2=I2[ii];
 } 
 for (unsigned long int ii = 0; ii< nxyz; ii++){
  double val1 = 0;
  double val2 = 0;
  if (max1>min1)
   val1 = (I1Target[ii] - min1)*(numBins-1)/(max1-min1);
  if (max2>min2)
   val2 = (I2[ii] - min2)*(numBins-1)/(max2-min2);
  int index1 = floor(val1);
  int index2 = floor(val2);
  resultI1[ii]=pow(probabilityVector[index1]-probabilityVector[index2],2);
  resultI2[ii]=pow(I1Target[ii]-I2[ii],2);
 }

 double entropy1=computeEntropy(resultI1, nxyz);
 double entropy2=computeEntropy(resultI2, nxyz);
 std::vector<float> results;
 results.push_back(entropy1);
 results.push_back(entropy2);
 delete [] H11;
 delete [] H12;
 return results;
}




template<typename T>
double computeDerivatives (T* I1, T* IReproj, T* IMaskReproj, unsigned long int nx, unsigned long int ny, int depth, double sigma){
    const double _MIN_VAL=0.0f;//0.000001;
    const unsigned long int nxy = nx * ny;
    double CC1=crossCorrelationDistance(I1, IReproj, nxy, IMaskReproj);
    //double CC1 = NormalisedMutualInformation(I1, IReproj, nxy, 50, IMaskReproj); 

    //first
    T* I1_dx = new T [nxy];
    T* I1_dy = new T [nxy];
    T* IReproj_dx = new T [nxy];
    T* IReproj_dy = new T [nxy];

    //second
    T* I1_dxx = new T [nxy];
    T* I1_dxy = new T [nxy];
    T* I1_dyx = new T [nxy];
    T* I1_dyy = new T [nxy];
    T* IReproj_dxx = new T [nxy];
    T* IReproj_dxy = new T [nxy];
    T* IReproj_dyx = new T [nxy];
    T* IReproj_dyy = new T [nxy];


    //third
    T* I1_dxyx = new T [nxy];
    T* I1_dxxx = new T [nxy];
    T* I1_dyxx = new T [nxy];
    T* I1_dyyx = new T [nxy];
    T* I1_dxyy = new T [nxy];
    T* I1_dxxy = new T [nxy];
    T* I1_dyxy = new T [nxy];
    T* I1_dyyy = new T [nxy];

    T* IReproj_dxyx = new T [nxy];
    T* IReproj_dxxx = new T [nxy];
    T* IReproj_dyxx = new T [nxy];
    T* IReproj_dyyx = new T [nxy];
    T* IReproj_dxyy = new T [nxy];
    T* IReproj_dxxy = new T [nxy];
    T* IReproj_dyxy = new T [nxy];
    T* IReproj_dyyy = new T [nxy];


    //forth
    T* I1_dxyxx = new T [nxy];
    T* I1_dxxxx = new T [nxy];
    T* I1_dyxxx = new T [nxy];
    T* I1_dyyxx = new T [nxy];
    T* I1_dxyyx = new T [nxy];
    T* I1_dxxyx = new T [nxy];
    T* I1_dyxyx = new T [nxy];
    T* I1_dyyyx = new T [nxy];
    T* I1_dxyxy = new T [nxy];
    T* I1_dxxxy = new T [nxy];
    T* I1_dyxxy = new T [nxy];
    T* I1_dyyxy = new T [nxy];
    T* I1_dxyyy = new T [nxy];
    T* I1_dxxyy = new T [nxy];
    T* I1_dyxyy = new T [nxy];
    T* I1_dyyyy = new T [nxy];


    T* IReproj_dxyxx = new T [nxy];
    T* IReproj_dxxxx = new T [nxy];
    T* IReproj_dyxxx = new T [nxy];
    T* IReproj_dyyxx = new T [nxy];
    T* IReproj_dxyyx = new T [nxy];
    T* IReproj_dxxyx = new T [nxy];
    T* IReproj_dyxyx = new T [nxy];
    T* IReproj_dyyyx = new T [nxy];
    T* IReproj_dxyxy = new T [nxy];
    T* IReproj_dxxxy = new T [nxy];
    T* IReproj_dyxxy = new T [nxy];
    T* IReproj_dyyxy = new T [nxy];
    T* IReproj_dxyyy = new T [nxy];
    T* IReproj_dxxyy = new T [nxy];
    T* IReproj_dyxyy = new T [nxy];
    T* IReproj_dyyyy = new T [nxy];



   //fifth
    T* I1_dxyxxx = new T [nxy];
    T* I1_dxxxxx = new T [nxy];
    T* I1_dyxxxx = new T [nxy];
    T* I1_dyyxxx = new T [nxy];
    T* I1_dxyyxx = new T [nxy];
    T* I1_dxxyxx = new T [nxy];
    T* I1_dyxyxx = new T [nxy];
    T* I1_dyyyxx = new T [nxy];
    T* I1_dxyxyx = new T [nxy];
    T* I1_dxxxyx = new T [nxy];
    T* I1_dyxxyx = new T [nxy];
    T* I1_dyyxyx = new T [nxy];
    T* I1_dxyyyx = new T [nxy];
    T* I1_dxxyyx = new T [nxy];
    T* I1_dyxyyx = new T [nxy];
    T* I1_dyyyyx = new T [nxy];

    T* I1_dxyxxy = new T [nxy];
    T* I1_dxxxxy = new T [nxy];
    T* I1_dyxxxy = new T [nxy];
    T* I1_dyyxxy = new T [nxy];
    T* I1_dxyyxy = new T [nxy];
    T* I1_dxxyxy = new T [nxy];
    T* I1_dyxyxy = new T [nxy];
    T* I1_dyyyxy = new T [nxy];
    T* I1_dxyxyy = new T [nxy];
    T* I1_dxxxyy = new T [nxy];
    T* I1_dyxxyy = new T [nxy];
    T* I1_dyyxyy = new T [nxy];
    T* I1_dxyyyy = new T [nxy];
    T* I1_dxxyyy = new T [nxy];
    T* I1_dyxyyy = new T [nxy];
    T* I1_dyyyyy = new T [nxy];

    T* IReproj_dxyxxx = new T [nxy];
    T* IReproj_dxxxxx = new T [nxy];
    T* IReproj_dyxxxx = new T [nxy];
    T* IReproj_dyyxxx = new T [nxy];
    T* IReproj_dxyyxx = new T [nxy];
    T* IReproj_dxxyxx = new T [nxy];
    T* IReproj_dyxyxx = new T [nxy];
    T* IReproj_dyyyxx = new T [nxy];
    T* IReproj_dxyxyx = new T [nxy];
    T* IReproj_dxxxyx = new T [nxy];
    T* IReproj_dyxxyx = new T [nxy];
    T* IReproj_dyyxyx = new T [nxy];
    T* IReproj_dxyyyx = new T [nxy];
    T* IReproj_dxxyyx = new T [nxy];
    T* IReproj_dyxyyx = new T [nxy];
    T* IReproj_dyyyyx = new T [nxy];
    T* IReproj_dxyxxy = new T [nxy];
    T* IReproj_dxxxxy = new T [nxy];
    T* IReproj_dyxxxy = new T [nxy];
    T* IReproj_dyyxxy = new T [nxy];
    T* IReproj_dxyyxy = new T [nxy];
    T* IReproj_dxxyxy = new T [nxy];
    T* IReproj_dyxyxy = new T [nxy];
    T* IReproj_dyyyxy = new T [nxy];
    T* IReproj_dxyxyy = new T [nxy];
    T* IReproj_dxxxyy = new T [nxy];
    T* IReproj_dyxxyy = new T [nxy];
    T* IReproj_dyyxyy = new T [nxy];
    T* IReproj_dxyyyy = new T [nxy];
    T* IReproj_dxxyyy = new T [nxy];
    T* IReproj_dyxyyy = new T [nxy];
    T* IReproj_dyyyyy = new T [nxy];





    //first derivative
    if (depth>0){
      GaussianFirstDerivative2D(I1, I1_dx, sigma, 0, nx, ny);
      GaussianFirstDerivative2D(I1, I1_dy, sigma, 1, nx, ny);
      GaussianFirstDerivative2D(IReproj, IReproj_dx, sigma, 0, nx, ny);
      GaussianFirstDerivative2D(IReproj, IReproj_dy, sigma, 1, nx, ny);
      double CC1_X=crossCorrelationDistance(I1_dx, IReproj_dx, nxy, IMaskReproj);
      double CC1_Y=crossCorrelationDistance(I1_dy, IReproj_dy, nxy, IMaskReproj);
      CC1*=0.5*(CC1_X+CC1_Y);
    }

    //second derivative
    if (depth>1){
      GaussianFirstDerivative2D(I1_dx, I1_dxy, sigma, 1, nx, ny);
      GaussianFirstDerivative2D(I1_dx, I1_dxx, sigma, 0, nx, ny);
      GaussianFirstDerivative2D(I1_dy, I1_dyx, sigma, 1, nx, ny);
      GaussianFirstDerivative2D(I1_dy, I1_dyy, sigma, 0, nx, ny);
      GaussianFirstDerivative2D(IReproj_dx, IReproj_dxy, sigma, 1, nx, ny);
      GaussianFirstDerivative2D(IReproj_dx, IReproj_dxx, sigma, 0, nx, ny);
      GaussianFirstDerivative2D(IReproj_dy, IReproj_dyx, sigma, 1, nx, ny);
      GaussianFirstDerivative2D(IReproj_dy, IReproj_dyy, sigma, 0, nx, ny);
      
      double CC1_XY=crossCorrelationDistance(I1_dxy, IReproj_dxy, nxy, IMaskReproj);
      double CC1_XX=crossCorrelationDistance(I1_dxx, IReproj_dxx, nxy, IMaskReproj);
      double CC1_YX=crossCorrelationDistance(I1_dyx, IReproj_dyx, nxy, IMaskReproj);
      double CC1_YY=crossCorrelationDistance(I1_dyy, IReproj_dyy, nxy, IMaskReproj);
      CC1*= (1.0 / 4.0) * (CC1_XY + CC1_XX + CC1_YX + CC1_YY);
    }


    //third derivative
    if (depth>2){
    GaussianFirstDerivative2D(I1_dxy, I1_dxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxx, I1_dxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyx, I1_dyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyy, I1_dyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxy, I1_dxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxx, I1_dxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyx, I1_dyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyy, I1_dyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxy, IReproj_dxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxx, IReproj_dxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyx, IReproj_dyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyy, IReproj_dyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxy, IReproj_dxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxx, IReproj_dxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyx, IReproj_dyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyy, IReproj_dyyy, sigma, 1, nx, ny);
    double CC1_XYX=crossCorrelationDistance(I1_dxyx, IReproj_dxyx, nxy, IMaskReproj);
    double CC1_XXX=crossCorrelationDistance(I1_dxxx, IReproj_dxxx, nxy, IMaskReproj);
    double CC1_YXX=crossCorrelationDistance(I1_dyxx, IReproj_dyxx, nxy, IMaskReproj);
    double CC1_YYX=crossCorrelationDistance(I1_dyyx, IReproj_dyyx, nxy, IMaskReproj);
    double CC1_XYY=crossCorrelationDistance(I1_dxyy, IReproj_dxyy, nxy, IMaskReproj);
    double CC1_XXY=crossCorrelationDistance(I1_dxxy, IReproj_dxxy, nxy, IMaskReproj);
    double CC1_YXY=crossCorrelationDistance(I1_dyxy, IReproj_dyxy, nxy, IMaskReproj);
    double CC1_YYY=crossCorrelationDistance(I1_dyyy, IReproj_dyyy, nxy, IMaskReproj);
    CC1 += (1.0/8.0) * ( CC1_XYX + CC1_XXX + CC1_YXX + CC1_YYX + CC1_XYY + CC1_XXY + CC1_YXY + CC1_YYY);
    }



    //forth derivative
    if (depth>3){
    GaussianFirstDerivative2D(I1_dxyx, I1_dxyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxxx, I1_dxxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyxx, I1_dyxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyyx, I1_dyyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxyy, I1_dxyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxxy, I1_dxxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyxy, I1_dyxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyyy, I1_dyyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxyx, I1_dxyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxxx, I1_dxxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyxx, I1_dyxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyyx, I1_dyyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxyy, I1_dxyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxxy, I1_dxxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyxy, I1_dyxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyyy, I1_dyyyy, sigma, 1, nx, ny);

    GaussianFirstDerivative2D(IReproj_dxyx, IReproj_dxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxx, IReproj_dxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxx, IReproj_dyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyx, IReproj_dyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyy, IReproj_dxyy, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxy, IReproj_dxxy, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxy, IReproj_dyxy, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyy, IReproj_dyyy, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyx, IReproj_dxyx, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxx, IReproj_dxxx, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxx, IReproj_dyxx, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyx, IReproj_dyyx, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyy, IReproj_dxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxy, IReproj_dxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxy, IReproj_dyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyy, IReproj_dyyy, sigma, 1, nx, ny);


    double CC1_XYXX=crossCorrelationDistance(I1_dxyxx, IReproj_dxyxx, nxy, IMaskReproj);
    double CC1_XXXX=crossCorrelationDistance(I1_dxxxx, IReproj_dxxxx, nxy, IMaskReproj);
    double CC1_YXXX=crossCorrelationDistance(I1_dyxxx, IReproj_dyxxx, nxy, IMaskReproj);
    double CC1_YYXX=crossCorrelationDistance(I1_dyyxx, IReproj_dyyxx, nxy, IMaskReproj);
    double CC1_XYYX=crossCorrelationDistance(I1_dxyyx, IReproj_dxyyx, nxy, IMaskReproj);
    double CC1_XXYX=crossCorrelationDistance(I1_dxxyx, IReproj_dxxyx, nxy, IMaskReproj);
    double CC1_YXYX=crossCorrelationDistance(I1_dyxyx, IReproj_dyxyx, nxy, IMaskReproj);
    double CC1_YYYX=crossCorrelationDistance(I1_dyyyx, IReproj_dyyyx, nxy, IMaskReproj);
    double CC1_XYXY=crossCorrelationDistance(I1_dxyxy, IReproj_dxyxy, nxy, IMaskReproj);
    double CC1_XXXY=crossCorrelationDistance(I1_dxxxy, IReproj_dxxxy, nxy, IMaskReproj);
    double CC1_YXXY=crossCorrelationDistance(I1_dyxxy, IReproj_dyxxy, nxy, IMaskReproj);
    double CC1_YYXY=crossCorrelationDistance(I1_dyyxy, IReproj_dyyxy, nxy, IMaskReproj);
    double CC1_XYYY=crossCorrelationDistance(I1_dxyyy, IReproj_dxyyy, nxy, IMaskReproj);
    double CC1_XXYY=crossCorrelationDistance(I1_dxxyy, IReproj_dxxyy, nxy, IMaskReproj);
    double CC1_YXYY=crossCorrelationDistance(I1_dyxyy, IReproj_dyxyy, nxy, IMaskReproj);
    double CC1_YYYY=crossCorrelationDistance(I1_dyyyy, IReproj_dyyyy, nxy, IMaskReproj);


    CC1 += (1.0/16.0) * ( CC1_XYXX + CC1_XXXX + CC1_YXXX + CC1_YYXX + CC1_XYYX + CC1_XXYX + CC1_YXYX + CC1_YYYX
                          +CC1_XYXY + CC1_XXXY + CC1_YXXY + CC1_YYXY + CC1_XYYY + CC1_XXYY + CC1_YXYY + CC1_YYYY);
    }



    //FIFTH derivative
    if (depth>4){
    GaussianFirstDerivative2D(I1_dxyxx, I1_dxyxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxxxx, I1_dxxxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyxxx, I1_dyxxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyyxx, I1_dyyxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxyyx, I1_dxyyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxxyx, I1_dxxyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyxyx, I1_dyxyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyyyx, I1_dyyyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxyxy, I1_dxyxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxxxy, I1_dxxxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyxxy, I1_dyxxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyyxy, I1_dyyxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxyyy, I1_dxyyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxxyy, I1_dxxyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyxyy, I1_dyxyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dyyyy, I1_dyyyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(I1_dxyxx, I1_dxyxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxxxx, I1_dxxxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyxxx, I1_dyxxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyyxx, I1_dyyxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxyyx, I1_dxyyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxxyx, I1_dxxyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyxyx, I1_dyxyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyyyx, I1_dyyyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxyxy, I1_dxyxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxxxy, I1_dxxxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyxxy, I1_dyxxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyyxy, I1_dyyxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxyyy, I1_dxyyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dxxyy, I1_dxxyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyxyy, I1_dyxyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(I1_dyyyy, I1_dyyyyy, sigma, 1, nx, ny);



    GaussianFirstDerivative2D(IReproj_dxyx, IReproj_dxyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxx, IReproj_dxxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxx, IReproj_dyxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyx, IReproj_dyyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyy, IReproj_dxyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxy, IReproj_dxxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxy, IReproj_dyxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyy, IReproj_dyyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyx, IReproj_dxyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxx, IReproj_dxxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxx, IReproj_dyxxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyx, IReproj_dyyxx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyy, IReproj_dxyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxy, IReproj_dxxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxy, IReproj_dyxyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyy, IReproj_dyyyx, sigma, 0, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyx, IReproj_dxyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxx, IReproj_dxxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxx, IReproj_dyxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyx, IReproj_dyyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyy, IReproj_dxyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxy, IReproj_dxxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxy, IReproj_dyxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyy, IReproj_dyyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyx, IReproj_dxyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxx, IReproj_dxxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxx, IReproj_dyxxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyx, IReproj_dyyxy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxyy, IReproj_dxyyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dxxy, IReproj_dxxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyxy, IReproj_dyxyy, sigma, 1, nx, ny);
    GaussianFirstDerivative2D(IReproj_dyyy, IReproj_dyyyy, sigma, 1, nx, ny);


    double CC1_XYXXX=crossCorrelationDistance(I1_dxyxxx, IReproj_dxyxxx, nxy, IMaskReproj);
    double CC1_XXXXX=crossCorrelationDistance(I1_dxxxxx, IReproj_dxxxxx, nxy, IMaskReproj);
    double CC1_YXXXX=crossCorrelationDistance(I1_dyxxxx, IReproj_dyxxxx, nxy, IMaskReproj);
    double CC1_YYXXX=crossCorrelationDistance(I1_dyyxxx, IReproj_dyyxxx, nxy, IMaskReproj);
    double CC1_XYYXX=crossCorrelationDistance(I1_dxyyxx, IReproj_dxyyxx, nxy, IMaskReproj);
    double CC1_XXYXX=crossCorrelationDistance(I1_dxxyxx, IReproj_dxxyxx, nxy, IMaskReproj);
    double CC1_YXYXX=crossCorrelationDistance(I1_dyxyxx, IReproj_dyxyxx, nxy, IMaskReproj);
    double CC1_YYYXX=crossCorrelationDistance(I1_dyyyxx, IReproj_dyyyxx, nxy, IMaskReproj);
    double CC1_XYXYX=crossCorrelationDistance(I1_dxyxyx, IReproj_dxyxyx, nxy, IMaskReproj);
    double CC1_XXXYX=crossCorrelationDistance(I1_dxxxyx, IReproj_dxxxyx, nxy, IMaskReproj);
    double CC1_YXXYX=crossCorrelationDistance(I1_dyxxyx, IReproj_dyxxyx, nxy, IMaskReproj);
    double CC1_YYXYX=crossCorrelationDistance(I1_dyyxyx, IReproj_dyyxyx, nxy, IMaskReproj);
    double CC1_XYYYX=crossCorrelationDistance(I1_dxyyyx, IReproj_dxyyyx, nxy, IMaskReproj);
    double CC1_XXYYX=crossCorrelationDistance(I1_dxxyyx, IReproj_dxxyyx, nxy, IMaskReproj);
    double CC1_YXYYX=crossCorrelationDistance(I1_dyxyyx, IReproj_dyxyyx, nxy, IMaskReproj);
    double CC1_YYYYX=crossCorrelationDistance(I1_dyyyyx, IReproj_dyyyyx, nxy, IMaskReproj);


    double CC1_XYXXY=crossCorrelationDistance(I1_dxyxxy, IReproj_dxyxxy, nxy, IMaskReproj);
    double CC1_XXXXY=crossCorrelationDistance(I1_dxxxxy, IReproj_dxxxxy, nxy, IMaskReproj);
    double CC1_YXXXY=crossCorrelationDistance(I1_dyxxxy, IReproj_dyxxxy, nxy, IMaskReproj);
    double CC1_YYXXY=crossCorrelationDistance(I1_dyyxxy, IReproj_dyyxxy, nxy, IMaskReproj);
    double CC1_XYYXY=crossCorrelationDistance(I1_dxyyxy, IReproj_dxyyxy, nxy, IMaskReproj);
    double CC1_XXYXY=crossCorrelationDistance(I1_dxxyxy, IReproj_dxxyxy, nxy, IMaskReproj);
    double CC1_YXYXY=crossCorrelationDistance(I1_dyxyxy, IReproj_dyxyxy, nxy, IMaskReproj);
    double CC1_YYYXY=crossCorrelationDistance(I1_dyyyxy, IReproj_dyyyxy, nxy, IMaskReproj);
    double CC1_XYXYY=crossCorrelationDistance(I1_dxyxyy, IReproj_dxyxyy, nxy, IMaskReproj);
    double CC1_XXXYY=crossCorrelationDistance(I1_dxxxyy, IReproj_dxxxyy, nxy, IMaskReproj);
    double CC1_YXXYY=crossCorrelationDistance(I1_dyxxyy, IReproj_dyxxyy, nxy, IMaskReproj);
    double CC1_YYXYY=crossCorrelationDistance(I1_dyyxyy, IReproj_dyyxyy, nxy, IMaskReproj);
    double CC1_XYYYY=crossCorrelationDistance(I1_dxyyyy, IReproj_dxyyyy, nxy, IMaskReproj);
    double CC1_XXYYY=crossCorrelationDistance(I1_dxxyyy, IReproj_dxxyyy, nxy, IMaskReproj);
    double CC1_YXYYY=crossCorrelationDistance(I1_dyxyyy, IReproj_dyxyyy, nxy, IMaskReproj);
    double CC1_YYYYY=crossCorrelationDistance(I1_dyyyyy, IReproj_dyyyyy, nxy, IMaskReproj);

    CC1 += (1.0/32.0) * ( CC1_XYXXX + CC1_XXXXX + CC1_YXXXX + CC1_YYXXX + CC1_XYYXX + CC1_XXYXX + CC1_YXYXX + CC1_YYYXX
                          +CC1_XYXYX + CC1_XXXYX + CC1_YXXYX + CC1_YYXYX + CC1_XYYYX + CC1_XXYYX + CC1_YXYYX + CC1_YYYYX
                          +CC1_XYXXY + CC1_XXXXY + CC1_YXXXY + CC1_YYXXY + CC1_XYYXY + CC1_XXYXY + CC1_YXYXY + CC1_YYYXY
                          +CC1_XYXYY + CC1_XXXYY + CC1_YXXYY + CC1_YYXYY + CC1_XYYYY + CC1_XXYYY + CC1_YXYYY + CC1_YYYYY);
    }





    double returnVal = 0;
    if (!(CC1!=CC1)){
        returnVal =  CC1;
        if (CC1>_MIN_VAL){
           returnVal =  CC1;
        }else{
            returnVal =  _MIN_VAL;
        }
    }


  //first derivative
  delete [] I1_dx;
  delete [] I1_dy;
  delete [] IReproj_dx;
  delete [] IReproj_dy;

  //second derivative
  delete [] I1_dxy;
  delete [] I1_dxx;
  delete [] I1_dyx;
  delete [] I1_dyy;

  delete [] IReproj_dxy;
  delete [] IReproj_dxx;
  delete [] IReproj_dyx;
  delete [] IReproj_dyy;

  //third derivative
  delete [] I1_dxyx;
  delete [] I1_dxxx;
  delete [] I1_dyxx;
  delete [] I1_dyyx;
  delete [] I1_dxyy;
  delete [] I1_dxxy;
  delete [] I1_dyxy;
  delete [] I1_dyyy;


  delete [] IReproj_dxyx;
  delete [] IReproj_dxxx;
  delete [] IReproj_dyxx;
  delete [] IReproj_dyyx;
  delete [] IReproj_dxyy;
  delete [] IReproj_dxxy;
  delete [] IReproj_dyxy;
  delete [] IReproj_dyyy;


  //forth derivative
  delete [] I1_dxyxx;
  delete [] I1_dxxxx;
  delete [] I1_dyxxx;
  delete [] I1_dyyxx;
  delete [] I1_dxyyx;
  delete [] I1_dxxyx;
  delete [] I1_dyxyx;
  delete [] I1_dyyyx;
  delete [] I1_dxyxy;
  delete [] I1_dxxxy;
  delete [] I1_dyxxy;
  delete [] I1_dyyxy;
  delete [] I1_dxyyy;
  delete [] I1_dxxyy;
  delete [] I1_dyxyy;
  delete [] I1_dyyyy;


  delete [] IReproj_dxyxx;
  delete [] IReproj_dxxxx;
  delete [] IReproj_dyxxx;
  delete [] IReproj_dyyxx;
  delete [] IReproj_dxyyx;
  delete [] IReproj_dxxyx;
  delete [] IReproj_dyxyx;
  delete [] IReproj_dyyyx;
  delete [] IReproj_dxyxy;
  delete [] IReproj_dxxxy;
  delete [] IReproj_dyxxy;
  delete [] IReproj_dyyxy;
  delete [] IReproj_dxyyy;
  delete [] IReproj_dxxyy;
  delete [] IReproj_dyxyy;
  delete [] IReproj_dyyyy;



  //fifth derivative
  delete [] I1_dxyxxx;
  delete [] I1_dxxxxx;
  delete [] I1_dyxxxx;
  delete [] I1_dyyxxx;
  delete [] I1_dxyyxx;
  delete [] I1_dxxyxx;
  delete [] I1_dyxyxx;
  delete [] I1_dyyyxx;
  delete [] I1_dxyxyx;
  delete [] I1_dxxxyx;
  delete [] I1_dyxxyx;
  delete [] I1_dyyxyx;
  delete [] I1_dxyyyx;
  delete [] I1_dxxyyx;
  delete [] I1_dyxyyx;
  delete [] I1_dyyyyx;
  delete [] I1_dxyxxy;
  delete [] I1_dxxxxy;
  delete [] I1_dyxxxy;
  delete [] I1_dyyxxy;
  delete [] I1_dxyyxy;
  delete [] I1_dxxyxy;
  delete [] I1_dyxyxy;
  delete [] I1_dyyyxy;
  delete [] I1_dxyxyy;
  delete [] I1_dxxxyy;
  delete [] I1_dyxxyy;
  delete [] I1_dyyxyy;
  delete [] I1_dxyyyy;
  delete [] I1_dxxyyy;
  delete [] I1_dyxyyy;
  delete [] I1_dyyyyy;


  delete [] IReproj_dxyxxx;
  delete [] IReproj_dxxxxx;
  delete [] IReproj_dyxxxx;
  delete [] IReproj_dyyxxx;
  delete [] IReproj_dxyyxx;
  delete [] IReproj_dxxyxx;
  delete [] IReproj_dyxyxx;
  delete [] IReproj_dyyyxx;
  delete [] IReproj_dxyxyx;
  delete [] IReproj_dxxxyx;
  delete [] IReproj_dyxxyx;
  delete [] IReproj_dyyxyx;
  delete [] IReproj_dxyyyx;
  delete [] IReproj_dxxyyx;
  delete [] IReproj_dyxyyx;
  delete [] IReproj_dyyyyx;
  delete [] IReproj_dxyxxy;
  delete [] IReproj_dxxxxy;
  delete [] IReproj_dyxxxy;
  delete [] IReproj_dyyxxy;
  delete [] IReproj_dxyyxy;
  delete [] IReproj_dxxyxy;
  delete [] IReproj_dyxyxy;
  delete [] IReproj_dyyyxy;
  delete [] IReproj_dxyxyy;
  delete [] IReproj_dxxxyy;
  delete [] IReproj_dyxxyy;
  delete [] IReproj_dyyxyy;
  delete [] IReproj_dxyyyy;
  delete [] IReproj_dxxyyy;
  delete [] IReproj_dyxyyy;
  delete [] IReproj_dyyyyy;


  return returnVal;
}


// **********************************
//  largest eigen value of symmetric matrix in the form
//    | a  b |
//    | b  c |
//
double largestEigenvalue(double a, double b, double c){
    double delta = pow(a-c,2.0)+4.0*b*b;
    if (delta<0) return 0;
    delta=pow(delta,0.5);
//    return ((a+c+delta)/2.0);
    return (a+c+delta);
}

//
template<typename T>
double computeEigenScore (T* I1, T* IReproj, T* IMaskReproj, unsigned long int nx, unsigned long int ny, int depth, double sigma){
    const double _MIN_VAL=-1.0;//0.0f;//0.000001;
    const unsigned long int nxy = nx * ny;
    double CC1=0;//crossCorrelationDistance(I1, IReproj, nxy, IMaskReproj);

    double angpix = 1;
    int padding = 3;

    float * Ixx = new float[nxy];
    float * Iyy = new float[nxy];
    float * Ixy = new float[nxy];
    float * IReproj_xx= new float[nxy];
    float * IReproj_yy= new float[nxy];
    float * IReproj_xy= new float[nxy];
    float * I1_Eigen = new float[nxy];
    float * IReproj_Eigen = new float[nxy];


    if (depth>1){
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, I1, Ixx);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 1, depth, I1, Iyy);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth-1, I1, Ixy);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 2, depth-1, Ixy);

      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, IReproj, IReproj_xx);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 1, depth, IReproj, IReproj_yy);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth-1, IReproj, IReproj_xy);
      gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 2, depth-1, IReproj_xy);

      for (unsigned long int ii=0; ii<nxy; ii++){
          double tmpEigen_I1=largestEigenvalue(Ixx[ii],Ixy[ii],Iyy[ii]);
          double tmpEigen_Reproj=largestEigenvalue(IReproj_xx[ii],IReproj_xy[ii],IReproj_yy[ii]);
          I1_Eigen[ii]=tmpEigen_I1;
          IReproj_Eigen[ii]=tmpEigen_Reproj;
      }
      double CC_eigen=crossCorrelationDistance(I1_Eigen, IReproj_Eigen, nxy, IMaskReproj);
      if (CC_eigen<0) CC_eigen = 0;
      CC1=CC_eigen;
    }

    double returnVal = 0;
    if (!(CC1!=CC1)){
        returnVal =  CC1;
        if (CC1>_MIN_VAL){
           returnVal =  CC1;
        }else{
            returnVal = _MIN_VAL;
        }
    }


    delete [] Ixx;
    delete [] Iyy;
    delete [] Ixy;

    delete [] IReproj_xx;
    delete [] IReproj_yy;
    delete [] IReproj_xy;
    delete [] I1_Eigen;
    delete [] IReproj_Eigen;

  return returnVal;

}



template<typename T>
double computeSingleScores (T* I1, T* IReproj, T* IMaskReproj, unsigned long int nx, unsigned long int ny, int depth, double sigma){
  //return computeDeepDerivativeScore ( I1, IReproj, IMaskReproj, nx, ny, depth, sigma);
  return computeDerivatives ( I1, IReproj, IMaskReproj, nx, ny, depth, sigma);
  //return computeEigenScore (I1, IReproj, IMaskReproj, nx, ny, depth, sigma);
}

template<typename T>
std::vector<double> computeMultipleScores (T* I1, T* IReproj, T* IMaskReproj, unsigned long int nx, unsigned long int ny, double sigma){
  //const double _MIN_VAL=0.0f;//0.000001;
  const unsigned long int nxy = nx * ny;
  double angpix = 1;
  int padding = 3;
  std::vector<double> result;
  const double minVal = -1;
  double tmpVal = crossCorrelationDistance(I1, IReproj, nxy, IMaskReproj);
  if (tmpVal!=tmpVal) tmpVal =  minVal;
  result.push_back(tmpVal);

  //derivatives
  float * Ix = new float[nxy];
  float * Iy = new float[nxy];
  float * IReproj_x= new float[nxy];
  float * IReproj_y= new float[nxy];
  int depth = 1;
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, I1, Ix);
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, IReproj, IReproj_x);
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, I1, Iy);
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, IReproj, IReproj_y);
  tmpVal=crossCorrelationDistance(Ix, IReproj_x, nxy, IMaskReproj);
  if (tmpVal!=tmpVal) tmpVal =  minVal;
  result.push_back(tmpVal);
  tmpVal=crossCorrelationDistance(Iy, IReproj_y, nxy, IMaskReproj);
  if (tmpVal!=tmpVal) tmpVal =  minVal;
  result.push_back(tmpVal);

  depth = 2;
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, I1, Ix);
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, IReproj, IReproj_x);
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, I1, Iy);
  gaussRecursiveDerivatives1D (sigma,  nx,  ny,  1, angpix, padding, 0, depth, IReproj, IReproj_y);
  tmpVal=crossCorrelationDistance(Ix, IReproj_x, nxy, IMaskReproj);
  if (tmpVal!=tmpVal) tmpVal =  minVal;
  result.push_back(tmpVal);
  tmpVal=crossCorrelationDistance(Iy, IReproj_y, nxy, IMaskReproj);
  if (tmpVal!=tmpVal) tmpVal =  minVal;
  result.push_back(tmpVal);
  tmpVal=computeEigenScore (I1, IReproj, IMaskReproj, nx, ny, depth, sigma);
  if (tmpVal!=tmpVal) tmpVal =  minVal;
  result.push_back(tmpVal);

  delete [] Ix;
  delete [] Iy;
  delete [] IReproj_x;
  delete [] IReproj_y;

  return result;
}





#endif
