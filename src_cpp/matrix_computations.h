#ifndef ___MATRIX_COMPUTE___H___
#define ___MATRIX_COMPUTE___H___
#include <vector>
#include "mrcIO.h"
// ****************************
// matrixMultiplication (useful for computing the transform matrix)
//   space of matrixes already allocated
//   size of outM=x2*y1
//
void matrixMultiplication(double M1[], double M2[], double outM[],  unsigned int x1,  unsigned int y1, unsigned int x2,  unsigned int y2, bool verbose=false)
{
    unsigned int xOut = x2;
    unsigned int yOut = y1;
    unsigned int c, d, i;
    double sum =0;

   if(verbose){
        printf("M1 is:\n");
        for ( c = 0 ; c < y1 ; c++ ){
            for ( d = 0 ; d < x1 ; d++ )
                printf("%f\t", M1[d+x1*c]);
            printf("\n");
        }
        printf("M2 is:\n");
        for ( c = 0 ; c < y2 ; c++ ){
            for ( d = 0 ; d < x2 ; d++ )
                printf("%f\t", M2[d+x2*c]);
            printf("\n");
        }
    }


   if (x1!=y2){
     if(verbose)
      std::cerr<<"not possible to multiply\n";
     return;
   }
   for (i=0;i<y1*x2;i++){
     outM[i]=0;
   }    

    // matrix multiplication operation
    for ( unsigned int jj = 0 ; jj < y1 ; jj++ ){
           unsigned int jjx1=jj*x1;
           unsigned int jjx2=jj*x2;
           for ( unsigned int ii = 0 ; ii< x2 ; ii++ ){
                sum = 0;
                for ( unsigned int kk = 0 ; kk < x1 ; kk++ ){
                    sum = sum + M1[kk+jjx1]*M2[ii+kk*x2];
                    if(verbose) std::cerr<<M1[kk+jjx1]<<"*"<<M2[ii+kk*x2]<<", ";
                }
                if(verbose) std::cerr<<"="<<sum<<"\n";
                outM[ii+jjx2] = sum;
                
            }
    }


        //Printing the final product matrix
    if(verbose){
        printf("\nThe product of entered matrices is:\n");
        for (unsigned int jj = 0 ; jj < yOut ; jj++ ){
            for ( unsigned int ii = 0 ; ii < xOut ; ii++ ) {
                printf("%f\t", outM[ii+jj*xOut]);
            }
            printf("\n");
        }
    }
 
}


// ****************************
// transformMatrix
// returns a 3x3 transform matrix, rigid 
void transformMatrixRelion(double M[], double Phi, double Theta, double Psi, double tx=0, double ty=0, double tz=0){
 
    double ca, sa, cb, sb, cg, sg;
    double cc, cs, sc, ss;
    //Phi=fmod(Phi,360.0);
    //Theta=fmod(Theta,180.0);
    //Psi=fmod(Psi,360.0);
    double pi_180=M_PI/double(180.0);
    double alpha = Phi*pi_180;
    double beta  = Theta*pi_180;
    double gamma = Psi*pi_180;

    ca = cos(alpha);
    cb = cos(beta);
    cg = cos(gamma);
    sa = sin(alpha);
    sb = sin(beta);
    sg = sin(gamma);
    cc = cb * ca;
    cs = cb * sa;
    sc = sb * ca;
    ss = sb * sa;

  M[0]=cg * cc - sg * sa;    M[4]=cg * cs + sg * ca;    M[8]=-cg * sb;    M[3]=tx;
  M[1]=-sg * cc - cg * sa;   M[5]=-sg * cs + cg * ca;   M[9]=sg * sb;     M[7]=ty;
  M[2]=sc;                   M[6]=ss;                   M[10]=cb;         M[11]=tz;
  M[12]=0;                   M[13]=0;                   M[14]=0;          M[15]=1;

}





//    [a11 a12 a13 a14]
//    [a21 a22 a23 a24]
// M= [a31 a32 a33 a34]
//    [a41 a42 a43 a44]
void inverseMatrix(double M[]){
/*
    double A11=M[5]*M[10]*M[15] + M[9]*M[14]*M[7] + M[13]*M[6]*M[11]
                -M[1+4*3]*M[2+4*2]*M[3+4*1] - M[1+4*2]*M[2+4*1]*M[3+4*3] - M[1+4*1]*M[2+4*3]*M[3+4*2];
    double A12= -M[0+4*1]*M[2+4*2]*M[3+4*3] - M[0+4*2]*M[2+4*3]*M[3+4*1] - M[0+4*3]*M[2+4*1]*M[3+4*2]
                +M[0+4*3]*M[2+4*2]*M[3+4*1] + M[0+4*2]*M[2+4*1]*M[3+4*3] + M[0+4*1]*M[2+4*3]*M[3+4*2];
    double A13=M[0+4*1]*M[1+4*2]*M[3+4*3] + M[0+4*2]*M[1+4*3]*M[3+4*1] + M[0+4*3]*M[1+4*1]*M[3+4*2]
                -M[0+4*3]*M[1+4*2]*M[3+4*1] - M[0+4*2]*M[1+4*1]*M[3+4*3] - M[0+4*1]*M[1+4*3]*M[3+4*2];
    double A14= -M[0+4*1]*M[1+4*2]*M[2+4*3] - M[0+4*2]*M[1+4*3]*M[2+4*1] - M[0+4*3]*M[1+4*1]*M[2+4*2]
                +M[0+4*3]*M[1+4*2]*M[2+4*1] + M[0+4*2]*M[1+4*1]*M[2+4*3] + M[0+4*1]*M[1+4*3]*M[2+4*2];

    double A21= -M[1+4*0]*M[2+4*2]*M[3+4*3] - M[1+4*2]*M[2+4*3]*M[3+4*0] - M[1+4*3]*M[2+4*0]*M[3+4*2]
                +M[1+4*3]*M[2+4*2]*M[3+4*0] + M[1+4*2]*M[2+4*0]*M[3+4*3] + M[1+4*0]*M[2+4*3]*M[3+4*2];
    double A22=M[0+4*0]*M[2+4*2]*M[3+4*3] + M[0+4*2]*M[2+4*3]*M[3+4*0] + M[0+4*3]*M[2+4*0]*M[3+4*2]
                -M[0+4*3]*M[2+4*2]*M[3+4*0] - M[0+4*2]*M[2+4*0]*M[3+4*3] - M[0+4*0]*M[2+4*3]*M[3+4*2];
    double A23= -M[0+4*0]*M[1+4*2]*M[3+4*3] - M[0+4*2]*M[1+4*3]*M[3+4*0] - M[0+4*3]*M[1+4*0]*M[3+4*2]
                +M[0+4*3]*M[1+4*2]*M[3+4*0] + M[0+4*2]*M[1+4*0]*M[3+4*3] + M[0+4*0]*M[1+4*3]*M[3+4*2];
    double A24=M[0+4*0]*M[1+4*2]*M[2+4*3] + M[0+4*2]*M[1+4*3]*M[2+4*0] + M[0+4*3]*M[1+4*0]*M[2+4*2]
                -M[0+4*3]*M[1+4*2]*M[2+4*0] - M[0+4*2]*M[1+4*0]*M[2+4*3] - M[0+4*0]*M[1+4*3]*M[2+4*2];


    double A31=M[1+4*0]*M[2+4*1]*M[3+4*3] + M[1+4*1]*M[2+4*3]*M[3+4*0] + M[1+4*3]*M[2+4*0]*M[3+4*1]
                -M[1+4*3]*M[2+4*1]*M[3+4*0] - M[1+4*1]*M[2+4*0]*M[3+4*3] - M[1+4*0]*M[2+4*3]*M[3+4*1];
    double A32= -M[0+4*0]*M[2+4*1]*M[3+4*3] - M[0+4*1]*M[2+4*3]*M[3+4*0] - M[0+4*3]*M[2+4*0]*M[3+4*1]
                +M[0+4*3]*M[2+4*1]*M[3+4*0] + M[0+4*1]*M[2+4*0]*M[3+4*3] + M[0+4*0]*M[2+4*3]*M[3+4*1];
    double A33=M[0+4*0]*M[1+4*1]*M[3+4*3] + M[0+4*1]*M[1+4*3]*M[3+4*0] + M[0+4*3]*M[1+4*0]*M[3+4*1]
                -M[0+4*3]*M[1+4*1]*M[3+4*0] - M[0+4*1]*M[1+4*0]*M[3+4*3] - M[0+4*0]*M[1+4*3]*M[3+4*1];
    double A34= -M[0+4*0]*M[1+4*1]*M[2+4*3] - M[0+4*1]*M[1+4*3]*M[2+4*0] - M[0+4*3]*M[1+4*0]*M[2+4*1]
                +M[0+4*3]*M[1+4*1]*M[2+4*0] + M[0+4*1]*M[1+4*0]*M[2+4*3] + M[0+4*0]*M[1+4*3]*M[2+4*1];


    double A41= -M[1+4*0]*M[2+4*1]*M[3+4*2] - M[1+4*1]*M[2+4*2]*M[3+4*0] - M[1+4*2]*M[2+4*0]*M[3+4*1]
                +M[1+4*2]*M[2+4*1]*M[3+4*0] + M[1+4*1]*M[2+4*0]*M[3+4*2] + M[1+4*0]*M[2+4*2]*M[3+4*1];
    double A42=M[0+4*0]*M[2+4*1]*M[3+4*2] + M[0+4*1]*M[2+4*2]*M[3+4*0] + M[0+4*2]*M[2+4*0]*M[3+4*1]
                -M[0+4*2]*M[2+4*1]*M[3+4*0] - M[0+4*1]*M[2+4*0]*M[3+4*2] - M[0+4*0]*M[2+4*2]*M[3+4*1];
    double A43= -M[0+4*0]*M[1+4*1]*M[3+4*2] - M[0+4*1]*M[1+4*2]*M[3+4*0] - M[0+4*2]*M[1+4*0]*M[3+4*1]
                +M[0+4*2]*M[1+4*1]*M[3+4*0] + M[0+4*1]*M[1+4*0]*M[3+4*2] + M[0+4*0]*M[1+4*2]*M[3+4*1];
    double A44=M[0+4*0]*M[1+4*1]*M[2+4*2] + M[0+4*1]*M[1+4*2]*M[2+4*0] + M[0+4*2]*M[1+4*0]*M[2+4*1]
                -M[0+4*2]*M[1+4*1]*M[2+4*0] - M[0+4*1]*M[1+4*0]*M[2+4*2] - M[0+4*0]*M[1+4*2]*M[2+4*1];

*/
    double A11=M[5]*M[10]*M[15] + M[9]*M[14]*M[7] + M[13]*M[6]*M[11]
                -M[13]*M[10]*M[7] - M[9]*M[6]*M[15] - M[5]*M[14]*M[11];
    double A12= -M[4]*M[10]*M[15] - M[8]*M[14]*M[7] - M[12]*M[6]*M[11]
                +M[12]*M[10]*M[7] + M[8]*M[6]*M[15] + M[4]*M[14]*M[11];
    double A13=M[4]*M[9]*M[15] + M[8]*M[13]*M[7] + M[12]*M[5]*M[11]
                -M[12]*M[9]*M[7] - M[8]*M[5]*M[15] - M[4]*M[13]*M[11];
    double A14= -M[4]*M[9]*M[14] - M[8]*M[13]*M[6] - M[12]*M[5]*M[10]
                +M[12]*M[9]*M[6] + M[8]*M[5]*M[14] + M[4]*M[13]*M[10];

    double A21= -M[1]*M[10]*M[15] - M[9]*M[14]*M[3] - M[13]*M[2]*M[11]
                +M[13]*M[10]*M[3] + M[9]*M[2]*M[15] + M[1]*M[14]*M[11];
    double A22=M[0]*M[10]*M[15] + M[8]*M[14]*M[3] + M[12]*M[2]*M[11]
                -M[12]*M[10]*M[3] - M[8]*M[2]*M[15] - M[0]*M[14]*M[11];
    double A23= -M[0]*M[9]*M[15] - M[8]*M[13]*M[3] - M[12]*M[1]*M[11]
                +M[12]*M[9]*M[3] + M[8]*M[1]*M[15] + M[0]*M[13]*M[11];
    double A24=M[0]*M[9]*M[14] + M[8]*M[13]*M[2] + M[12]*M[1]*M[10]
                -M[12]*M[9]*M[2] - M[8]*M[1]*M[14] - M[0]*M[13]*M[10];


    double A31=M[1]*M[6]*M[15] + M[5]*M[14]*M[3] + M[13]*M[2]*M[7]
                -M[13]*M[6]*M[3] - M[5]*M[2]*M[15] - M[1]*M[14]*M[7];
    double A32= -M[0]*M[6]*M[15] - M[4]*M[14]*M[3] - M[12]*M[2]*M[7]
                +M[12]*M[6]*M[3] + M[4]*M[2]*M[15] + M[0]*M[14]*M[7];
    double A33=M[0]*M[5]*M[15] + M[4]*M[13]*M[3] + M[12]*M[1]*M[7]
                -M[12]*M[5]*M[3] - M[4]*M[1]*M[15] - M[0]*M[13]*M[7];
    double A34= -M[0]*M[5]*M[14] - M[4]*M[13]*M[2] - M[12]*M[1]*M[6]
                +M[12]*M[5]*M[2] + M[4]*M[1]*M[14] + M[0]*M[13]*M[6];


    double A41= -M[1]*M[6]*M[11] - M[5]*M[10]*M[3] - M[9]*M[2]*M[7]
                +M[9]*M[6]*M[3] + M[5]*M[2]*M[11] + M[1]*M[10]*M[7];
    double A42=M[0]*M[6]*M[11] + M[4]*M[10]*M[3] + M[8]*M[2]*M[7]
                -M[8]*M[6]*M[3] - M[4]*M[2]*M[11] - M[0]*M[10]*M[7];
    double A43= -M[0]*M[5]*M[11] - M[4]*M[9]*M[3] - M[8]*M[1]*M[7]
                +M[8]*M[5]*M[3] + M[4]*M[1]*M[11] + M[0]*M[9]*M[7];
    double A44=M[0]*M[5]*M[10] + M[4]*M[9]*M[2] + M[8]*M[1]*M[6]
                -M[8]*M[5]*M[2] - M[4]*M[1]*M[10] - M[0]*M[9]*M[6];

    M[0]=A11;M[1]=A12;M[2]=A13;M[3]=A14;
    M[4]=A21;M[5]=A22;M[6]=A23;M[7]=A24;
    M[8]=A31;M[9]=A32;M[10]=A33;M[11]=A34;
    M[12]=A41;M[13]=A42;M[14]=A43;M[15]=A44;
}






//SORT
template<typename T>
std::vector<T> bubbleSortAscendingValues(const std::vector<T> valuesIn){
  std::vector<T> valuesOut = valuesIn;
  long int size = valuesOut.size();
  for (long int i = (size - 1); i > 0; i--)
    {
      for (long int j = 1; j <= i; j++)
	{
	  if (valuesOut[j - 1] > valuesOut[j])
	    {
	      T temp = valuesOut[j - 1];
	      valuesOut[j - 1] = valuesOut[j];
	      valuesOut[j] = temp;
	    }
	}
    }
  return valuesOut;
}


void invertRigidTransformInPlace(double M[16]) {
    // Temporary buffers for rotation and translation
    double R[3][3] = {
        {M[0], M[1], M[2]},
        {M[4], M[5], M[6]},
        {M[8], M[9], M[10]}
    };
    double T[3] = {M[12], M[13], M[14]};
    double tempT[3];

    // Compute the transposed rotation matrix (in-place)
    for (int i = 0; i < 3; i++) {
        for (int j = i+1; j < 3; j++) {
            double temp = R[i][j];
            R[i][j] = R[j][i];
            R[j][i] = temp;
        }
    }

    // Compute the translated vector
    for (int i = 0; i < 3; i++) {
        tempT[i] = -(R[0][i] * T[0] + R[1][i] * T[1] + R[2][i] * T[2]);
    }

    // Update the input matrix M
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            M[i*4 + j] = R[i][j];
        }
    }

    M[12] = tempT[0];
    M[13] = tempT[1];
    M[14] = tempT[2];

    // The bottom row remains unchanged ([0 0 0 1])
}


//Matrix this shape:
//[a11 a12 a13 a14]
//[a21 a22 a23 a24]
//[a31 a32 a33 a34]
//[0   0   0   1  ]
void matrixMultiplyRigid(double M1[], double M2[], double outM[], unsigned int x1, unsigned int y1, unsigned int x2, unsigned int y2, bool verbose=false) {
    // Validate matrix dimensions
    if (x1 != y2) {
        if (verbose) {
            std::cerr << "Not possible to multiply\n";
        }
        return;
    }

    for (unsigned int i = 0; i < 3; i++) {
        for (unsigned int j = 0; j < x2; j++) {
            outM[j + x2 * i] = M1[4 * i] * M2[j] + M1[4 * i + 1] * M2[j + x2] + M1[4 * i + 2] * M2[j + 2 * x2] + M1[4 * i + 3] * M2[j + 3 * x2];
        }
    }

    for (unsigned int j = 0; j < x2; j++) {
        outM[j + 3 * x2] = M2[j + 3 * x2];
    }

    // Printing for verbose mode
    if (verbose) {
        printf("\nThe product of entered matrices is:\n");
        for (unsigned int i = 0; i < y1; i++) {
            for (unsigned int j = 0; j < x2; j++) {
                printf("%f\t", outM[j + x2 * i]);
            }
            printf("\n");
        }
    }
}






template<typename T,typename U>
void backprojectToVolumeRealSpace (T* inProjection, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, double Phi, double Theta, double Psi, double tx, double ty, T* mask3D=NULL, T* externalReferenceMap=NULL, T* locresMap=NULL, double maskThreshold = 0.9) {

    //Phi=0; Theta=0; Psi=0;

    double M[16];
    double C[4];
    double CT[4];
    long int nxy = (long int)nx * ny;
    long int nxyz = nxy*nz;
    transformMatrixRelion(M, Psi, Theta, Phi, 0, 0, 0);
    inverseMatrix(M);
    float * tmpI =new float[nxyz];
    float * newX =new float[nxy];
    float * newY =new float[nxy];

    for (long int ij=0; ij<nxy; ij++) {
        //outI[ij] = 0;
        tmpI[ij] = 0;
        newX[ij]=0;
        newY[ij]=0;
    }
    double sum=0;
    const double avgPixNumber=(nx+ny+nz)/3;
    double nz2 = ((double)nz) / 2.0;
    double ny2 = ((double)ny) / 2.0;
    double nx2 = ((double)nx) / 2.0;
    C[3] = 1;
    C[2]=nz/2.0;
    for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              long int nyjj = ny*jj;
              for ( long int ii=0;ii<(long int)nx;ii++){
                        long int ij=ii+nyjj;
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;


                        /*
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                          tmpI[newX0+newY0*nx]+=c*d; //

                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                          tmpI[newX1+newY0*nx]+=a*d;
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                          tmpI[newX0+newY1*nx]+=b*c;
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                          tmpI[newX1+newY1*nx]+=a*b;
                        }*/
              }
      }

/*
      for (int ij=0; ij<nxy;ij++){
        if (tmpI[ij]>0.00001)
          outI[ij]=outI[ij]*avgPixNumber/tmpI[ij];
      }
*/

/*
    for ( long int kk=0;kk<(long int)nz;kk++){
      long int kknxy = nxy*kk;
      for ( long int jj=0;jj<(long int)ny;jj++){
        long int jjnx = jj*nx;
        long int jjnx_kknxy = jjnx+kknxy;
        for ( long int ii=0;ii<(long int)nx;ii++){
          long int ii_jjnx = ii+jjnx;
          long int ii_jjnx_kknxy = ii+jjnx_kknxy;
          long int tmpX=floor(newX[ii_jjnx]);
          long int tmpY=floor(newY[ii_jjnx]);
          if (tmpX>=0&&tmpX<nx && tmpY>=0&&tmpY<ny){
            outI[ii_jjnx_kknxy]= inProjection[tmpX+tmpY*nx];
          }
        }
      }
    }

*/

    //writeMrcImage("counter.mrc", tmpI, nx,ny,nz);
    delete []tmpI; 
    delete []newX; 
    delete []newY; 


}



template<typename T,typename U>
void ProjectVolumeRealSpace(T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, 
                            double Phi, double Theta, double Psi, double tx, double ty, 
                            T* mask3D=NULL, double maskThreshold = 0.9) {


//Phi=0; Theta=0; Psi=0;

    double M[16];
    double C[4];
    double CT[4];
    long int nxy = (long int)nx * ny;
    long int nxyz = nxy*nz;
    transformMatrixRelion(M, Psi, Theta, Phi, 0, 0, 0);
    inverseMatrix(M);
    float * tmpI =new float[nxy];

    for (long int ij=0; ij<nxy; ij++) {
        outI[ij] = 0;
        tmpI[ij] = 0;
    }
    double sum=0;
    const double avgPixNumber=(nx+ny+nz)/3;


    double nz2 = ((double)nz) / 2.0;
    double ny2 = ((double)ny) / 2.0;
    double nx2 = ((double)nx) / 2.0;
    C[3] = 1;


          for ( long int kk=0, ijk=0;kk<(long int)nz;kk++){
            C[2]=kk-(nz)/2.0;
            for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              for ( long int ii=0;ii<(long int)nx;ii++,ijk++){
                if (!mask3D || (mask3D && mask3D[ijk] > maskThreshold)) {
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                          tmpI[newX0+newY0*nx]+=c*d; //

                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                          tmpI[newX1+newY0*nx]+=a*d;
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                          tmpI[newX0+newY1*nx]+=b*c;
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                          tmpI[newX1+newY1*nx]+=a*b;
                        }
               }
              }
            }
          }
          for (int ij=0; ij<nxy;ij++){
            if (tmpI[ij]>0.00001)
              outI[ij]=outI[ij]*avgPixNumber/tmpI[ij];
          }

    //writeMrcImage("counter.mrc", tmpI, nx,ny,nz);
    delete []tmpI; 

}








template<typename T,typename U>
void ProjectVolumeRealSpace_original (T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, double Phi, double Theta, double Psi, double tx, double ty, T* mask3D=NULL, double maskThreshold = 0.9){

  double M[16];
  double C[4];
  double CT[4];
  //long int nxyz=(long int)nx*ny*nz;
  long int nxy=(long int)nx*ny;
  transformMatrixRelion(M, Psi, Theta, Phi, tx, ty, 0);
  inverseMatrix(M);


  //initialize
  for ( long int ij=0; ij<nxy; ij++){
   outI[ij]=0;
  }

  // ********
  //  3D MASK
  if (mask3D){
          double nz2=  (double)nz/2.0;
          double ny2=  (double)ny/2.0;
          double nx2=  (double)nx/2.0;
          C[3]=1;
          for ( long int kk=0, ijk=0;kk<(long int)nz;kk++){
            C[2]=kk-(nz)/2.0;
            for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              for ( long int ii=0;ii<(long int)nx;ii++,ijk++){
                if ( mask3D[ijk] > maskThreshold ){
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                        }
               }
              }
            }
          }
  }
  // *******************
  // NO MASK
  else if(!mask3D){
          double nz2=  (double)nz/2.0;
          double ny2=  (double)ny/2.0;
          double nx2=  (double)nx/2.0;
          C[3]=1;
          for ( long int kk=0, ijk=0;kk<(long int)nz;kk++){
            C[2]=kk-(nz)/2.0;
            for ( long int jj=0;jj<(long int)ny;jj++){
              C[1]=jj-(ny)/2.0;
              for ( long int ii=0;ii<(long int)nx;ii++,ijk++){
                        C[0]=ii-(nx)/2.0;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        double a = X-floor(X);
                        double b = Y-floor(Y);
                        double c = 1.0-a;
                        double d = 1.0-b;                
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]+=c*d*inI[ijk]; //
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]+=a*d*inI[ijk];//
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]+=b*c*inI[ijk];//
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]+=a*b*inI[ijk];//
                        }
              }
            }
          }
  }
}








template<typename T,typename U>
void ProjectMaskRealSpace(T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz, double Phi, double Theta, double Psi, double tx, double ty, double maskThreshold = 1){

  double M[16];
  double C[4];
  double CT[4];
  //long int nxyz=(long int)nx*ny*nz;
  long int nxy=(long int)nx*ny;
  transformMatrixRelion(M, Psi, Theta, Phi, tx, ty, 0);
  inverseMatrix(M);

  //initialize
  for (long int ij=0; ij<nxy; ij++){
   outI[ij]=0;
  }
  T* mask3D=inI;

  // ********
  //  3D MASK
          double nz2=  (double)nz/2.0;
          double ny2=  (double)ny/2.0;
          double nx2=  (double)nx/2.0;
          C[3]=1;
          for (unsigned long int kk=0, ijk=0;kk<nz;kk++){
            C[2]=kk-nz2;
            for ( unsigned long int jj=0;jj<ny;jj++){
              C[1]=jj-ny2;
              for (unsigned long int ii=0;ii<nx;ii++,ijk++){
                if (mask3D[ijk] > maskThreshold ){
                        C[0]=ii-nx2;
                        matrixMultiplication(M,C,CT,4,4,1,4);
                        
                        //bilinear interpolation
                        double X = (nx2+CT[0])+tx;
                        double Y = (ny2+CT[1])+ty;
                        long int newX0 = (long int)floor(X);
                        long int newX1 = (long int)ceil(X);
                        long int newY0 = (long int)floor(Y);
                        long int newY1 = (long int)ceil(Y);
                        
                        if(newX0>=0 && newX0<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX0+newY0*nx]=1; //
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY0>=0 && newY0<(long int)ny){
                          outI[newX1+newY0*nx]=1;//
                        }
                        if(newX0>=0 && newX0<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX0+newY1*nx]=1;//
                        }
                        if(newX1>=0 && newX1<(long int)nx && newY1>=0 && newY1<(long int)ny){
                          outI[newX1+newY1*nx]=1;//
                        }
               }
              }
            }
          }
}





#endif