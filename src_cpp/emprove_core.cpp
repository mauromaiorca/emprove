#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <cstdio>
#include <string.h>
#include <errno.h>
#include <stdlib.h>
#include <vector>
#include "mrcIO.h"
#include "matrix_computations.h"
#include "ctfLibs.h"
#include "derivative_libs.h"
#include "scores.h"
#include "imageProcessingLib.h"
#include "euler_libs.h"


// ******************************************
// adapted from: https://gist.github.com/rjzak/5681680
// INTERFACE FUNCTION
// PyObject -> Vector
std::vector<float> listTupleToVector(PyObject* incoming) {
	std::vector<float> data;
	if (PyTuple_Check(incoming)) {
		for(Py_ssize_t i = 0; i < PyTuple_Size(incoming); i++) {
			PyObject *value = PyTuple_GetItem(incoming, i);
			data.push_back( PyFloat_AsDouble(value) );
		}
	} else {
		if (PyList_Check(incoming)) {
			for(Py_ssize_t i = 0; i < PyList_Size(incoming); i++) {
				PyObject *value = PyList_GetItem(incoming, i);
				data.push_back( PyFloat_AsDouble(value) );
			}
		} else {
			throw std::logic_error("Passed PyObject pointer was not a list or tuple!");
		}
	}
	return data;
}



template<typename T>
PyObject* vectorToList(const std::vector<T> &data) {
  PyObject* listObj = PyList_New( data.size() );
	if (!listObj) throw std::logic_error("Unable to allocate memory for Python list");
	for (Py_ssize_t i = 0; i < data.size(); i++) {
		PyObject *num = PyFloat_FromDouble( (double) data[i]);
		if (!num) {
			Py_DECREF(listObj);
			throw std::logic_error("Unable to allocate memory for Python list");
		}
		PyList_SET_ITEM(listObj, i, num);
	}
	return listObj;
}






// ******************************************
// MAIN INTERFACE

// ******************************************
// projectMap
static PyObject * emprove_projectMap(PyObject *self, PyObject *args){

    PyObject * data_map3D_in;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    double phi;
    double theta;
    double psi;
    double centerX;
    double centerY;
    double centerZ;
    //printf("Pre Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);

    if (!PyArg_ParseTuple(args, "Okkkdddddd", &data_map3D_in, &nx,&ny,&nz, &phi,&theta,&psi,&centerX,&centerY,&centerZ))
        return NULL;
 
      Py_ssize_t nxyz=nx*ny*nz;
      Py_ssize_t nxy=nx*ny;

      if ( PyList_Check(data_map3D_in) ){
        if ( nxyz == PyList_Size(data_map3D_in) ){
            float * mapI = new float [nxyz];
            float * RI = new float [nxy];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *value = PyList_GetItem(data_map3D_in, i);
                mapI[i] = PyFloat_AsDouble(value);
            }
            //writeMrcImage("test.mrc", mapI,  nx,  ny,  nz,  1);
            ProjectVolumeRealSpace(mapI, RI,  nx,  ny,  nz, -phi, -theta, -psi, -centerX, -centerY);
            //ProjectVolumeRealSpace(mapI, RI,  nx,  ny,  nz, fmod(360.0 - phi, 360.0),   180.0 - (theta + 90.0) - 90.0, fmod(360.0 - psi, 360.0), -centerX, -centerY);

            //writeMrcImage("test.mrc", RI,  nx,  ny,  1,  1);
            delete [] mapI;
            PyObject *lst = PyList_New(nxy);
            if (!lst)
                return NULL;
            for (Py_ssize_t i = 0; i < nxy; i++) {
                PyObject *num = PyFloat_FromDouble(RI[i]);
                if (!num) {
                    Py_DECREF(lst);
                    throw std::logic_error("Unable to allocate memory for Python list");
                    return NULL;
                }
                PyList_SET_ITEM(lst, i, num);   // reference to num stolen
            }
            delete [] RI;
            return lst;
        }else{
            //printf('ERROR: mismatch dimensions');
            return NULL; //mismatch dimensions
        }
      }else{
            //printf('ERROR: map not found');
            return NULL; //map not found
        }
}


static PyObject * emprove_backprojectParticles(PyObject *self, PyObject *args){

    PyObject * data_particle2D_in;
    PyObject * data_map3D_out;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    double phi;
    double theta;
    double psi;
    double centerX;
    double centerY;
    double centerZ;
    //printf("Pre Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);

    if (!PyArg_ParseTuple(args, "OOkkkdddddd", &data_particle2D_in, &data_map3D_out, &nx,&ny,&nz, &phi,&theta,&psi,&centerX,&centerY,&centerZ))
        return NULL;
 
      Py_ssize_t nxyz=nx*ny*nz;
      Py_ssize_t nxy=nx*ny;

      if ( PyList_Check(data_particle2D_in) && PyList_Check(data_map3D_out)){
        if ( nxy == PyList_Size(data_particle2D_in) && nxyz == PyList_Size(data_map3D_out) ){
            float * ProjI = new float [nxy];
            float * mapI = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxy; i++) {
                PyObject *value = PyList_GetItem(data_map3D_out, i);
                ProjI[i] = PyFloat_AsDouble(value);
            }
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *value = PyList_GetItem(data_map3D_out, i);
                mapI[i] = PyFloat_AsDouble(value);
            }

            //writeMrcImage("test.mrc", mapI,  nx,  ny,  nz,  1);
            backprojectToVolumeRealSpace(ProjI,  mapI, nx,  ny,  nz, phi, theta, psi, centerX, centerY);
            //ProjectVolumeRealSpace(mapI, ProjI,  nx,  ny,  nz, fmod(360.0 - phi, 360.0),   180.0 - (theta + 90.0) - 90.0, fmod(360.0 - psi, 360.0), -centerX, -centerY);

            //writeMrcImage("test.mrc", ProjI,  nx,  ny,  1,  1);
            delete [] ProjI;
            PyObject *lst = PyList_New(nxy);
            if (!lst)
                return NULL;
            for (Py_ssize_t i = 0; i < nxy; i++) {
                PyObject *num = PyFloat_FromDouble(mapI[i]);
                if (!num) {
                    Py_DECREF(lst);
                    throw std::logic_error("Unable to allocate memory for Python list");
                    return NULL;
                }
                PyList_SET_ITEM(lst, i, num);   // reference to num stolen
            }
            delete [] mapI;
            return lst;
        }else{
            //printf('ERROR: mismatch dimensions');
            return NULL; //mismatch dimensions
        }
      }else{
            //printf('ERROR: map not found');
            return NULL; //map not found
        }
        return NULL;
}



// ******************************************
// rotAverage
static PyObject * emprove_rotAverage(PyObject *self, PyObject *args){

    PyObject * data_map3D_in;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    if (!PyArg_ParseTuple(args, "Okkk", &data_map3D_in, &nx,&ny,&nz))
        return NULL;
 
    Py_ssize_t nxyz=nx*ny*nz;
    Py_ssize_t nxy=nx*ny;

    std::vector<double> returnVector;
    if ( PyList_Check(data_map3D_in) ){
        if ( nxyz == PyList_Size(data_map3D_in) ){
            float * mapI = new float [nxyz];
                for(Py_ssize_t i = 0; i < nxyz; i++) {
                    PyObject *value = PyList_GetItem(data_map3D_in, i);
                    mapI[i] = PyFloat_AsDouble(value); 
                }
            returnVector=rotAverage(mapI,nx,ny,nz);
            delete []mapI;
        }
    }
    //std::cerr<<"vector ("<<returnVector.size()<<")   ="<<returnVector[0]<<","<<returnVector[1]<<"\n";
    //return values
    PyObject *lst = PyList_New(returnVector.size());
    if (!lst)
        return NULL;
    for (Py_ssize_t i = 0; i < returnVector.size(); i++) {
        PyObject *num = PyFloat_FromDouble(returnVector[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    return lst;
}



// ******************************************
// amplitudeNormalization
static PyObject * emprove_amplitudeNormalization(PyObject *self, PyObject *args){

    PyObject * data_map3D_mapInOut;
    PyObject * data_map3D_amplitudeToReplace;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    if (!PyArg_ParseTuple(args, "OOkkk", &data_map3D_mapInOut, &data_map3D_amplitudeToReplace, &nx,&ny,&nz))
        return NULL;
 
    Py_ssize_t nxyz=nx*ny*nz;
    Py_ssize_t nxy=nx*ny;
    float * ampI = new float [nxyz];

    std::vector<double> returnVector;
    if ( PyList_Check(data_map3D_mapInOut) && PyList_Check(data_map3D_amplitudeToReplace) ){
        if ( nxyz == PyList_Size(data_map3D_mapInOut) && nxyz == PyList_Size(data_map3D_mapInOut) ){
            float * ampReplacementI = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *valueI = PyList_GetItem(data_map3D_mapInOut, i);
                ampI[i] = PyFloat_AsDouble(valueI); 
                PyObject *valueR = PyList_GetItem(data_map3D_amplitudeToReplace, i);
                ampReplacementI[i] = PyFloat_AsDouble(valueR);
            }
            amplitudeNormalization(ampI,ampReplacementI,nx,ny,nz);
            delete [] ampReplacementI;
        }
    }
    //std::cerr<<"vector ("<<returnVector.size()<<")   ="<<returnVector[0]<<","<<returnVector[1]<<"\n";
    //return values
    PyObject *lst = PyList_New(nxyz);
    if (!lst)
        return NULL;
    for (Py_ssize_t i = 0; i < nxyz; i++) {
        PyObject *num = PyFloat_FromDouble(ampI[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] ampI;
    return lst;
}





// ******************************************
// projectMap
static PyObject * emprove_projectMask(PyObject *self, PyObject *args){

    PyObject * data_map3D_in;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    double phi;
    double theta;
    double psi;
    double centerX;
    double centerY;
    double centerZ;
    double maskThreshold=0.9;
    //printf("Pre Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);

    if (!PyArg_ParseTuple(args, "Okkkdddddd|d", &data_map3D_in, &nx,&ny,&nz, &phi,&theta,&psi,&centerX,&centerY,&centerZ,&maskThreshold))
        return NULL;
 
      Py_ssize_t nxyz=nx*ny*nz;
      Py_ssize_t nxy=nx*ny;

      if ( PyList_Check(data_map3D_in) ){
        if ( nxyz == PyList_Size(data_map3D_in) ){
            float * mapI = new float [nxyz];
            float * RI = new float [nxy];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *value = PyList_GetItem(data_map3D_in, i);
                mapI[i] = PyFloat_AsDouble(value);
            }
            //writeMrcImage("test.mrc", mapI,  nx,  ny,  nz,  1);
            //ProjectVolumeRealSpace(mapI, RI,  nx,  ny,  nz, -phi, -theta, -psi, -centerX, -centerY);
            ProjectMaskRealSpace(mapI, RI,  nx,  ny,  nz, -phi, -theta, -psi, -centerX, -centerY, maskThreshold);

            //writeMrcImage("test.mrc", RI,  nx,  ny,  1,  1);
            delete [] mapI;
            PyObject *lst = PyList_New(nxy);
            if (!lst)
                return NULL;
            for (Py_ssize_t i = 0; i < nxy; i++) {
                PyObject *num = PyFloat_FromDouble(RI[i]);
                if (!num) {
                    Py_DECREF(lst);
                    return NULL;
                }
                PyList_SET_ITEM(lst, i, num);   // reference to num stolen
            }
            delete [] RI;
            return lst;
        }else{
            //printf('ERROR: mismatch dimensions');
            return NULL; //mismatch dimensions
        }
      }else{
            //printf('ERROR: map not found');
            return NULL; //map not found
        }

}



// *********************************
// project (with file)
//
static PyObject * emprove_project(PyObject *self, PyObject *args){

    const char *map3D_in;
    const char *map3D_out;
    double phi;
    double theta;
    double psi;
    double centerX;
    double centerY;
    double centerZ;
    //printf("Pre Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);

    if (!PyArg_ParseTuple(args, "ssdddddd", &map3D_in, &map3D_out, &phi,&theta,&psi,&centerX,&centerY,&centerZ))
        return NULL;
    printf("Map_in=%s Map_out=%s phi=%f theta=%f psi=%f centerX=%f centerY=%f centerZ=%f\n", map3D_in, map3D_out, phi,theta,psi,centerX,centerY,centerZ);
 

      MRCHeader mapHeader;
      readHeaderMrc(map3D_in, mapHeader);
      unsigned long int nx=mapHeader.nx;
      unsigned long int ny=mapHeader.ny;
      unsigned long int nz=mapHeader.nz;
      unsigned long int nxyz=nx*ny*nz;
      unsigned long int nxy=nx*ny;

      float * referenceMapI = new float[nxyz];
      float * projectionI = new float[nxy];
      readMrcImage(map3D_in, referenceMapI, mapHeader);
      ProjectVolumeRealSpace(referenceMapI, projectionI,  nx,  ny,  nz, -phi, -theta, -psi, -centerX, -centerY);
      // (referenceMapI, projectionI,  nx,  ny,  nz, fmod(360.0 - phi, 360.0), fmod(180.0 - theta, 180.0), fmod(360.0 - psi, 360.0), -centerX, -centerY);
      writeMrcImage(map3D_out, projectionI, nx, ny, 1);
      delete [] referenceMapI;
      delete [] projectionI;
    //gets a map and projects it according to coordinates, from center
//    printf("Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);
    return Py_BuildValue("d",0);
}

// ****************************
// WriteMRC
static PyObject * emprove_WriteMRC(PyObject *self, PyObject *args){

    PyObject * dataIn;
    char *map3D_out;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    double aPix;

    if (!PyArg_ParseTuple(args, "Oskkkd", &dataIn, &map3D_out, &nx, &ny, &nz, &aPix)){
        std::cerr<<"WARNING: unable to write "<< map3D_out << " (incorrect input parameters)\n";
       return NULL;
    }
    Py_ssize_t nxyz = nx * ny * nz;
    std::cerr<<"Writing MRC file: "<<  map3D_out << "   size:"<<nxyz;
    std::cerr<<"  (nx="<<nx;
    std::cerr<<"  ny="<<ny;
    std::cerr<<"  nz="<<nz<<")\n";
    //std::cerr<<"  PyList_Size(dataIn)="<<PyList_Size(dataIn)<<"\n";
    

	if (PyList_Check(dataIn)) {
        if (nxyz == PyList_Size(dataIn)){
            float * dataI = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *value = PyList_GetItem(dataIn, i);
                dataI[i] = PyFloat_AsDouble(value);
            }
            writeMrcImage(map3D_out, dataI,  nx,  ny,  nz,  aPix);
            delete [] dataI;
        }else{
            std::cerr<<"WARNING: unable to write "<< map3D_out << " (incorrect size)\n";
            if (PyList_Size(dataIn)==nx){
                std::cerr<<"        Hint: input should be a flattened 1D array\n";

            }
        }
	}else{
        std::cerr<<"WARNING: unable to write "<< map3D_out << "  (incorrect input datatype, hint:  it should be a python list, e.g. use .tolist() fuction)\n";
    }
    return Py_BuildValue("d",0);
}


// ****************************
// emprove_MaskedComparison
static PyObject * emprove_EqualizedParticlesRank(PyObject *self, PyObject *args){
    PyObject * object_phiListParticle;
    PyObject * object_thetaListParticle;
    PyObject * object_scores;
    PyObject * object_randomSubset;
    unsigned long int numViews = 200;
    if (!PyArg_ParseTuple(args, "OOOO|k", &object_phiListParticle, &object_thetaListParticle, &object_scores, &object_randomSubset, &numViews))
        return NULL;
    //printf("object_phiListParticle size=%d\n",PyList_Size(object_phiListParticle));
    //printf("object_thetaListParticle size=%d\n",PyList_Size(object_thetaListParticle));
    //printf("object_scores size=%d\n",PyList_Size(object_scores));
    //printf("object_randomSubset size=%d\n",PyList_Size(object_randomSubset));
    std::vector<double> phiListParticle;
    std::vector<double> thetaListParticle;
    std::vector<double> scores;
    std::vector<double> randomSubset;
    std::vector<double> result;
	if (PyList_Check(object_phiListParticle) && PyList_Check(object_thetaListParticle) && PyList_Check(object_scores) && PyList_Check(object_randomSubset) ) {
            for(Py_ssize_t i = 0; i < PyList_Size(object_phiListParticle); i++) {
                phiListParticle.push_back( PyFloat_AsDouble(PyList_GetItem(object_phiListParticle, i) ));                
                thetaListParticle.push_back( PyFloat_AsDouble(PyList_GetItem(object_thetaListParticle, i) ));
                scores.push_back( PyFloat_AsDouble(PyList_GetItem(object_scores, i) ));
                randomSubset.push_back( PyFloat_AsDouble(PyList_GetItem(object_randomSubset, i) ));
            }
            result = scoreNormalizationHomogeneousViews_halfMaps(  phiListParticle, thetaListParticle, scores, randomSubset, numViews);
    }else{
        printf("ERROR: passing object\n");
    }

    PyObject *lst = PyList_New(result.size());
    if (!lst){
        printf("ERROR: memory error, creating lst\n");
        return NULL;
    }
    for (Py_ssize_t i = 0; i < (Py_ssize_t)result.size(); i++) {
        PyObject *rank = PyFloat_FromDouble( result[i] );
        if (!rank) {
            Py_DECREF(lst);
            printf("ERROR: memory error, creating Py_DECREF\n");
            return NULL;
        }
        PyList_SET_ITEM ( lst, i, PyFloat_FromDouble( result[i] ) );   // reference to num stolen
    }


    return lst;
}


// ****************************
// emprove_MaskedComparison
static PyObject * emprove_GetEulerClassGroup(PyObject *self, PyObject *args){
    PyObject * object_phiListParticle;
    PyObject * object_thetaListParticle;
    unsigned long int numViews;
    if (!PyArg_ParseTuple(args, "OOk", &object_phiListParticle, &object_thetaListParticle, &numViews))
        return NULL;

    std::vector<double> phiListParticle;
    std::vector<double> thetaListParticle;
    std::vector<std::string> result;
	if (PyList_Check(object_phiListParticle) && PyList_Check(object_thetaListParticle)  ) {
            for(Py_ssize_t i = 0; i < PyList_Size(object_phiListParticle); i++) {
                phiListParticle.push_back( PyFloat_AsDouble(PyList_GetItem(object_phiListParticle, i) ));                
                thetaListParticle.push_back( PyFloat_AsDouble(PyList_GetItem(object_thetaListParticle, i) ));
            }
            result = getEulerClassGroup(  phiListParticle, thetaListParticle, numViews);
    }else{
        printf("ERROR: passing object\n");
    }

    PyObject *lst = PyList_New(result.size());
    if (!lst){
        printf("ERROR: memory error, creating lst\n");
        return NULL;
    }
    for (Py_ssize_t i = 0; i < (Py_ssize_t)result.size(); i++) {
        
        PyObject *rank = PyUnicode_FromString( result[i].c_str() );
        if (!rank) {
            Py_DECREF(lst);
            printf("ERROR: memory error, creating Py_DECREF\n");
            return NULL;
        }
        PyList_SET_ITEM(lst, i, rank);   // reference to num stolen
        
    }
    return lst;
}

// ****************************
// emprove_MaskedComparison
static PyObject * emprove_SCI_distance(PyObject *self, PyObject *args){
    PyObject * object_I1;
    PyObject * object_I2;
    PyObject * object_M  = NULL;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    char * sigmaStr = (char *)"1";


    if (!PyArg_ParseTuple(args, "OOkkks|O", &object_I1, &object_I2, &nx, &ny, &nz, &sigmaStr, &object_M ))
        return NULL;
    double sigma=atof(sigmaStr);
    Py_ssize_t nxyz = nx * ny * nz;
    double score = 0;

    
	if ( PyList_Check(object_I1) && PyList_Check(object_I2) ) {
        if (nxyz == PyList_Size(object_I1) && nxyz == PyList_Size(object_I2) ){
            float * I1 = new float [nxyz];
            float * I2 = new float [nxyz];
            float * M = NULL;// = new float [nxyz];
            if (object_M) if (PyList_Size(object_M) > 0 && PyList_Check(object_M) )
                 M = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *valueI1 = PyList_GetItem(object_I1, i);
                I1[i] = PyFloat_AsDouble(valueI1);

                PyObject *valueI2 = PyList_GetItem(object_I2, i);
                I2[i] = PyFloat_AsDouble(valueI2);
                
                if ( M ){
                    PyObject *valueM = PyList_GetItem(object_M, i);
                    M[i] = PyFloat_AsDouble(valueM);
                }
            }

            double minMaskValue = 0.1;
            score=SCI((float * )I1, (float *) I2, nx, ny, nz, sigma, (float *)M, minMaskValue);

            if ( M ){
                //writeMrcImage("mask_before.mrc",  M, nx,  ny,  1, 1);
                blurMaskEdgeOnImage(I1, M, I1, 0.8, 3, nx,  ny);
                blurMaskEdgeOnImage(I2, M, I2, 0.8, 3, nx,  ny);
            }

            if (I1) delete [] I1;
            if (I2) delete [] I2;
            if (M) delete [] M;
        }
	}
    
    return Py_BuildValue("d",score);

}




// ****************************
// emprove_MaskedComparison
static PyObject * emprove_SDIM(PyObject *self, PyObject *args){
    PyObject * object_I1;
    PyObject * object_I2;
    PyObject * object_M  = NULL;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    char * sigmaStr = (char *)"1";


    if (!PyArg_ParseTuple(args, "OOkkks|O", &object_I1, &object_I2, &nx, &ny, &nz, &sigmaStr, &object_M ))
        return NULL;

    double sigma=atof(sigmaStr);
    Py_ssize_t nxyz = nx * ny * nz;
    double score = 0;

    PyObject *lst = PyList_New(nxyz);
    if (!lst)
        return NULL;


	if ( PyList_Check(object_I1) && PyList_Check(object_I2) ) {
        if (nxyz == PyList_Size(object_I1) && nxyz == PyList_Size(object_I2) ){
            float * I1 = new float [nxyz];
            float * I2 = new float [nxyz];
            float * outI = new float [nxyz];
            float * M = NULL;// = new float [nxyz];
            if (object_M) if (PyList_Size(object_M) > 0 && PyList_Check(object_M) )
                 M = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                outI[i]=0;
                PyObject *valueI1 = PyList_GetItem(object_I1, i);
                I1[i] = PyFloat_AsDouble(valueI1);

                PyObject *valueI2 = PyList_GetItem(object_I2, i);
                I2[i] = PyFloat_AsDouble(valueI2);
                
                if ( M ){
                    PyObject *valueM = PyList_GetItem(object_M, i);
                    M[i] = PyFloat_AsDouble(valueM);
                }
            }

            double minMaskValue = 0.1;
            SDIM((float * )I1, (float *) I2, (float *) outI, nx, ny, nz, sigma, (float *)M, minMaskValue);
            //lst = PyList_New(nxyz);
            for (Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *num = PyFloat_FromDouble(outI[i]);
                if (!num) {
                    Py_DECREF(lst);
                    return NULL;
                }
                PyList_SET_ITEM(lst, i, num);   // reference to num stolen
            }

            if ( M ){
                //writeMrcImage("mask_before.mrc",  M, nx,  ny,  1, 1);
                blurMaskEdgeOnImage(I1, M, I1, 0.8, 3, nx,  ny);
                blurMaskEdgeOnImage(I2, M, I2, 0.8, 3, nx,  ny);
            }

            if (I1) delete [] I1;
            if (I2) delete [] I2;
            if (M) delete [] M;
            if (outI) delete [] outI;
        }
	}
    

    return lst;

//   return Py_BuildValue("d",0);

}



// ****************************
// emprove_MaskedComparison
static PyObject * emprove_maskedNormalizeInfo(PyObject *self, PyObject *args){
    PyObject * object_I1;
    PyObject * object_I2;
    PyObject * object_M;
    Py_ssize_t nxy;
    double maskThreshold = 0.1;
    if (!PyArg_ParseTuple(args, "OOOk", &object_I1, &object_I2, &object_M, &nxy))
        return NULL;

    std::vector<double> output;
            float * I1 = new float [nxy];
            float * I2 = new float [nxy];
            float * M = new float [nxy];

            for(Py_ssize_t i = 0; i < nxy; i++) {

                PyObject *valueM = PyList_GetItem(object_M, i);
                M[i] = PyFloat_AsDouble(valueM);
                if (M[i] > maskThreshold){
                    PyObject *valueI1 = PyList_GetItem(object_I1, i);
                    I1[i] = PyFloat_AsDouble(valueI1);
                    PyObject *valueI2 = PyList_GetItem(object_I2, i);
                    I2[i] = PyFloat_AsDouble(valueI2);
                }else{
                    I1[i] = 0;
                    I2[i] = 0;
                }
            }
            output=normalize(I1, I2, M, nxy);
    delete [] I1;
    delete [] I2;
    delete [] M;

    PyObject *lst = PyList_New(4);
    if (output.size()==4){
        PyList_SET_ITEM(lst, 0, PyFloat_FromDouble(output[0]));
        PyList_SET_ITEM(lst, 1, PyFloat_FromDouble(output[1]));
        PyList_SET_ITEM(lst, 2, PyFloat_FromDouble(output[2]));
        PyList_SET_ITEM(lst, 3, PyFloat_FromDouble(output[3]));
    }else{
        PyList_SET_ITEM(lst, 0, PyFloat_FromDouble(0));
        PyList_SET_ITEM(lst, 1, PyFloat_FromDouble(0));
        PyList_SET_ITEM(lst, 2, PyFloat_FromDouble(0));
        PyList_SET_ITEM(lst, 3, PyFloat_FromDouble(0));
    }
    return lst;
}



// ****************************
// emprove_derivative
static PyObject * emprove_derivative1D(PyObject *self, PyObject *args){
    PyObject * object_I1;
    double blurSigma;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    int DerivativeOrder=0;
    int direction=0;

    if (!PyArg_ParseTuple(args, "Odkkk|ii", &object_I1, &blurSigma, &nx, &ny, &nz, &DerivativeOrder, &direction))
        return NULL;
    if (nz!=nz) nz=1; //trick to check if nz is NaN (self-comparison *NOT* always evaluates to false)
    if (nz<1) nz=1;

    Py_ssize_t nxyz=nx*ny*nz;
    std::vector<double> output;
            float * dataI = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                    PyObject *valueI1 = PyList_GetItem(object_I1, i);
                    dataI[i] = PyFloat_AsDouble(valueI1);
            }
//    delete [] I1;

//flipImage( dataI, dataI,  nx,  ny,  nz, direction);
gaussRecursiveDerivatives1D ( blurSigma, nx, ny, nz, 1,  2, direction, DerivativeOrder, dataI);
//flipImage( dataI, dataI,  nx,  ny,  nz, direction);


    PyObject *lst = PyList_New(nxyz);
    if (!lst)
        return NULL;
    for (Py_ssize_t i = 0; i < nxyz; i++) {
        PyObject *num = PyFloat_FromDouble(dataI[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] dataI;
    return lst;
}


// ****************************
// emprove_createHalfMaps
static PyObject * emprove_createHalfMaps(PyObject *self, PyObject *args){
    PyObject * object_I1;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    unsigned long int mapNumber;

    if (!PyArg_ParseTuple(args, "Okkkk", &object_I1, &nx, &ny, &nz, &mapNumber))
        return NULL;
    Py_ssize_t nxyz=nx*ny*nz;
    float * dataI = new float [nxyz];
    for(Py_ssize_t i = 0; i < nxyz; i++) {
        PyObject *valueI1 = PyList_GetItem(object_I1, i);
        dataI[i] = PyFloat_AsDouble(valueI1);
    }
//    delete [] I1;
//    gaussRecursiveDerivatives2D ( blurSigma, nx, ny,  1, 1,  2, DerivativeOrder, dataI);
    float * maskI = new float [nxyz];
    float * outMapI1 = new float [nxyz];
    float * outMapI2 = new float [nxyz];
    unsigned long int nxy=nx*ny;
    for (unsigned long int kk= 0; kk < nz; kk++){
        long int nxykk=nxy*kk;
        for (unsigned long int jj = 0; jj < ny; jj++){
            for (unsigned long int ii = 0; ii < nx; ii++){
                 long int idx = ii+nx*jj+ nxykk;
                 long int idxN = ii+nx*(jj-1)+nxykk;
                 long int idxS = ii+nx*(jj+1)+nxykk;
                 long int idxE = idx-1;
                 long int idxW = idx+1;
                 long int idxB = ii+nx*jj+nxy*(kk-1);
                 long int idxF = ii+nx*jj+nxy*(kk+1);


                if (jj==0) idxN = idxS;
                if (jj==ny-1) idxS = idxN;
                if (ii==0) idxE = idxW;
                if (ii==nx-1) idxW = idxE;
                if (kk==0) idxB = idxF;
                if (kk==nz-1) idxF = idxB;

                if (idxN<0)idxN=0;
                if (idxS<0)idxS=0;
                if (idxW<0)idxW=0;
                if (idxE<0)idxE=0;
                if (idxB<0)idxB=0;
                if (idxF<0)idxF=0;

                if (idxN>=nxyz)idxN=nxyz-1;
                if (idxS>=nxyz)idxS=nxyz-1;
                if (idxW>=nxyz)idxW=nxyz-1;
                if (idxE>=nxyz) idxE=nxyz-1;
                if (idxB>=nxyz)idxB=nxyz-1;
                if (idxF>=nxyz)idxF=nxyz-1;

                double value;
                if (nz<3){
                    value=dataI[idxN]+dataI[idxS]+dataI[idxE]+dataI[idxW];
                    value/=4;
                }else{
                    value=dataI[idxN]+dataI[idxS]+dataI[idxE]+dataI[idxW]+dataI[idxF]+dataI[idxB];
                    value/=6;
                }
                maskI[idx]=(ii+jj+kk+mapNumber)%2;
                if (maskI[idx]==0){
                    outMapI1[idx]=value;
                //    outMapI2[idx]=dataI[idx];
                }else{
                    outMapI1[idx]=dataI[idx];
                 //   outMapI2[idx]=value;
               }
            }
        }
    }


    PyObject *lst = PyList_New(nxyz);
    if (!lst)
        return NULL;
    for (Py_ssize_t i = 0; i < nxyz; i++) {
        PyObject *num = PyFloat_FromDouble(outMapI1[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    /*
    for (Py_ssize_t i = 0; i < nxyz; i++) {
        PyObject *num = PyFloat_FromDouble(outMapI2[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, nxyz+i-1, num);   // reference to num stolen
    }*/
    
    delete [] dataI;
    delete [] maskI;
    delete [] outMapI1;
    delete [] outMapI2;
    return lst;
}




// ****************************
// emprove_derivative
static PyObject * emprove_derivative(PyObject *self, PyObject *args){
    PyObject * object_I1;
    double blurSigma;
    unsigned long int nx;
    unsigned long int ny;
    int DerivativeOrder=0;

    if (!PyArg_ParseTuple(args, "Odkk|i", &object_I1, &blurSigma, &nx, &ny, &DerivativeOrder))
        return NULL;
    Py_ssize_t nxy=nx*ny;
    std::vector<double> output;
            float * dataI = new float [nxy];

            for(Py_ssize_t i = 0; i < nxy; i++) {
                    PyObject *valueI1 = PyList_GetItem(object_I1, i);
                    dataI[i] = PyFloat_AsDouble(valueI1);
            }
//    delete [] I1;
    gaussRecursiveDerivatives2D ( blurSigma, nx, ny,  1, 1,  2, DerivativeOrder, dataI);

    PyObject *lst = PyList_New(nxy);
    if (!lst)
        return NULL;
    for (Py_ssize_t i = 0; i < nxy; i++) {
        PyObject *num = PyFloat_FromDouble(dataI[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] dataI;
    return lst;
}


// ****************************
// emprove_maskedMapContrast
////    {"maskedMapContrast",  (PyCFunction)emprove_maskedMapContrast, METH_VARARGS, "image, mask"},
static PyObject * emprove_maskedMapContrast(PyObject *self, PyObject *args){
    PyObject * ImageObject; //GT
    PyObject * MaskObject; //GT
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    if (!PyArg_ParseTuple(args, "OOkkk", &ImageObject, &MaskObject, &nx, &ny, &nz))
        return NULL;
    Py_ssize_t nxyz = nx * ny * nz;
    std::vector<double> result;
	if (PyList_Check(ImageObject) && PyList_Check(MaskObject) ) {
        if (nxyz == PyList_Size(ImageObject) && nxyz == PyList_Size(MaskObject) ){
            float * I = new float [nxyz];
            float * M = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *valueI1 = PyList_GetItem(ImageObject, i);
                I[i] = PyFloat_AsDouble(valueI1);
                PyObject *valueM = PyList_GetItem(MaskObject, i);
                M[i] = PyFloat_AsDouble(valueM);
            }
            result=computeMaskedAnalysis ( I, M, nx, ny, nz);
            delete [] I;
            delete [] M;
        }
    }
    return Py_BuildValue("dd",result[0],result[1]);
}

// ****************************
// emprove_thresholdMap
//    {"thresholdMap",  (PyCFunction)emprove_thresholdMap, METH_VARARGS, "image, map"},
static PyObject * emprove_AutomaticThresholdMap(PyObject *self, PyObject *args){
    PyObject * ImageObject; //GT
    PyObject * MaskObject; //GT
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
     int numTreshold = 2;
    if (!PyArg_ParseTuple(args, "OOkkk|i", &ImageObject, &MaskObject, &nx, &ny, &nz, &numTreshold))
        return NULL;
    Py_ssize_t nxyz = nx * ny * nz;
    float * outputI = NULL;
	if (PyList_Check(ImageObject) && PyList_Check(MaskObject) ) {
        if (nxyz == PyList_Size(ImageObject) && nxyz == PyList_Size(MaskObject) ){
            float * I = new float [nxyz];
            float * M = new float [nxyz];
            outputI = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *valueI1 = PyList_GetItem(ImageObject, i);
                I[i] = PyFloat_AsDouble(valueI1);
                PyObject *valueM = PyList_GetItem(MaskObject, i);
                M[i] = PyFloat_AsDouble(valueM);
            }
            maskedThreshold ( I, M, outputI, nx, ny, nz,numTreshold);
            delete [] I;
            delete [] M;
        }
    }

   
    PyObject *lst = PyList_New(nxyz);
    if (!lst)
        return NULL;
    if (outputI){
        
        for (Py_ssize_t i = 0; i < nxyz; i++) {
            PyObject *num = PyFloat_FromDouble(outputI[i]);
            if (!num) {
                Py_DECREF(lst);
                return NULL;
            }
            PyList_SET_ITEM(lst, i, num);   // reference to num stolen
        }
        delete [] outputI;
        return lst;
    }else{
        return NULL;
    }


}


// ****************************
// emprove_MaskedComparison
static PyObject * emprove_MaskedImageComparison(PyObject *self, PyObject *args){

    PyObject * object_I1; //GT
    PyObject * object_I2;
    PyObject * object_M;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
    char * comparisonMethod =  (char *)"CC";
    char * preprocessingMethod =  (char *)"unprocessed";
    char * sigmaStr = (char *)"1";
    
    //score_labels=['_LRA_CC_unprocessed_simple','_LRA_CC_firstDerivativeX','_LRA_CC_firstDerivativeY','_LRA_CC_largestEigenvalueHessian']

    double score = 0;
    if (!PyArg_ParseTuple(args, "OOOkkk|sss", &object_I1, &object_I2, &object_M, &nx, &ny, &nz,&comparisonMethod, &preprocessingMethod, &sigmaStr))
        return NULL;
    double sigma=atof(sigmaStr);
    Py_ssize_t nxyz = nx * ny * nz;
    
    //printf("nx= %d  ny=%d  nz=%d   nxyz=%d", nx, ny, nz, nxyz);
	if (PyList_Check(object_I1) && PyList_Check(object_I2) && PyList_Check(object_M)) {
        if (nxyz == PyList_Size(object_I1) && nxyz == PyList_Size(object_I2) && nxyz == PyList_Size(object_M) ){
            float * I1 = new float [nxyz];
            float * I2 = new float [nxyz];
            float * M = new float [nxyz];
            for(Py_ssize_t i = 0; i < nxyz; i++) {
                PyObject *valueI1 = PyList_GetItem(object_I1, i);
                I1[i] = PyFloat_AsDouble(valueI1);

                PyObject *valueI2 = PyList_GetItem(object_I2, i);
                I2[i] = PyFloat_AsDouble(valueI2);

                PyObject *valueM = PyList_GetItem(object_M, i);
                M[i] = PyFloat_AsDouble(valueM);
            }

        if ( strcmp(comparisonMethod,"SCI") !=0 ){
            if ( strcmp(preprocessingMethod,"firstDerivativeX")==0 ){
                //double sigma=1;
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 0, 1, I1);
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 0, 1, I2);
            } else if ( strcmp(preprocessingMethod,"firstDerivativeY")==0 ){
                //double sigma=1;
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 1, 1, I1);
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 1, 1, I2);
            } else if ( strcmp(preprocessingMethod,"secondDerivativeX")==0 ){
                //double sigma=1;
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 0, 2, I1);
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 0, 2, I2);
            } else if ( strcmp(preprocessingMethod,"secondDerivativeY")==0 ){
                //double sigma=1;
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 1, 2, I1);
                gaussRecursiveDerivatives1D(sigma, nx, ny, nz, 1, 3, 1, 2, I2);
            } else if ( strcmp(preprocessingMethod,"blur")==0 ){
                //double sigma=1;
                gaussRecursiveDerivatives2D(sigma, nx, ny, nz, 1, 3, 0, I1);
                gaussRecursiveDerivatives2D(sigma, nx, ny, nz, 1, 3, 0, I2);
            } else if ( strcmp(preprocessingMethod,"largestEigenvalueHessian")==0 ){
                //double sigma=1;
                computeEigenImage (I1,  nx,  ny, 2, sigma);
                computeEigenImage (I2,  nx,  ny, 2, sigma);
            }
        }


            /*
            //debug:
            if ( strcmp(comparisonMethod,"CC")==0 && strcmp(preprocessingMethod,"unprocessed")==0 ){
                blurMaskEdgeOnImage(I1, M, I1, 0.8, 3, nx,  ny);
                blurMaskEdgeOnImage(I2, M, I2, 0.8, 3, nx,  ny);
                writeMrcImage("CC_I1.mrc", I1, nx, ny, 1);
                writeMrcImage("CC_I2.mrc", I2, nx, ny, 1);
                writeMrcImage("CC_M.mrc", M, nx, ny, 1);
            }
            */

            //writeMrcImage("mask_before.mrc",  M, nx,  ny,  1, 1);
            /* OK
            blurMaskEdgeOnImage(I1, M, I1, 0.6, 5, nx,  ny);
            blurMaskEdgeOnImage(I2, M, I2, 0.6, 5, nx,  ny);
            */

            //writeMrcImage(map3D_out, dataI,  nx,  ny,  nz,  aPix);
            //printf("options: %s   %s   sigma=%f\n",comparisonMethod,preprocessingMethod,sigma);
            if ( strcmp(comparisonMethod,"CC")==0 ){
                //printf("Cross correlation\n");
                score=crossCorrelationDistance(I1, I2,  nxyz, M, 0.1);
            }else if ( strcmp(comparisonMethod,"CC")==0 ){
                //printf("Cross correlation\n");
                score=crossCorrelationDistance(I1, I2,  nxyz, M, 0.1);
            }else if ( strcmp(comparisonMethod,"SSIM")==0 ){
                //printf(">>>>SSIM\n");
                //double minMaskValue=0.1;
                score=SSIM(I1, I2,  nxyz, M, 0.1);
            }else if ( strcmp(comparisonMethod,"PSNR")==0 ){
                 //printf("PSNR\n");
                //double minMaskValue=0.1;
                score=PSNR(I1, I2,  nxyz, M, 0.1);
            }else if ( strcmp(comparisonMethod,"MI")==0 ){
                //printf("Mutual Information\n");
                unsigned long int numBins = 30;
                double minMaskValue=0.1;
                score = NormalisedMutualInformation(I1, I2, nxyz, numBins, M, minMaskValue);
            }else if ( strcmp(comparisonMethod,"SCI")==0 ){
                //printf("Mutual Information\n");
                //unsigned long int numBins = 30;
                double minMaskValue=0.1;
                //score = emprove_SCI_distance(I1, I2, nx, ny, nz, sigmaStr, M, minMaskValue);
                score=SCI((float * )I1, (float *) I2, nx, ny, nz, sigma, (float *)M, minMaskValue);
                //std::cerr<<"\n**********\n";
            }
            if (score != score){
                score = 0;
            }

            //std::cerr<<"comparisonMethod="<<comparisonMethod<<", preprocessingMethod="<<preprocessingMethod<<",  sigma="<<sigma<<"  score="<<score<<"\n";

            delete [] I1;
            delete [] I2;
            delete [] M;
        }
	}
    return Py_BuildValue("d",score);
}


// ****************************
// infoMrc
static PyObject * emprove_sizeMRC(PyObject *self, PyObject *args){

    char *map3D_in;

    if (!PyArg_ParseTuple(args, "s", &map3D_in))
        return NULL;
 
      MRCHeader mapHeader;
      readHeaderMrc(map3D_in, mapHeader);
      unsigned long int nx=mapHeader.nx;
      unsigned long int ny=mapHeader.ny;
      unsigned long int nz=mapHeader.nz;


    //gets a map and projects it according to coordinates, from center
//    char *msg="ZProjection!\n";
//    printf("Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);
    return Py_BuildValue("kkk",nx,ny,nz);
}

// ****************************
// createMRC_emptyImage
static PyObject * emprove_writeEmptyMRC(PyObject *self, PyObject *args){
    char *map3D;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int nz;
//    double aPix = 1.0;
    if (!PyArg_ParseTuple(args, "skkk", &map3D,&nx,&ny,&nz))
        return NULL;
    MRCHeader header(MRCHeader(nx, ny, nz));
    writeEmptyMrcImage(map3D, header);
    return Py_BuildValue("d",0);
}

// ****************************
// ReadMrcSlice
static PyObject * emprove_ReadMrcSlice(PyObject *self, PyObject *args){
//int readMrcSlice(const char * filenameMRC, T * I, MRCHeader & header, long int sliceNumber){
    char *map3D_in;
    unsigned long int sliceNumber;
    if (!PyArg_ParseTuple(args, "sk", &map3D_in, &sliceNumber)){
        printf("ERROR: reading parameters");
        return NULL;
    }

    MRCHeader header;
    readHeaderMrc(map3D_in, header);

    if ( sliceNumber >= (unsigned long int)header.nz){
        printf("ERROR: out of range sliceNumber\n");
        return NULL;
    }

    unsigned long int nxy = header.nx * header.ny;
    float * dataI = new float [nxy];
    //readMrcImage(map3D_in, dataI, header);
    readMrcSlice(map3D_in, dataI, header, sliceNumber);

    PyObject *lst = PyList_New(nxy);
    if (!lst){
        printf("ERROR: memory error, creating lst\n");
        return NULL;
    }
    for (unsigned long i = 0; i < nxy; i++) {
        PyObject *num = PyFloat_FromDouble(dataI[i]);
        if (!num) {
            Py_DECREF(lst);
            printf("ERROR: memory error, creating Py_DECREF\n");
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] dataI;
    return lst;
}



// ****************************
// replaceMrcSlice
static PyObject * emprove_replaceMrcSlice(PyObject *self, PyObject *args){
    PyObject * dataIn;
    char *map3D_out;
    unsigned long int nx;
    unsigned long int ny;
    unsigned long int sliceNumber;
    double aPix;

    if (!PyArg_ParseTuple(args, "Oskkk|d", &dataIn, &map3D_out, &nx, &ny, &sliceNumber, &aPix))
        return NULL;
    unsigned long int nxy = nx * ny;
	if (PyList_Check(dataIn)) {
        if ((Py_ssize_t)nxy == PyList_Size(dataIn)){
            float * dataI = new float [nxy];
            for(Py_ssize_t i = 0; i < (Py_ssize_t)nxy; i++) {
                PyObject *value = PyList_GetItem(dataIn, i);
                dataI[i] = PyFloat_AsDouble(value);
            }
            MRCHeader mapHeader;
            readHeaderMrc(map3D_out, mapHeader);
            //printf("nx =%f  ny=%f  mapHeader.nx=%f  mapHeader.ny=%f\n",nx,ny,mapHeader.nx,mapHeader.ny);
            if ((unsigned long int)mapHeader.nx == nx  && (unsigned long int)mapHeader.ny==ny){
                replaceMrcSlice(map3D_out, dataI, mapHeader, sliceNumber);
            //    writeMrcImage("foo.mrc",  dataI, nx,  ny,  1, 1);
            }
            delete [] dataI;
        }
	}
    return Py_BuildValue("d",0);
}


// ****************************
// ReadMRC
static PyObject * emprove_ReadMRC(PyObject *self, PyObject *args){

    char *map3D_in;
    if (!PyArg_ParseTuple(args, "s", &map3D_in))
        return NULL;

    MRCHeader header;
    readHeaderMrc(map3D_in, header);
    unsigned long int nxyz = header.nx * header.ny * header.nz;
    float * dataI = new float [nxyz];
    readMrcImage(map3D_in, dataI, header);

    PyObject *lst = PyList_New(nxyz);
    if (!lst)
        return NULL;
    for (unsigned long i = 0; i < nxyz; i++) {
        PyObject *num = PyFloat_FromDouble(dataI[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] dataI;
    return lst;
}





// ****************************
//ZPROJECT
static PyObject * emprove_ZProject(PyObject *self, PyObject *args){

    char *map3D_in;
    char *map3D_out;
    char *mode;

    if (!PyArg_ParseTuple(args, "sss", &map3D_in, &map3D_out, &mode))
        return NULL;
 
      MRCHeader mapHeader;
      readHeaderMrc(map3D_in, mapHeader);
      unsigned long int nx=mapHeader.nx;
      unsigned long int ny=mapHeader.ny;
      unsigned long int nz=mapHeader.nz;
      unsigned long int nxyz=nx*ny*nz;
      unsigned long int nxy=nx*ny;

      float * referenceMapI = new float[nxyz];
      float * projectionI = new float[nxy];
      readMrcImage(map3D_in, referenceMapI, mapHeader);

//template<typename T,typename U>
//void ZProject(T* inI, U* outI, unsigned long int nx, unsigned long int ny, unsigned long int nz){

    for (unsigned long int jj=0; jj<ny; jj++){
       unsigned long int jjnx = jj*nx;
       for (unsigned long int ii=0; ii<nx; ii++){
            std::vector<double> values;
            double meanValue = 0;
            for (unsigned long int kk=0; kk<nz; kk++){
                meanValue+=referenceMapI[ii + jjnx + kk*nxy];
                values.push_back( referenceMapI[ii + jjnx + kk*nxy] );
            }
            std::vector<double> valuesSorted = bubbleSortAscendingValues(values);
            if (strcmp(mode,"min")==0){
                projectionI[ii+jj*nx]=valuesSorted[0];
            } else if (strcmp(mode,"max")==0 ){
                projectionI[ii+jj*nx]=valuesSorted[valuesSorted.size()-1];
            }else if (strcmp(mode,"mean")==0 ){
                projectionI[ii+jj*nx]=meanValue/values.size();
            }else if (strcmp(mode,"median")==0 ){
                int meanIdx=floor(valuesSorted.size()/2.0);
                projectionI[ii+jj*nx]=valuesSorted[meanIdx];
            }
       }
    }
    writeMrcImage(map3D_out, projectionI, nx, ny, 1);
    delete [] referenceMapI;
    delete [] projectionI;

    //gets a map and projects it according to coordinates, from center
    //char *msg="ZProjection!\n";
//    printf("Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);
    return Py_BuildValue("d",0);
}


// ***********************************************
// compute Ctf Centered Image (to multiply in CTF)
/*
static PyObject * emprove_CtfCenteredImage(PyObject *self, PyObject *args){
    unsigned long int nx; //k is for unsigned long, https://docs.python.org/3/c-api/arg.html
    unsigned long int ny;
    double angpix;

    double SphericalAberration;
    double voltage;
    double DefocusAngle;
    double DefocusU;
    double DefocusV;
    double AmplitudeContrast;
    double Bfac;
    double phase_shift;
    


    //printf("Pre Map=%s phi=%d theta=%d psi=%d centerX=%d centerY=%d\n", map3D, phi, theta, psi, centerX, centerY);

    if (!PyArg_ParseTuple(args, "kk", &nx, &ny))
        return NULL;
    unsigned long int nxy=nx*ny;

    char *msg="CtfCenteredImage!!!!\n";
    //PyFloat_FromDouble

//int computeCtfCenteredImage(T * ctfCenteredImage2D, CTFParameters & ctf_parameters, unsigned long int nx, unsigned long int ny, double angpix){


    float * outImage=new float [nxy];
    for (unsigned long int ii=0;ii<nxy;ii++){
        outImage[ii]=ii;
    }

    PyObject *lst = PyList_New(nxy);
    if (!lst)
        return NULL;
    for (unsigned long i = 0; i < nxy; i++) {
        PyObject *num = PyFloat_FromDouble(outImage[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] outImage;
    return lst;
}
*/

// **************************
// compute ctf centered image
static PyObject * emprove_CtfCenteredImage(PyObject *self, PyObject *args){
    unsigned long int nx;
    unsigned long int ny;
    double angpix;
    double SphericalAberration;
    double voltage;
    double DefocusAngle;
    double DefocusU; 
    double DefocusV;
    double AmplitudeContrast; 
    double Bfac;
    double phase_shift;
    double _scale=1;

    if (!PyArg_ParseTuple(args, "kkddddddddd|d", &nx, &ny, &angpix, &SphericalAberration, &voltage, &DefocusAngle, &DefocusU, &DefocusV, &AmplitudeContrast, &Bfac, &phase_shift, & _scale))
            return NULL;

    unsigned long int nxy = nx * ny;
    float * dataI = new float [nxy];
    CTFParameters ctf_parameters (SphericalAberration, voltage, DefocusAngle, DefocusU, DefocusV, AmplitudeContrast, Bfac, phase_shift, _scale);
    computeCtfCenteredImage(dataI, ctf_parameters,  nx,  ny,  angpix);

    PyObject *lst = PyList_New(nxy);
    if (!lst)
        return NULL;
    for (unsigned long i = 0; i < nxy; i++) {
        PyObject *num = PyFloat_FromDouble(dataI[i]);
        if (!num) {
            Py_DECREF(lst);
            return NULL;
        }
        PyList_SET_ITEM(lst, i, num);   // reference to num stolen
    }
    delete [] dataI;
    return lst;
}

// ****************************
// replaceMrcHeader
static PyObject * emprove_replaceMrcHeader(PyObject *self, PyObject *args){

    char *map3D_in;
    char *map3D_out;
    if (!PyArg_ParseTuple(args, "ss", &map3D_in, &map3D_out))
        return NULL;
    MRCHeader header;
    readHeaderMrc(map3D_in, header);
    copyHeaderMrcImage(header, map3D_out);
    return Py_BuildValue("d",0);


}

//static char emprove_docs[]="emprove procedure to project images\n";
static PyMethodDef emprove_funcs[] = {
    {"project",  (PyCFunction)emprove_project, METH_VARARGS, "program that writes project"},
    {"projectMap",  (PyCFunction)emprove_projectMap, METH_VARARGS, "program that writes project"},
    {"backprojectParticles",  (PyCFunction)emprove_backprojectParticles, METH_VARARGS, "program that writes backprojections"},
    {"projectMask",  (PyCFunction)emprove_projectMask, METH_VARARGS, "program that writes project"},
    {"WriteMRC",  (PyCFunction)emprove_WriteMRC, METH_VARARGS, "write mrc File"},
    {"ReadMRC",  (PyCFunction)emprove_ReadMRC, METH_VARARGS, "read mrc File"},
    {"WriteEmptyMRC",  (PyCFunction)emprove_writeEmptyMRC, METH_VARARGS, "read mrc File"},
    {"sizeMRC",  (PyCFunction)emprove_sizeMRC, METH_VARARGS, "size mrc File"},
    {"ReadMrcSlice",  (PyCFunction)emprove_ReadMrcSlice, METH_VARARGS, "read mrc slice from File"},
    {"ReplaceMrcSlice",  (PyCFunction)emprove_replaceMrcSlice, METH_VARARGS, "write mrc slice on File"},
    {"CtfCenteredImage",  (PyCFunction)emprove_CtfCenteredImage, METH_VARARGS, "produce CtfCenteredImage"},
    {"ZProject",  (PyCFunction)emprove_ZProject, METH_VARARGS, "zprojection images (min,max,mean,median)"},
    {"MaskedImageComparison",  (PyCFunction)emprove_MaskedImageComparison, METH_VARARGS, "masked comparison between images (CC,)"},
    {"EqualizedParticlesRank",  (PyCFunction)emprove_EqualizedParticlesRank, METH_VARARGS, "equalized rank of particles"},
    {"GetEulerClassGroup",  (PyCFunction)emprove_GetEulerClassGroup, METH_VARARGS, "Get Euler Class Group"},
    {"maskedNormalizeInfo",  (PyCFunction)emprove_maskedNormalizeInfo, METH_VARARGS, "(mean1, max1-min1), (mean2, max2-min2)"},
    {"derivative",  (PyCFunction)emprove_derivative, METH_VARARGS, "image, blur, nx, ny, derivative=0..1..2.."},
    {"derivative1D",  (PyCFunction)emprove_derivative1D, METH_VARARGS, "image, blur, nx, ny, nz, derivativeOrder=0..1..2.., direction==0..1..2.. (direction 0 is X, direction 1 is Y, etc)"},
    {"SCI_distance",  (PyCFunction)emprove_SCI_distance, METH_VARARGS, "image1, image2, mask, sigma_blur, nx, ny, nz"},

    {"SDIM",  (PyCFunction)emprove_SDIM, METH_VARARGS, "image1, image2, mask, sigma_blur, nx, ny, nz"},
    {"rotAverage",  (PyCFunction)emprove_rotAverage, METH_VARARGS, "image, nx, ny, nz"},

    {"createHalfMaps",  (PyCFunction)emprove_createHalfMaps, METH_VARARGS, "image, nx, ny, nz"},
    {"amplitudeNormalization",  (PyCFunction)emprove_amplitudeNormalization, METH_VARARGS, "AmplitudesImage, amplitudesToReplace, nx,ny,xz"},
    {"maskedMapContrast",  (PyCFunction)emprove_maskedMapContrast, METH_VARARGS, "image, mask"},
    {"AutomaticThresholdMap",  (PyCFunction)emprove_AutomaticThresholdMap, METH_VARARGS, "image, map"},
    {"replaceMrcHeader",  (PyCFunction)emprove_replaceMrcHeader, METH_VARARGS, "imageSource, imageDst"},


    {NULL, NULL, 0, NULL}
};


static struct PyModuleDef emprove_module = {
    PyModuleDef_HEAD_INIT,
    "emprove_core",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    emprove_funcs
};

PyMODINIT_FUNC PyInit_emprove_core(void){
    return PyModule_Create(&emprove_module);
}

