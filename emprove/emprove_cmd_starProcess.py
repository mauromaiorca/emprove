#!/usr/local/bin/python3.9

# ! /usr/bin/python3.8

import matplotlib.pyplot as plt
import numpy as np

#import emprove
#from emprove import starHandler
#from emprove import starDisplay

#from emprove import assessParticles


from os import PathLike
import pandas as pd
from emprove import starHandler
import emprove_core
import numpy as np
#import metrics
#import tensorflow as tf
#import matplotlib.pyplot as plt
from skimage.metrics import structural_similarity as ssim
from skimage.metrics import peak_signal_noise_ratio
import operator
import timeit

import multiprocessing

def transformCtfImage(I,nx,ny,angpix,Voltage,DefocusU, DefocusV, DefocusAngle, SphericalAberration, CtfBfactor, PhaseShift, AmplitudeContrast, DetectorPixelSize):
    #print (Voltage,'   angpix=',angpix)
    ctfI=emprove_core.CtfCenteredImage(nx,ny,angpix,SphericalAberration,Voltage,DefocusAngle,DefocusU,DefocusV,AmplitudeContrast,CtfBfactor,PhaseShift)
    if (not len(ctfI) == nx*ny) or  (not len(I) == nx*ny):
	    #print ('\n\n size error  nx=',nx,'   ny=',ny,'   lenCtf=',len(ctfI) ,'   lenI=',len(I) ,'\n')
	    return np.zeros(nx*ny).tolist()
    else:
	    mapI=np.array(I).reshape((ny, nx))
	    ctfI=np.array(ctfI).reshape((ny, nx))
	    mapFFT=np.fft.fftshift(np.fft.fft2(mapI))
	    mapIFFT=np.real(np.fft.ifft2(np.fft.ifftshift(mapFFT*ctfI))).flatten().tolist()
	    #plt.imshow(np.array(mapIFFT).reshape((ny, nx)))
	    #plt.show()
	    return mapIFFT
 
def scoreBlockParticles (particleIdxStart, particleIdxEnd, listScoresTags, mapI, maskI, sizeMap, angpix, pandasStarHeader, subsetCtfParameters, procnum, return_dict):
#      print ('ecchimi:')
#      print (subsetCtfParameters)
#      print ('ecchimi detail:')
      tmpScore=np.zeros([int(particleIdxEnd)-int(particleIdxStart),len(listScoresTags)+1])
      #print ('idx=',particleIdxStart, ' ', particleIdxEnd)
      counter=-1
      for ii in range(int(particleIdxStart),int(particleIdxEnd)):
        counter+=1
        if procnum==0:
            percentage=100*ii/float(particleIdxEnd)
            print('percentage= ','{:10.4f}'.format(percentage),'%  ',end='\r')
            #print('percentage= ','{:10.4f}'.format(percentage),'%  ')
        tmpLine=(pandasStarHeader['_rlnImageName'][ii])
        atPosition=tmpLine.find('@')
        imageNo=int(tmpLine[:atPosition])
        stackName=tmpLine[atPosition+1:]
        phi=pandasStarHeader.at[ii, '_rlnAngleRot']
        theta=pandasStarHeader.at[ii, '_rlnAngleTilt']
        psi=pandasStarHeader.at[ii, '_rlnAnglePsi']
        tx=pandasStarHeader.at[ii, '_rlnOriginX']
        ty=pandasStarHeader.at[ii, '_rlnOriginY']
        #outImageNames.append(str(str(ii+1).zfill(7)+'@'+stackName+'.mrcs'))
        I=emprove_core.ReadMrcSlice(stackName,imageNo-1)
        #nx=sizeMap[0]
        #ny=sizeMap[1]
        #if not len(I) == nx*ny:
        #    print ('\n\n size error  nx=',nx,'   ny=',ny,'   lenI=',len(I) )
        #    print ('   stackName=',stackName,'   imageNo=',imageNo-1,'\n')
        #    for iii in range(0,5):
        #    	I=emprove_core.ReadMrcSlice(stackName,imageNo-1)
        #    	print ('attempt ',iii,',  size=',len(I),'   (expected=',nx*ny,')')
        #    print('\n\n')
        #    I=np.zeros(nx*ny).tolist()
        RI=emprove_core.projectMap(mapI,sizeMap[0],sizeMap[1],sizeMap[2],phi,theta,psi,tx,ty,0)
        MI=emprove_core.projectMask(maskI,sizeMap[0],sizeMap[1],sizeMap[2],phi,theta,psi,tx,ty,0,0.5)

        if not subsetCtfParameters.empty:
            #print ("doing CTF...")
            Voltage=subsetCtfParameters.iloc[counter]['_rlnVoltage']
            DefocusU=subsetCtfParameters.iloc[counter]['_rlnDefocusU']
            DefocusV=subsetCtfParameters.iloc[counter]['_rlnDefocusV']
            DefocusAngle=subsetCtfParameters.iloc[counter]['_rlnDefocusAngle']
            SphericalAberration=subsetCtfParameters.iloc[counter]['_rlnSphericalAberration']
            CtfBfactor=subsetCtfParameters.iloc[counter]['_rlnCtfBfactor']
            PhaseShift=subsetCtfParameters.iloc[counter]['_rlnPhaseShift']
            AmplitudeContrast=subsetCtfParameters.iloc[counter]['_rlnAmplitudeContrast']
            DetectorPixelSize=subsetCtfParameters.iloc[counter]['_rlnDetectorPixelSize']
            I=transformCtfImage(I,sizeMap[0],sizeMap[1], angpix, Voltage, DefocusU, DefocusV, DefocusAngle, SphericalAberration, CtfBfactor, PhaseShift, AmplitudeContrast, DetectorPixelSize)

        tmpScore[counter][0]=ii
        for kk in range(0,len(listScoresTags)):
            comparisonMethod = listScoresTags[kk].split('_')[2] if len(listScoresTags[kk].split('_'))>2 else "CC"
            preprocessingMethod = listScoresTags[kk].split('_')[3] if len(listScoresTags[kk].split('_'))>3 else "unprocessed"
            sigmaBlur = str(listScoresTags[kk].split('_')[4]) if len(listScoresTags[kk].split('_'))>4 else "1"
            #print ( 'comparing methods=',comparisonMethod,' ', preprocessingMethod)
            #print ('comparisonMethod=',comparisonMethod)
            if comparisonMethod=="CC" or comparisonMethod=="MI": 
                #I_fft=np.fft.fft(np.reshape(I,[sizeMap[1],sizeMap[0]]))
                #RI_fft=np.fft.fft(np.reshape(RI,[sizeMap[1],sizeMap[0]]))
                #I_abs_fft=np.abs(I_fft)+0.0000001
                #RI_abs_fft=np.abs(RI_fft)+0.0000001
                #ampAvg=0.5*(I_abs_fft+RI_abs_fft)
                #I_out=np.fft.ifft(I_fft*(ampAvg/I_abs_fft)).real.flatten().tolist()
                #RI_out=np.fft.ifft(RI_fft*(ampAvg/RI_abs_fft)).real.flatten().tolist()
                #tmpScore.append(emprove_core.MaskedImageComparison(I_out, RI_out, MI, sizeMap[0], sizeMap[1], 1,comparisonMethod,preprocessingMethod,sigmaBlur))
                tmpScore[counter][kk+1]=(emprove_core.MaskedImageComparison(I, RI, MI, sizeMap[0], sizeMap[1], 1,comparisonMethod,preprocessingMethod,sigmaBlur))
            elif comparisonMethod=="SSIM":
                #if preprocessingMethod.lower()=='blur':
                #    RI=emprove_core.derivative(I,float(sigmaBlur),sizeMap[0], sizeMap[1],0)
                #mean1,diff1,mean2,diff2= emprove_core.maskedNormalizeInfo(I,RI,MI,sizeMap[1]*sizeMap[0] )
                #img1 = ((np.array(I,dtype=np.float32)-mean1)/diff1).reshape(sizeMap[1],sizeMap[0])
                #img2 = ((np.array(RI,dtype=np.float32)-mean2)/diff2).reshape(sizeMap[1],sizeMap[0])
                #scoreSSIM=ssim(img1,img2)
                scoreSSIM=emprove_core.MaskedImageComparison(RI, I, MI, sizeMap[0], sizeMap[1], 1,comparisonMethod,preprocessingMethod,sigmaBlur)
                tmpScore[counter][kk+1]=(scoreSSIM)
                #plt.imshow(img2)
                #plt.show()
            elif comparisonMethod=="PSNR":
                #img1 = ((np.array(I,dtype=np.float32)-mean1)/diff1).reshape(sizeMap[1],sizeMap[0])
                #img_GT = ((np.array(RI,dtype=np.float32)-mean2)/diff2).reshape(sizeMap[1],sizeMap[0])
                #scorePSNR=peak_signal_noise_ratio(img_GT, img1)
                scorePSNR=emprove_core.MaskedImageComparison(RI, I, MI, sizeMap[0], sizeMap[1], 1,comparisonMethod,preprocessingMethod,sigmaBlur)
                tmpScore[counter][kk+1]=(scorePSNR)
            elif comparisonMethod=="SCI":
                I_fft=np.fft.fftn(np.reshape(I,[sizeMap[1],sizeMap[0]]))
                RI_fft=np.fft.fftn(np.reshape(RI,[sizeMap[1],sizeMap[0]]))
                I_abs_fft=np.abs(I_fft)+0.0000001
                RI_abs_fft=np.abs(RI_fft)+0.0000001
                ampAvg=0.5*(I_abs_fft+RI_abs_fft)
                I_out=np.fft.ifftn(I_fft*(ampAvg/I_abs_fft)).real.flatten().tolist()
                RI_out=np.fft.ifftn(RI_fft*(ampAvg/RI_abs_fft)).real.flatten().tolist()

                scoreSCI=emprove_core.MaskedImageComparison(RI_out, I_out, MI, sizeMap[0], sizeMap[1], 1,comparisonMethod,preprocessingMethod,sigmaBlur)
                #scoreSCI=emprove_core.MaskedImageComparison(RI_out, I_out, MI, sizeMap[0], sizeMap[1], 1,comparisonMethod,preprocessingMethod,sigmaBlur)
                tmpScore[counter][kk+1]=(scoreSCI)
      return_dict[procnum] = tmpScore
      #print (tmpScore)
      return tmpScore


def ParticleVsReprojectionScores(particlesStarFile: PathLike, scoredParticlesStarFile: PathLike, referenceMap: PathLike, referenceMask: PathLike, angpix, listScoresTags = "", numProcesses = 1, numViews=[50,200,400], doCTF=False):
    if listScoresTags == "":
        listScoresTags = ['_emprove_SCI__1','_emprove_SCI__1.5']
    print ('listScoresTags=',listScoresTags)

    mapI=emprove_core.ReadMRC(referenceMap)
    sizeMap=emprove_core.sizeMRC(referenceMap)
    maskI=emprove_core.ReadMRC(referenceMask)
#    sizeMask=emprove_core.sizeMRC(referenceMask)



    version=starHandler.infoStarFile(particlesStarFile)[2]
    if version=="relion_v31":
        coordinatesFULL = starHandler.readColumns(particlesStarFile, ['_rlnRandomSubset','_rlnImageName','_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi','_rlnOriginXAngst','_rlnOriginYAngst','_rlnOpticsGroup'])
        idx=[x for x in range(0, len(coordinatesFULL))]
        coordinatesFULL['idx']=idx
        coordinatesDataOptics = starHandler.dataOptics(particlesStarFile)[['_rlnOpticsGroup','_rlnImagePixelSize']]
        coordinates =  pd.merge(coordinatesFULL, coordinatesDataOptics,  on=['_rlnOpticsGroup']).sort_values(['idx'])
        coordinates['_rlnOriginX']=coordinates['_rlnOriginXAngst']/coordinates['_rlnImagePixelSize']
        coordinates['_rlnOriginY']=coordinates['_rlnOriginYAngst']/coordinates['_rlnImagePixelSize']
        coordinates=coordinates.drop(['_rlnOpticsGroup'],axis=1).reindex()
        coordinates=coordinates.set_index('idx')
        coordinates=coordinates.drop(['_rlnOriginXAngst','_rlnOriginYAngst'],axis=1)
    else:
        coordinates = starHandler.readColumns(particlesStarFile, ['_rlnRandomSubset','_rlnImageName','_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi','_rlnOriginX','_rlnOriginY'])

    ctfParameters = pd.DataFrame( [] )
    if doCTF:
        columns=starHandler.header_columns(particlesStarFile)
        print (len(coordinates))
        if not '_rlnPhaseShift' in columns:
            PhaseShift=pd.DataFrame(np.zeros(len(coordinates)))
        else:
            PhaseShift = starHandler.readColumns(particlesStarFile, ['_rlnPhaseShift'])

        #print (columns)

        print ("doing CTF...")
        if version=="relion_v31":
            print('READ parameters')
            parametersFULL = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnOpticsGroup','_rlnCtfBfactor'])
            print('GOT the parameters')
            idx=[x for x in range(0, len(parametersFULL))]
            parametersFULL['idx']=idx
            parametersDataOptics = starHandler.dataOptics(particlesStarFile)[['_rlnImagePixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnOpticsGroup']]
            ctfParameters =  pd.merge(parametersFULL, parametersDataOptics,  on=['_rlnOpticsGroup']).sort_values(['idx'])
            ctfParameters=ctfParameters.drop(['_rlnOpticsGroup'],axis=1).reindex()
            ctfParameters=ctfParameters.set_index('idx')
            ctfParameters.rename(columns={'_rlnImagePixelSize':'_rlnDetectorPixelSize'},inplace=True)
        else:
            (ctfParameters) = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnDetectorPixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnCtfBfactor'])
        ctfParameters['_rlnPhaseShift']=PhaseShift
    numParticles=len(coordinates['_rlnImageName'])
    print('num particles=',numParticles,'  num Processes=',numProcesses)


    if numProcesses > multiprocessing.cpu_count():
        numProcesses = multiprocessing.cpu_count()
    elif numProcesses < 1:
        numProcesses = 1
    if numProcesses > numParticles:
        numProcesses = numParticles
    
    #if numProcesses == 1:
    #manager = multiprocessing.Manager()
    #return_dict = manager.dict()
#OK   subsetCtfParameters=ctfParameters[0:numParticles]
#OK   #print (subsetCtfParameters)
#OK    scoresOut = scoreBlockParticles(0,numParticles,listScoresTags, mapI, maskI, sizeMap, angpix, coordinates, subsetCtfParameters, 0,0)
#OK   df_full_scores = pd.DataFrame(scoresOut)
#OK    df_full_scores = df_full_scores.iloc[: , 1:]
#OK    df_full_scores.columns=listScoresTags
#OK    starHandler.removeColumnsTagsStartingWith(particlesStarFile, scoredParticlesStarFile, "_emprove_")
#OK    starHandler.addDataframeColumns(scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores)
#OK    return 

    blockSize=(np.floor(numParticles/numProcesses))
    idxMatrix=np.zeros([numProcesses,2])
    for mm in range (0, numProcesses):
        idxMatrix[mm][0]=(mm*blockSize)
        idxMatrix[mm][1]=((mm+1)*(blockSize))
        if mm == numProcesses-1:
            idxMatrix[mm][1]=int(numParticles)


    manager = multiprocessing.Manager()
    return_dict = manager.dict()
    jobs=[]
    for ii in range(0, len(idxMatrix) ):
        #print ('**************\n***********\ndebug=',idxMatrix[ii][0]),'   ',int(idxMatrix[ii][1])
        subsetCtfParameters=ctfParameters[int(idxMatrix[ii][0]):int(idxMatrix[ii][1])]
        p=multiprocessing.Process(target=scoreBlockParticles, args=(idxMatrix[ii][0],idxMatrix[ii][1],listScoresTags, mapI, maskI, sizeMap, angpix, coordinates, subsetCtfParameters, ii,return_dict,))
        jobs.append(p)
        p.start()

    for proc in jobs:
        proc.join()

    scoresOut=[]
    for ii in range(0, len(idxMatrix) ):
        scoresOut=scoresOut+return_dict.values()[ii].tolist()
    scoresOut = np.array(scoresOut).reshape( [numParticles,len(listScoresTags)+1] )
    listScoreOut = np.array(sorted(scoresOut, key=operator.itemgetter(0), reverse=False)).reshape( [numParticles,len(listScoresTags)+1] )[:,1:].tolist()
    df_full_scores = pd.DataFrame(data=listScoreOut,columns=listScoresTags)
    starHandler.removeColumnsTagsStartingWith(particlesStarFile, scoredParticlesStarFile, "_emprove_")
    starHandler.addDataframeColumns(scoredParticlesStarFile, scoredParticlesStarFile, listScoresTags, df_full_scores)


def createDiffStack(particlesStarFile: PathLike, outputBasename, referenceMap: PathLike, referenceMask: PathLike, useCTF=False):

    mapI=emprove_core.ReadMRC(referenceMap)
    sizeMap=emprove_core.sizeMRC(referenceMap)
    maskI=emprove_core.ReadMRC(referenceMask)
#    sizeMask=emprove_core.sizeMRC(referenceMask)
    imageNameTag='_rlnImageName'
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    tmpLine= (imageNames[imageNameTag][0])
    stackName=tmpLine[tmpLine.find('@')+1:]
    sizeI=emprove_core.sizeMRC(stackName)
    emprove_core.WriteEmptyMRC(outputBasename+'.mrcs',sizeI[0],sizeI[1],len(imageNames[imageNameTag]))

    version=starHandler.infoStarFile(particlesStarFile)[2]
    if version=="relion_v31":
        coordinatesFULL = starHandler.readColumns(particlesStarFile, ['_rlnRandomSubset','_rlnImageName','_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi','_rlnOriginXAngst','_rlnOriginYAngst','_rlnOpticsGroup'])
        idx=[x for x in range(0, len(coordinatesFULL))]
        coordinatesFULL['idx']=idx
        coordinatesDataOptics = starHandler.dataOptics(particlesStarFile)[['_rlnOpticsGroup','_rlnImagePixelSize']]
        coordinates =  pd.merge(coordinatesFULL, coordinatesDataOptics,  on=['_rlnOpticsGroup']).sort_values(['idx'])
        coordinates['_rlnOriginX']=coordinates['_rlnOriginXAngst']/coordinates['_rlnImagePixelSize']
        coordinates['_rlnOriginY']=coordinates['_rlnOriginYAngst']/coordinates['_rlnImagePixelSize']
        coordinates=coordinates.drop(['_rlnOpticsGroup'],axis=1).reindex()
        coordinates=coordinates.set_index('idx')
        coordinates=coordinates.drop(['_rlnOriginXAngst','_rlnOriginYAngst'],axis=1)
    else:
        coordinates = starHandler.readColumns(particlesStarFile, ['_rlnRandomSubset','_rlnImageName','_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi','_rlnOriginX','_rlnOriginY'])    

    
    nx=sizeI[0]
    ny=sizeI[1]
    outImageNames = []
    numIterations=len(imageNames[imageNameTag])
    for ii in range(0,len(coordinates['_rlnImageName'])):
        print('iteration ',ii,' out of ', numIterations,end='\r')
        
        #print('iteration ',ii,' out of ', numIterations)
        tmpLine= (coordinates['_rlnImageName'][ii])
        atPosition=tmpLine.find('@')
        imageNo=int(tmpLine[:atPosition])
        stackName=tmpLine[atPosition+1:]
        newFileNames=tmpLine[:atPosition]+"@"

        phi=coordinates.at[ii, '_rlnAngleRot']
        theta=coordinates.at[ii, '_rlnAngleTilt']
        psi=coordinates.at[ii, '_rlnAnglePsi']
        tx=coordinates.at[ii, '_rlnOriginX']
        ty=coordinates.at[ii, '_rlnOriginY']
        outImageNames.append(str(tmpLine[:atPosition]+'@'+outputBasename+'.mrcs'))
        I=emprove_core.ReadMrcSlice(stackName,imageNo-1)
        RI=emprove_core.projectMap(mapI,sizeMap[0],sizeMap[1],sizeMap[2],phi,theta,psi,tx,ty,0)
        MI=emprove_core.projectMask(maskI,sizeMap[0],sizeMap[1],sizeMap[2],phi,theta,psi,tx,ty,0,0.5)
        
        
        I_out = np.array(I)
        meanValue1=np.mean( I_out )
        stdValue1=np.std( I_out )
        if stdValue1 > 0:
            I_out = ((I_out-meanValue1)/(stdValue1)).tolist()

        RI_out = np.array(RI)
        meanValue2=np.mean( RI_out )
        stdValue2=np.std( RI_out )
        if stdValue2 > 0:
            RI_out = ((RI_out-meanValue2)/(stdValue2)).tolist()

# emprove_core.WriteMRC( I_out, "I1_preproc1.mrc",  shapeImg1[0], shapeImg1[1],  1, 1)
# emprove_core.WriteMRC( RI_out, "I2_preproc1.mrc", shapeImg1[0], shapeImg1[1],  1, 1)

 #shapeImg1 = shapeImg1 + (1,1,)
        scoreImage=emprove_core.SDIM(RI_out, I_out,nx, ny, 1,str(1), MI)

        Iout=np.array(I)*np.array(MI)*np.array(scoreImage)
        Iout[Iout < 0.0] = 0.0
        Iout=Iout.tolist()
        emprove_core.ReplaceMrcSlice(Iout,outputBasename+'.mrcs',sizeI[0],sizeI[1],ii)

    inputStarFile=starHandler.readStar(particlesStarFile)
    inputStarFile[imageNameTag]=outImageNames
    starHandler.writeDataframeToStar(particlesStarFile, outputBasename+'.star', inputStarFile)


#if __name__ == '__main__':
#    listScoresTags = ['_emprove_SCI__1','_emprove_SCI__1.25','_emprove_SCI__1.5']
#    ParticleVsReprojectionScores('few0.star', 'fewScored.star', 'J868_ctfMult_rec.mrc', 'mask_dilated5_closed.mrc', listScoresTags = listScoresTags, numProcesses =20)




