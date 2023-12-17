from emprove import starHandler
from os import PathLike
from os import path
import emprove_core
import pandas as pd
import numpy as np

#get_pixel_spacing
def get_MRC_map_pixel_spacing(filename):
    import struct
    with open(filename, 'rb') as f:
        # Read the number of columns, rows, and sections (bytes 0-11)
        nx = struct.unpack('i', f.read(4))[0]
        ny = struct.unpack('i', f.read(4))[0]
        nz = struct.unpack('i', f.read(4))[0]

        # Seek to the cell dimensions (bytes 40-51)
        f.seek(40)
        x_dim = struct.unpack('f', f.read(4))[0]
        y_dim = struct.unpack('f', f.read(4))[0]
        z_dim = struct.unpack('f', f.read(4))[0]

    # Calculate the actual pixel spacing by dividing dimensions by the number of voxels
    apix_x = x_dim / nx
    apix_y = y_dim / ny
    apix_z = z_dim / nz
    return apix_x, apix_y, apix_z

#############################
###  create ctfStack
def ctfStack(particlesStarFile: PathLike, outputStackBasename, saveOriginal=False):
    imageNameTag='_rlnImageName'
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    
#    if '_rlnPhaseShift' in imageNames:
#    	print ('we have _rlnPhaseShift')
#    else:
#    	print ('we DO NOT have _rlnPhaseShift')

    #CTFParameters ctf_parameters(rlnSphAberration[ii], rlnVoltage[ii], rlnDefocusAngle[ii], rlnDefocusU[ii], rlnDefocusV[ii], rlnAmplitudeContrast[ii], 0, 0);
    version=starHandler.infoStarFile(particlesStarFile)[2]
    if version=="relion_v31":
        
        if '_rlnPhaseShift' in imageNames:
        	parametersFULL = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnOpticsGroup','_rlnCtfBfactor', '_rlnPhaseShift'])
        else:
        	parametersFULL = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnOpticsGroup','_rlnCtfBfactor'])
    		
        idx=[x for x in range(0, len(parametersFULL))]
        parametersFULL['idx']=idx
        parametersDataOptics = starHandler.dataOptics(particlesStarFile)[['_rlnImagePixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnOpticsGroup']]
        ctfParameters =  pd.merge(parametersFULL, parametersDataOptics,  on=['_rlnOpticsGroup']).sort_values(['idx'])
        ctfParameters=ctfParameters.drop(['_rlnOpticsGroup'],axis=1).reindex()
        ctfParameters=ctfParameters.set_index('idx')
        ctfParameters.rename(columns={'_rlnImagePixelSize':'_rlnDetectorPixelSize'},inplace=True)
    else:
        if '_rlnPhaseShift' in imageNames:
        	ctfParameters = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnDetectorPixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnCtfBfactor', '_rlnPhaseShift'])
        else:
        	ctfParameters = starHandler.readColumns(particlesStarFile, ['_rlnImageName','_rlnDefocusU','_rlnDefocusV','_rlnDefocusAngle','_rlnDetectorPixelSize','_rlnVoltage','_rlnAmplitudeContrast','_rlnSphericalAberration','_rlnCtfBfactor'])

    if not '_rlnPhaseShift' in imageNames:
    	ctfParameters['_rlnPhaseShift']=np.zeros(len(ctfParameters))

    tmpLine= (imageNames[imageNameTag][0])
    stackName=tmpLine[tmpLine.find('@')+1:]
    sizeI=emprove_core.sizeMRC(stackName)
    emprove_core.WriteEmptyMRC(outputStackBasename+'.mrcs',sizeI[0],sizeI[1],len(imageNames[imageNameTag]))
    nx=sizeI[0]
    ny=sizeI[1]
    outImageNames = []
    numIterations=len(imageNames[imageNameTag])
    for ii in range(0,len(imageNames[imageNameTag])):
        print('iteration ',ii,' out of ', numIterations,end='\r')
        tmpLine= (imageNames[imageNameTag][ii])
        atPosition=tmpLine.find('@')
        imageNo=int(tmpLine[:atPosition])
        stackName=tmpLine[atPosition+1:]
        #print (ii,' >> ',imageNo,' >> ',imageNames[imageNameTag][ii])
        angpix=ctfParameters.at[ii, '_rlnDetectorPixelSize']
        SphericalAberration=ctfParameters.at[ii, '_rlnSphericalAberration']
        voltage=ctfParameters.at[ii, '_rlnVoltage']
        DefocusAngle=ctfParameters.at[ii, '_rlnDefocusAngle']
        DefocusU=ctfParameters.at[ii, '_rlnDefocusU']
        DefocusV=ctfParameters.at[ii, '_rlnDefocusV']
        AmplitudeContrast=ctfParameters.at[ii, '_rlnAmplitudeContrast']
        Bfac=ctfParameters.at[ii, '_rlnCtfBfactor']
        phase_shift=ctfParameters.at[ii, '_rlnPhaseShift']
        ctfI=emprove_core.CtfCenteredImage(nx,ny,angpix,SphericalAberration,voltage,DefocusAngle,DefocusU,DefocusV,AmplitudeContrast,Bfac,phase_shift)
        ctfI=np.array(ctfI).reshape((sizeI[1], sizeI[0]))
        mapI=np.array(emprove_core.ReadMrcSlice(stackName,imageNo-1))
        mapI=mapI.reshape((sizeI[1], sizeI[0]))
        mapFFT=np.fft.fftshift(np.fft.fft2(mapI))
        mapIFFT=np.real(np.fft.ifft2(np.fft.ifftshift(mapFFT*ctfI))).flatten()
        emprove_core.ReplaceMrcSlice(mapIFFT.tolist(),outputStackBasename+'.mrcs',sizeI[0],sizeI[1],ii)
        outImageNames.append(str(str(ii+1).zfill(7)+'@'+outputStackBasename+'.mrcs'))

#    for ii in range(0,len(outImageNames)):
#        print (outImageNames[ii])

    originalStarDataframe=starHandler.readStar(particlesStarFile)
    if saveOriginal==True:
    	print ('keeping the original')
    	#originalStarDataframe.rename(columns = {'_rlnImageName':'_original_rlnImageName'}, inplace = True)
    	originalStarDataframe['_original_rlnImageName']=imageNames
    originalStarDataframe['_rlnImageName']=outImageNames
    print (originalStarDataframe)
    starHandler.writeDataframeToStar(particlesStarFile, outputStackBasename+'.star', originalStarDataframe)


    
    #Phi=starHandler.readColumns(starFileIn, ['_rlnAngleRot'])

#############################
###  Read MRC
def readMRC(mrcFile: PathLike):
    print("read MRC file ", mrcFile)
    if (not path.exists(mrcFile)):
        print ('ERROR: file ',mrcFile,' does not exists')
        return None
    if not mrcFile.endswith('.mrc') or not mrcFile.endswith('.mrcs') or not mrcFile.endswith('.st') or not not mrcFile.endswith('.rawst'):
        print ('ERROR: file ',mrcFile,' does not have a recognized extension. If you sure it is a mrc file, just rename it as .mrc')
        return None

    #if (particlesStarFile.endswith(suffix))


#############################
###  maps Difference 
def mapsDifference(mrcFileList,mapout):
    for ii in range (0, len(mrcFileList)):
        mrcImage=mrcFileList[ii]
        if (not path.exists(mrcImage)):
            print ('ERROR: file ',mrcImage,' does not exists')
            return None
        #if not  mrcImage.endswith('.mrc') or not mrcImage.endswith('.mrcs') or not mrcImage.endswith('.st') or not not mrcImage.endswith('.rawst'):
        #    print ('ERROR: file ',mrcImage,' does not have a recognized extension. If you sure it is a mrc file, just rename it as .mrc')
        #    return None
    AvgMap=np.array(emprove_core.ReadMRC(mrcFileList[0]))
    sizeMap=emprove_core.sizeMRC(mrcFileList[0])
    for ii in range (1, len(mrcFileList)):
        AvgMap=AvgMap+np.array(emprove_core.ReadMRC(mrcFileList[ii]))
    AvgMap=AvgMap/len(mrcFileList)
    emprove_core.WriteMRC(AvgMap.tolist(), 'tmpAvg.mrc' ,sizeMap[0],sizeMap[1],sizeMap[2],1)
    emprove_core.replaceMrcHeader(mrcFileList[0],'tmpAvg.mrc')

    #amplitude comp
    AvgMap=np.reshape(AvgMap,sizeMap)
    #I_fft=np.fft.fftn(np.reshape(AvgMap,sizeMap))
    #I_abs_fft=np.abs(I_fft)+0.0000001
    diffMap=np.zeros(sizeMap)
    for ii in range (0, len(mrcFileList)):
        mapTmp=np.reshape(np.array(emprove_core.ReadMRC(mrcFileList[ii])),sizeMap)
        #tmpI_fft=np.fft.fftn(np.reshape(mapTmp,sizeMap))
        #tmpI=np.fft.ifftn(I_fft*(I_abs_fft/abs(tmpI_fft))).real
        #diffMap=diffMap+np.square(tmpI-AvgMap)
        diffMap=diffMap+np.square(mapTmp-AvgMap)

    #diffMap=diffMap/(len(mrcFileList)+np.square(AvgMap))

    #I_out=np.fft.ifftn(I_fft*(RI_abs_fft/I_abs_fft)).real.flatten().tolist()
    emprove_core.WriteMRC(diffMap.flatten().tolist(), mapout ,sizeMap[0],sizeMap[1],sizeMap[2],1)
    emprove_core.replaceMrcHeader(mrcFileList[0],mapout)






    #if (particlesStarFile.endswith(suffix))
