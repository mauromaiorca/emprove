#!/usr/bin/python3

import argparse
import os.path
from os import PathLike
from emprove import starHandler
from emprove import assessParticles

#from emprove import utils
import emprove_core
import numpy as np
import scipy.stats as stats
import math
import stat

emprove_parser = argparse.ArgumentParser(
    prog="emprove",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

command = emprove_parser.add_subparsers(dest="command")








#############################
###  create ctfStack
def createAverageStack(particlesStarFile, outputStackBasename, numViews):
    import pandas as pd
    from scipy import ndimage

    imageNameTag='_rlnImageName'
    imageNames = starHandler.readColumns(particlesStarFile, [imageNameTag])
    
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

    coordinates = starHandler.readColumns(particlesStarFile, ['_rlnAngleRot','_rlnAngleTilt','_rlnAnglePsi'])
    phiListParticle = coordinates['_rlnAngleRot'].tolist()
    thetaListParticle = coordinates['_rlnAngleTilt'].tolist()
    psiListParticle = coordinates['_rlnAnglePsi'].tolist()
    viewsParticles=emprove_core.GetEulerClassGroup(phiListParticle,thetaListParticle,int(numViews))
    #print("viewsParticles=",viewsParticles)
    tmpLine= (imageNames[imageNameTag][0])
    stackName=tmpLine[tmpLine.find('@')+1:]
    sizeI=emprove_core.sizeMRC(stackName)
    emprove_core.WriteEmptyMRC(outputStackBasename+'.mrcs',sizeI[0],sizeI[1],int(numViews))
    nx=sizeI[0]
    ny=sizeI[1]
    outImageNames = []
    numIterations=len(imageNames[imageNameTag])
    for ii in range(0,len(imageNames[imageNameTag])):
        #print('iteration ',ii,' out of ', numIterations,end='\r')
        tmpLine= (imageNames[imageNameTag][ii])
        atPosition=tmpLine.find('@')
        imageNo=int(tmpLine[:atPosition])
        stackName=tmpLine[atPosition+1:]
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
        mapIFFT=np.real(np.fft.ifft2(np.fft.ifftshift(mapFFT*ctfI)))
        #rotated_img = ndimage.rotate(mapIFFT, psiListParticle[ii]).flatten().tolist()
        mrcAveragedSlice=emprove_core.ReadMrcSlice(outputStackBasename+'.mrcs',int(viewsParticles[ii]))
        mrcAveragedSlice = mapIFFT.flatten().tolist()
        #print ("\n\n",len (mrcAveragedSlice))
        #return()
        #mrcAveragedSlice=[15 for item in mrcAveragedSlice]
        #mrcAveragedSlice=[item + 1 for item in mrcAveragedSlice]
        emprove_core.ReplaceMrcSlice(mrcAveragedSlice,outputStackBasename+'.mrcs',sizeI[0],sizeI[1],int(viewsParticles[ii]))
        #outImageNames.append(str(str(viewsParticles[ii]+1).zfill(7)+'@'+outputStackBasename+'.mrcs'))
    #originalStarDataframe=starHandler.readStar(particlesStarFile)
    #originalStarDataframe['_rlnImageName']=outImageNames
    #starHandler.writeDataframeToStar(particlesStarFile, outputStackBasename+'.star', originalStarDataframe)





##########################################
##########################################
##### emprove_rankParticles 
##########################################
emprove_rankParticles = command.add_parser (
    "rankParticles", description="Normalize Score in star file tag", help='normalize score'
)
emprove_rankParticles.add_argument("--i", required=True, type=str, help="input star file to normalize")
emprove_rankParticles.add_argument("--o", required=False,  type=str, help="output star file")
emprove_rankParticles.add_argument("--tag", required=True,  type=str, help="tag score to normalize")
emprove_rankParticles.add_argument("--tagOut", required=True,  type=str, help="output tag with normalized score")
emprove_rankParticles.add_argument("--avgViews", required=False,  default="", type=str, help="mrc stack with averaged views")
emprove_rankParticles.add_argument("--views", required=False, default=50,  type=int, help="number of views partitioned for normalization")
def rankParticles(args):
    views=args.views
    #tagPrefix=args.prefix
    if (not os.path.isfile(args.i)):
        print('ERROR: file \"',args.i,'\" not existing')
    if args.o is None:
        args.o = args.i
    elif not os.path.exists(os.path.dirname(args.o)):
        try:
            os.makedirs(os.path.dirname(args.o))
        except OSError as e:
            print(f"Error: {e.strerror}")

    columns=starHandler.header_columns(args.i)
    #columns=starHandler.readColumns(args.i, [args.tagOut], sortColumnsNameOrder=True)
    tagOut_view=args.tagOut
    if args.tagOut in columns:
            starHandler.removeColumnsTagsStartingWith(args.i, args.o, args.tagOut)
    if tagOut_view in columns:
            starHandler.removeColumnsTagsStartingWith(args.i, args.o, tagOut_view)


    #print ('args.tag=',args.tag)
    coordinates = starHandler.readColumns(args.i, ['_rlnAngleRot','_rlnAngleTilt','_rlnRandomSubset',args.tag])
    phiListParticle = coordinates['_rlnAngleRot'].tolist()
    thetaListParticle = coordinates['_rlnAngleTilt'].tolist()
    scores = coordinates[args.tag].tolist()
    randomSubset = coordinates['_rlnRandomSubset'].tolist()
    rankedParticles=emprove_core.EqualizedParticlesRank(phiListParticle,thetaListParticle,scores,randomSubset,int(views))
    starHandler.addColumns(args.i, args.o, [args.tagOut], [rankedParticles])
    viewsParticles=emprove_core.GetEulerClassGroup(phiListParticle,thetaListParticle,int(views))
    starHandler.addColumns(args.i, args.o, [tagOut_view], [viewsParticles])


##########################################
##########################################
##### emprove_scoreParticles 
##########################################
emprove_scoreParticles = command.add_parser (
    "scoreParticles", description="score Particles", help='score Particles described in star file'
)
emprove_scoreParticles.add_argument("--i", required=True, type=str, help="input star file with particles to score")
emprove_scoreParticles.add_argument("--mask", required=True, type=str, help="mrc volumetric file with mask")
emprove_scoreParticles.add_argument("--map", required=True, type=str, help="mrc volumetric file with reference map (or first half reference map if map2 is given)")
emprove_scoreParticles.add_argument("--map2", required=False, type=str, help="mrc volumetric file with second half map as reference")
emprove_scoreParticles.add_argument("--apix", required=True, type=float, help="angpix")
emprove_scoreParticles.add_argument("--sigma", required=False, type=float, default=1,  help="sigma blurring in pixels")
emprove_scoreParticles.add_argument("--selectionName", required=False, type=str, default="0",  help="Name for the selection")
emprove_scoreParticles.add_argument("--o", required=False, default=None,  type=str, help="output file with scores")
emprove_scoreParticles.add_argument("--mpi", required=False, default=4,  type=int, help="number of mpi parallel processes")
emprove_scoreParticles.add_argument("--rank", required=False, default="350", type=int, help="rank particles using a given number of Euler viewss")

def scoreParticles(args):
    inputFile=args.i
    outputFile=args.o
    templateMask=args.mask
    templateMap=args.map
    templateMap2=args.map2
    rankViews=args.rank
    #tagPrefix=args.prefix
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit()
    if (not os.path.isfile(templateMask)):
        print('ERROR: mask file \"',templateMask,'\" not existing')
        exit()
    if (not os.path.isfile(templateMap)):
        print('ERROR: map file \"',templateMap,'\" not existing')
        exit()
    if (outputFile==None):
        outputFile=inputFile
    tag='_emprove_SCI__'+"{:.2f}".format(args.sigma)+'_scored_selection_'+args.selectionName
    listScoresTags = [tag]
    print (listScoresTags)

    if args.map2 is None:
        assessParticles.ParticleVsReprojectionScores(inputFile,
            outputFile,
            templateMap,
            templateMask,
            angpix=args.apix,
            numProcesses=args.mpi,
            listScoresTags = listScoresTags,
            doCTF=True)
    else:
        if (not os.path.isfile(templateMap2)):
            print('ERROR: map file \"',templateMap2,'\" not existing')
            exit()
        assessParticles.ParticleVsReprojectionScores_HalfMaps(inputFile,
            outputFile,
            templateMap,
            templateMap2,
            templateMask,
            angpix=args.apix,
            numProcesses=args.mpi,
            listScoresTags = listScoresTags,
            doCTF=True)


    #ranking
    rankingTag=tag+"_norm"+str(rankViews)
    coordinates = starHandler.readColumns(outputFile, ['_rlnAngleRot','_rlnAngleTilt','_rlnRandomSubset',tag])
    phiListParticle = coordinates['_rlnAngleRot'].tolist()
    thetaListParticle = coordinates['_rlnAngleTilt'].tolist()
    scores = coordinates[tag].tolist()
    randomSubset = coordinates['_rlnRandomSubset'].tolist()
    rankedParticles=emprove_core.EqualizedParticlesRank(phiListParticle,thetaListParticle,scores,randomSubset,int(rankViews))
    columns=starHandler.header_columns(outputFile)
    if rankingTag in columns:
            starHandler.removeColumnsTagsStartingWith(outputFile, outputFile, rankingTag)
    starHandler.addColumns(outputFile, outputFile, [rankingTag], [rankedParticles])






##########################################
##########################################
##### scoreClassify 
##########################################
emprove_scoreClassify = command.add_parser (
    "scoreClassify", description="score classify", help='classify particles based on the measured score'
)
emprove_scoreClassify.add_argument("--i", required=True, type=str, help="input star file with scored particles")
emprove_scoreClassify.add_argument("--classLabelOut", required=False, default="_rlnClassNumber", type=str, help="label with the class to be assigned, default=_rlnClassNumber")
emprove_scoreClassify.add_argument("--classScoredLabels", required=False, nargs='+', type=str, help="labels with scores for each reference class")
emprove_scoreClassify.add_argument("--o", required=False, default=None,  type=str, help="output star file with updated classification")
def scoreClassify(args):
    inputFile=args.i
    outputFile=args.o
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit()
    if (outputFile==None):
        outputFile=inputFile
    print ("classLabelOut=",args.classLabelOut)
    print ("classScoredLabels=",args.classScoredLabels)
    class_columns=args.classScoredLabels
    columns=starHandler.header_columns(args.i)
    if class_columns == None:
        print ("WARNING: --classScoredLabels not properly given by the user, using the array labels ending with _class and a classID")
        class_columns = [col for col in columns if col.split('_')[-1].startswith('class') and col.split('_')[-1][5:].isdigit()]
    if class_columns == None:
        print ("ERROR: no suitable classes given by the user, please check the star file")
        exit()
    #print ("class_columns=",class_columns)
    import pandas as pd
    dataDF=starHandler.readColumns(args.i,class_columns)
    #print ("dataDF=",dataDF)
    max_column = dataDF.idxmax(axis=1)
    class_numbers = max_column.str.extract('class(\d+)$')[0]
    #if all the values are not positive, assign class -1, as the program can't decide
    all_non_positive = (dataDF <= 0).all(axis=1)
    class_numbers[all_non_positive] = '-1'
    resultDF = pd.DataFrame({args.classLabelOut: class_numbers})
    #print(resultDF)

    columnsToRemove= [item for item in columns if item in [args.classLabelOut]]
    if not columnsToRemove==[]:
        starHandler.removeColumns(inputFile, outputFile, columnsToRemove)
    starHandler.addDataframeColumns(outputFile, outputFile, columnsToRemove, resultDF)


#    headers_columns=starHandler.header_columns(particlesStarFile)
#    columnsToRemove= [item for item in headers_columns if item in listScoresTags]
#    starHandler.removeColumns(particlesStarFile, scoredParticlesStarFile, columnsToRemove)



##########################################
##########################################
##### emprove_produce_reconstructions_script 
##########################################
emprove_produce_reconstructions_script = command.add_parser (
    "produce_reconstructions_script", description="produce_reconstructions_script", help='Reconstruction Evaluation'
)
emprove_produce_reconstructions_script.add_argument("--i", required=True, type=str, help="input star file with particles to score")
emprove_produce_reconstructions_script.add_argument("--outDir", required=False, default='./', type=str, help="output Dir")
emprove_produce_reconstructions_script.add_argument("--tagRank", required=True, type=str, help="tag for the ranked particles")
emprove_produce_reconstructions_script.add_argument("--mask", required=True, type=str, help="mask to use for evaluation")
#emprove_produce_reconstructions_script.add_argument("--numReconstructions", required=False, type=int, default="5",  help="Number of Reconstructions")
emprove_produce_reconstructions_script.add_argument("--manualParticleSubsets", required=False, type=str, default=None,  help="Comma separated list of manual particles subsets")
emprove_produce_reconstructions_script.add_argument("--scriptName", required=False, default='script_reconstructions.sh', type=str, help="script name as output")
#emprove_produce_reconstructions_script.add_argument("--mask_XYZ_boxsize", required=False, default='FALSE', type=str, help="X,Y,Z,boxsize coordinates for masked locres")
emprove_produce_reconstructions_script.add_argument("--masked_crop", action='store_true', help="check if automatic masked crop")



def produce_reconstructions_script(args):
    if (not os.path.isfile(args.i)):
        print('ERROR: file \"',args.i,'\" not existing')
        exit()
    elif not os.path.exists(args.outDir):
        try:
            os.makedirs(args.outDir)
        except OSError as e:
            print(f"Error: {e.strerror}")

    #print("num_non_null_items=",num_non_null_items)

    if args.manualParticleSubsets is not None:
        print ("Manual Particle Subset Selection")
        numParticlesList = args.manualParticleSubsets.split(',')
 

#    elif args.automaticParticleSubsets:
        #print ("Automatic Particle Subset Selection,")
        #starFile=starHandler.readStar(args.i)
        #num_non_null_items = int(starFile.count()[0])
        #numParticleSubsetsSelected=int(args.automaticParticleSubsets[0])
        #expectedNumberOfPartilces=int(args.automaticParticleSubsets[1])
        #standardDeviation=int(args.automaticParticleSubsets[2])
        #print("numParticleSubsetsSelected=", args.automaticParticleSubsets[0])
        #print("expectedNumberOfPartilces=", args.automaticParticleSubsets[1])
        #print("num_non_null_items=", num_non_null_items)
        #numParticlesList = assessParticles.automaticParticleSubsetSelection(numParticleSubsetsSelected, expectedNumberOfPartilces, num_non_null_items, standardDeviation)


    print(numParticlesList)

    #######################
    #####RECONSTRUCTION 
    reconstruction_command = '''#!/bin/bash
\n\n
##############################
#######  RECONSTRUCTIONS
rec_subset() {
        fileIn=$1
        fileOut_basename=$2
        subset=$3
        if [ -f ${fileOut_basename}_recH${subset}.mrc ]; then
            echo "DOING NOTHING: Reconstructed file ${fileOut_basename}_recH${subset}.mrc exists"
        else
            echo "DOING Reconstruction for file ${fileOut_basename}_recH${subset}.mrc"
            #mpirun --np 28 relion_reconstruct --i ${fileIn} --o  ${fileOut_basename}_recH${subset}.mrc --subset ${subset} --ctf
            relion_reconstruct --i ${fileIn} --o  ${fileOut_basename}_recH${subset}.mrc --subset ${subset} --ctf &
            sleep 40
        fi
}
    \n\n'''
    numParticlesListStr = ','.join(map(str, numParticlesList))
    reconstruction_command+="numParticlesCsv="+numParticlesListStr+"\n"
    reconstruction_command += '''\n
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n'''
    outFile=os.path.join(args.outDir,"norm_"+str(os.path.split(args.outDir)[-1])+"_best${numParticles}")
    reconstruction_command += "    emprove selectBestRanked --i "+args.i+" --o "+outFile+".star --num  ${numParticles} --tag "+ args.tagRank+"\n"
    reconstruction_command += "    rec_subset  "+outFile+".star  "+outFile + " 1  \n"
    reconstruction_command += "    rec_subset  "+outFile+".star  "+outFile + " 2  \n"    
    reconstruction_command += '''\n
    echo ${scorelabelToNormalize}
done
    \n\n'''
    reconstruction_command+="wait\n"

    reconstruction_file_path=os.path.join(args.outDir,args.scriptName)
    with open(reconstruction_file_path, 'w') as f:
        f.write(reconstruction_command)
    os.chmod(reconstruction_file_path, os.stat(reconstruction_file_path).st_mode | stat.S_IXUSR)



    #######################
    #####LOCRES COMMAND 
    locres_command = '''
##############################
#######  LOCRES
locres() {
        file_basename=$1        
        if [ -f ${file_basename}_locres.mrc ]; then
            echo "DOING NOTHING: locres file ${file_basename}_locres.mrc exists"
        else
            echo "DOING locres for file ${file_basename}_locres.mrc"
            mpirun --np 28 relion_postprocess_mpi --i ${file_basename}_recH1.mrc --i2 ${file_basename}_recH2.mrc --o  ${file_basename} --locres --locres_thresholdFSC 0.5
            rm ${file_basename}_locres_fscs.star
            rm ${file_basename}_locres_filtered.mrc
        fi
}

crop_image() {
  local imageIn="$1"
  local imageOut="$2"
  trimvol -x 100,300 -y 100,300 -z 100,300  "$imageIn" "$imageOut"
}

    \n\n'''
    numParticlesListStr = ','.join(map(str, numParticlesList))
    locres_command+="numParticlesCsv="+numParticlesListStr+"\n"
    locres_command += '''
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n'''
    
    outFile=os.path.join(args.outDir,"norm_"+str(os.path.split(args.outDir)[-1])+"_best${numParticles}")
    if args.masked_crop:
        locres_command += "    emprove_utils maskedCrop --mask   "+args.mask + "  --padding 8 --i  " + outFile +"_recH1.mrc --o "+ outFile+"_crop_recH1.mrc\n" 
        locres_command += "    emprove_utils maskedCrop --mask   "+args.mask + "  --padding 8 --i  " + outFile +"_recH2.mrc --o "+ outFile+"_crop_recH2.mrc\n" 
        locres_command += "    locres  "+outFile + "_crop   \n"    
    else:
        locres_command += "    locres  "+outFile + "   \n"    
    locres_command += '''
    echo ${scorelabelToNormalize}
done
    \n\n'''
    reconstruction_file_path=os.path.join(args.outDir,'script_reconstructions.sh')
    with open(reconstruction_file_path, 'a') as f:
        f.write(locres_command)



    #######################
    #####ASSESS LOCRES 
    assess_command = '''
##############################
#######  ASSESS locres
'''
    maskLocresFilename=args.mask
    print ("mask location=",args.mask)
    locresSuffixfix="_locres"
    if  args.masked_crop:
        locresSuffixfix="_locres"
        print("DOING CROP!")
        maskLocresFilename= os.path.join(args.outDir,"mask_crop.mrc")
        assess_command += "emprove_utils maskedCrop --mask " +args.mask + " --padding 8 --i "+ args.mask+"  --o " +maskLocresFilename+"\n"
        locresSuffixfix="_crop_locres"
    else:
        print("Not doing DOING CROP!")

    numParticlesListStr = ','.join(map(str, numParticlesList))
    result_filename=os.path.join(args.outDir,"bestRanked_locres_values.csv")
    assess_command+="numParticlesCsv="+numParticlesListStr+"\n"
    assess_command+="echo numParticles,max,highQuartile,mean,lowQuartile,min >"+result_filename+"\n"
    assess_command += '''
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n'''
    assess_command+='printf "%s," ${numParticles} >> ' + result_filename +'\n'
    outFile=os.path.join(args.outDir,"norm_"+str(os.path.split(args.outDir)[-1])+"_best${numParticles}")
    assess_command += "emprove_app_meanMinMax  "+outFile +locresSuffixfix + ".mrc  "+maskLocresFilename+" >> "+ result_filename +"\n"
    assess_command += '''
    echo ${scorelabelToNormalize}
done
    \n\n'''
    reconstruction_file_path=os.path.join(args.outDir,'script_reconstructions.sh')
    with open(reconstruction_file_path, 'a') as f:
        f.write(assess_command)



    #######################
    #####SELECT OPTIMAL SUBSET 
    select_command = '''
##############################
#######  SELECT OPTIMAL SUBSET
'''
    numParticlesListStr = ','.join(map(str, numParticlesList))
    result_filename=os.path.join(args.outDir,"bestRanked_locres_values.csv")

#    with open(reconstruction_file_path, 'a') as f:
#        f.write(assess_command)




def find_emprove_ranking_tag(elements):
    """
    Search for a column name in the dataframe that starts with '_emprove_' and ends with '_norm' followed by a number.
    """
    import re
    for element in elements:
        if re.match(r'_emprove_.*_norm\d+$', element):
            return element
    return None



emprove_selectWorstRanked = command.add_parser (
    "selectWorstRanked", description="selectWorstRanked", help='select worst scoring particles'
)
emprove_selectWorstRanked.add_argument("--i", required=True, type=str, help="input file")
emprove_selectWorstRanked.add_argument("--num", required=True, type=str, help="number of ranked particles to select")
emprove_selectWorstRanked.add_argument("--tag", required=False, type=str, help="tag to select worst")
emprove_selectWorstRanked.add_argument("--o", required=False, default=None,  type=str, help="output file")
def selectWorstRanked(args):
    inputFile=args.i
    tag=str(args.tag)
    outputFile=args.o
    if (outputFile==None):
        outputFile=inputFile
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    columns=starHandler.header_columns(inputFile)
    if not args.tag:
        args.tag = args.tag if args.tag else find_emprove_ranking_tag(columns)
    if args.tag not in columns:
        print("target tag [", args.tag,"] not in ", inputFile)
        exit(0)
    print ("target tag=", args.tag)
    starHandler.extractWorst(inputFile, outputFile, int(args.num), args.tag)



emprove_selectBestRanked = command.add_parser (
    "selectBestRanked", description="selectBestRanked", help='select best particles'
)
emprove_selectBestRanked.add_argument("--i", required=True, type=str, help="input file")
emprove_selectBestRanked.add_argument("--num", required=True, type=str, help="number of ranked particles to select")
emprove_selectBestRanked.add_argument("--tag", required=False, type=str, help="tag to select best")
emprove_selectBestRanked.add_argument("--o", required=False, default=None,  type=str, help="output file")
def selectBestRanked(args):
    inputFile=args.i
    outputFile=args.o
    if (outputFile==None):
        outputFile=inputFile
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    columns=starHandler.header_columns(inputFile)
    if not args.tag:
        args.tag = args.tag if args.tag else find_emprove_ranking_tag(columns)
    if args.tag not in columns:
        print("target tag [", args.tag,"] not in ", inputFile)
        exit(0)
    print ("target tag=", args.tag)
    starHandler.extractBest(inputFile, outputFile, int(args.num), args.tag)


#def removeStarDuplicates(filenameIn):
#    version = starHandler.infoStarFile(filenameIn)[2]
#    target_section = "particles"
#    if version == 'relion_v30':
#        target_section = "images"
#    sections = starHandler.read_star_sections(filenameIn)
#    print(sections)

emprove_removeDuplicates = command.add_parser (
    "removeDuplicates", description="removeDuplicates", help='remove duplicate particles from star file'
)
emprove_removeDuplicates.add_argument("--i", required=True, type=str, help="input file")
emprove_removeDuplicates.add_argument("--o", required=False, default=None,  type=str, help="output file")
def removeDuplicates(args):
    inputFile=args.i
    outputFile=args.o
    if (outputFile==None):
        outputFile=inputFile
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    starHandler.removeStarDuplicates(inputFile,outputFile)



emprove_assignClassName = command.add_parser (
    "assignClassName", description="assignClassName", help='assign class name'
)
emprove_assignClassName.add_argument("--i", required=True, type=str, help="input file")
emprove_assignClassName.add_argument("--className", required=True, type=str, help="Class Name (usually a number)")
emprove_assignClassName.add_argument("--o", required=False, default=None,  type=str, help="output file")
def assignClassName(args):
    inputFile=args.i
    outputFile=args.o
    if (outputFile==None):
        outputFile=inputFile
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    starHandler.AssignStarClassName(inputFile,outputFile, args.className)



emprove_changeParamValue = command.add_parser (
    "changeParamValue", description="changeParamValue", help='change value of a certain column (parameter) in a star file'
)
emprove_changeParamValue.add_argument("--i", required=True, type=str, help="input file")
emprove_changeParamValue.add_argument("--columnName", required=True, type=str, help="Name of The column(parameter) to be changed")
emprove_changeParamValue.add_argument("--newValue", required=True, type=str, help="new value to be inserted")
emprove_changeParamValue.add_argument("--o", required=False, default=None,  type=str, help="output file")
def changeParamValue(args):
    inputFile=args.i
    outputFile=args.o
    if (outputFile==None):
        outputFile=inputFile
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    starHandler.changeStarFileParamValue(inputFile, outputFile, args.columnName, args.newValue)


emprove_removeParamValue = command.add_parser (
    "removeParamValue", description="removeParamValue", help='remove particles in a star file where a parameter is of a certain value'
)
emprove_removeParamValue.add_argument("--i", required=True, type=str, help="input file")
emprove_removeParamValue.add_argument("--columnName", required=True, type=str, help="Name of The column(parameter) to be removed")
emprove_removeParamValue.add_argument("--value", required=True, type=str, help="value to be removed")
emprove_removeParamValue.add_argument("--o", required=False, default=None,  type=str, help="output file")
def removeParamValue(args):
    inputFile=args.i
    outputFile=args.o
    if (outputFile==None):
        outputFile=inputFile
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    starHandler.deleteStarFileParamValue(inputFile, outputFile, args.columnName, args.value)


emprove_plotParamValue = command.add_parser (
    "plotParamValue", description="plotParamValue", help='plot values in a star file at a certain column'
)
emprove_plotParamValue.add_argument("--i", required=True, type=str, help="input file")
emprove_plotParamValue.add_argument("--columnName", required=True, type=str, help="Name of The column(parameter) to be changed")
def plotParamValue(args):
    inputFile=args.i
    if (not os.path.isfile(inputFile)):
        print('ERROR: file \"',inputFile,'\" not existing')
        exit(0)
    starHandler.plotStarFileParamValueInteractive(inputFile, args.columnName)




def main(command_line=None):
    args = emprove_parser.parse_args(command_line)
    if args.command == "scoreParticles":
        scoreParticles(args)
    elif args.command == "rankParticles":
        rankParticles(args)
    elif args.command == "scoreClassify":
        scoreClassify(args)
    elif args.command == "produce_reconstructions_script":
        produce_reconstructions_script(args)
    elif args.command == "analyse_reconstructions_outcome":
        analyse_reconstructions_outcome(args)
    elif args.command == "selectWorstRanked":
        selectWorstRanked(args)
    elif args.command == "selectBestRanked":
        selectBestRanked(args)
    elif args.command == "removeDuplicates":
        removeDuplicates(args)
    elif args.command == "assignClassName":
        assignClassName(args)    
    elif args.command == "changeParamValue":
        changeParamValue(args)  
    elif args.command == "removeParamValue":
        removeParamValue(args)  
    elif args.command == "plotParamValue":
        plotParamValue(args)  
    else:
        emprove_parser.print_help()







if __name__ == "__main__":
    main()
