#!/usr/bin/python3

import argparse
import os
import toml
import pandas as pd

#from emprove import utils
import numpy as np
import scipy.stats as stats
import stat
from emprove import starHandler

emprove_parser = argparse.ArgumentParser(
    prog="emprove_session_manager",
    usage="%(prog)s [command] [arguments]",
    formatter_class=argparse.RawDescriptionHelpFormatter,
)

command = emprove_parser.add_subparsers(dest="command")

def merge_and_sort_csv(file1, file2, output_file):
    # Read the CSV files
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Merge the two dataframes
    merged_df = pd.concat([df1, df2])

    # Sort by 'numParticles'
    sorted_df = merged_df.sort_values(by='numParticles')

    # Save the sorted dataframe to a new CSV file
    sorted_df.to_csv(output_file, index=False)

    return sorted_df





##########################################
##########################################
##### emprove_produce_reconstructions_script 
##########################################
emprove_produce_reconstructions_script = command.add_parser (
    "produce_reconstructions_script", description="produce_reconstructions_script", help='Reconstruction Evaluation'
)
emprove_produce_reconstructions_script.add_argument("--i", required=True, type=str, help="input star file with particles to score")
emprove_produce_reconstructions_script.add_argument("--outDir", required=False, default='./', type=str, help="output Dir")
emprove_produce_reconstructions_script.add_argument("--tagRank", required=False, type=str, help="tag for the ranked particles")
emprove_produce_reconstructions_script.add_argument("--mask", required=True, type=str, help="mask to use for evaluation")
emprove_produce_reconstructions_script.add_argument("--manualParticleSubsets", required=False, type=str, default=None,  help="Comma separated list of manual particles subsets")
emprove_produce_reconstructions_script.add_argument("--scriptName", required=False, default='script_reconstructions.sh', type=str, help="script name as output")
emprove_produce_reconstructions_script.add_argument("--resultFilename", required=False, default='bestRanked_locres_values.csv', type=str, help="filename with results")
emprove_produce_reconstructions_script.add_argument("--mode", required=False, default='bestRanked', type=str, help="mode for selecting particles choose: bestRanked(default), or random")
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

    file_selection_suffix="best"
    if args.mode=="random":
        file_selection_suffix="random"



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
    outFile=os.path.join(args.outDir,"norm_"+str(os.path.split(args.outDir)[-1])+"_"+file_selection_suffix+"${numParticles}")
    if args.mode=="bestRanked":
        reconstruction_command += "    emprove selectBestRanked --i "+args.i+" --o "+outFile+".star --num  ${numParticles} \n"
    elif args.mode=="random":
        reconstruction_command += "    emprove selectRandom --i "+args.i+" --o "+outFile+".star --num  ${numParticles} \n"
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



    ########################
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
            mpirun --np 28 relion_postprocess_mpi --i ${file_basename}_recH1.mrc --i2 ${file_basename}_recH2.mrc --o  ${file_basename} --locres #--locres_thresholdFSC 0.5
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
    
    outFile=os.path.join(args.outDir,"norm_"+str(os.path.split(args.outDir)[-1])+"_"+file_selection_suffix+"${numParticles}")
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
    reconstruction_file_path=os.path.join(args.outDir,args.scriptName)
    with open(reconstruction_file_path, 'a') as f:
        f.write(locres_command)



    ########################
    #####ASSESS LOCRES 
    assess_command = '''
############################
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
    result_filename=os.path.join(args.outDir,args.resultFilename)
    assess_command+="numParticlesCsv="+numParticlesListStr+"\n"
    assess_command+="echo numParticles,max,highQuartile,mean,lowQuartile,min >"+result_filename+"\n"
    assess_command += '''
for numParticles in $(echo $numParticlesCsv | sed "s/,/ /g")
do\n'''
    assess_command+='printf "%s," ${numParticles} >> ' + result_filename +'\n'
    outFile=os.path.join(args.outDir,"norm_"+str(os.path.split(args.outDir)[-1])+"_"+file_selection_suffix+"${numParticles}")
    assess_command += "emprove_app_meanMinMax  "+outFile +locresSuffixfix + ".mrc  "+maskLocresFilename+" >> "+ result_filename +"\n"
    assess_command += '''
    echo ${scorelabelToNormalize}
done
    \n\n'''
    reconstruction_file_path=os.path.join(args.outDir,args.scriptName)
    with open(reconstruction_file_path, 'a') as f:
        f.write(assess_command)





def find_emprove_ranking_tag(elements):
    """
    Search for a column name in the dataframe that starts with '_emprove_' and ends with '_norm' followed by a number.
    """
    import re
    for element in elements:
        if re.match(r'_emprove_.*_norm\d+$', element):
            return element
    return None




#######################################################
## GENERATE PARTICLE SELECTION run script
def generate_run_script(fileSettings):
    with open(fileSettings, 'r') as file:
        data = toml.load(fileSettings)
    run_script_cmd  = '#!/bin/bash\n\n'
    run_script_cmd += '''
ensure_directory() {
    local dir="$1"
    [ ! -d "$dir" ] && mkdir -p "$dir"
}

process_iteration() {
    local ite=$1
    local workingDir=$2
    local maskedCropOption=$3

    local tag="_emprove_SCI__${sigma}_scored_selection_${ite}"
    ensure_directory ${workingDir}
    ensure_directory "${workingDir}/${tag}"
    echo "#############" > ${workingDir}/${tag}/log.txt
    echo "ite=${ite}" >> ${workingDir}/${tag}/log.txt
    echo "directory=${tag}" >> ${workingDir}/${tag}/log.txt
    True=1
    False=0
    DO_WORK_WITH_SIGNAL_SUBTRACTION=$False  # or $False
    local fileWithNormalization=${workingDir}/${tag}/scored_selection_${ite}.star #for regular operations\n'''

    if data.get("map2", "None")=="None":
        run_script_cmd += '    emprove scoreParticles --i $targetStar --map $targetMap1  --mask $targetMask --apix ${apix} --sigma ${sigma} --selectionName ${ite} --o ${fileWithNormalization} --mpi '+ data.get("mpi", 5) + ' --rank '+ data.get("numViews", 350) + '\n'
    else:
        run_script_cmd += '    emprove scoreParticles --i $targetStar --map $targetMap1  --map $targetMap2  --mask $targetMask --apix ${apix} --sigma ${sigma} --selectionName ${ite} --o ${fileWithNormalization} --mpi '+ data.get("mpi", 5) + ' --rank '+ data.get("numViews", 350) + '\n'

    run_script_cmd += '''    if [ $DO_WORK_WITH_SIGNAL_SUBTRACTION -eq $True ]; then
    
    	echo "Working with signal subtraction"
    	emprove_app_starProcess --i ${fileWithNormalization}  --invertTagName _backup_rlnImageName _rlnImageName --o ${workingDir}/${tag}/tmp_inverted.star
    	local ${workingDir}/${tag}/tmp_inverted.star #for signal subtraction, replace the file with subtraction with the whole particles
    else
    	echo "Working without signal subtraction (regular)"
    fi


    #produce reconstruction script for selected subset of particles
    emprove_optimizer automaticParticleSubsets --starFile ${fileWithNormalization} --locres ${previousLocresFile} --numSamples '''+data.get("numRecs", 10)+'''  --save ${workingDir}/${tag}/reconstructions_list.csv
    listParticles=$(<${workingDir}/${tag}/reconstructions_list.csv)
    echo "listParticles=$listParticles" >>${workingDir}/${tag}/log.txt
    emprove_session_manager produce_reconstructions_script --i ${fileWithNormalization} --mask ${targetMask} --tagRank ${tag}_norm'''+data.get("numViews", 350)+''' --outDir ${workingDir}/${tag} --manualParticleSubsets $listParticles --scriptName script_reconstructions.sh ${maskedCropOption}
    ./${workingDir}/${tag}/script_reconstructions.sh
    ./${workingDir}/${tag}/script_reconstructions.sh

        
    #produce reconstruction script for the target
    echo "target=$(emprove_optimizer getNumParticles --locres ${workingDir}/${tag}/bestRanked_locres_values.csv --save ${workingDir}/${tag}/target_num_of_particles.csv)" >> ${workingDir}/${tag}/log.txt
    targetNumOfParticles=$(<${workingDir}/${tag}/target_num_of_particles.csv)
    emprove_session_manager produce_reconstructions_script --i ${fileWithNormalization} --mask ${targetMask} --tagRank ${tag}_norm'''+data.get("numViews", 350)+''' --outDir ${workingDir}/$tag --resultFilename target_locres_values.csv --manualParticleSubsets $targetNumOfParticles --scriptName targetParticlesScript.sh ${maskedCropOption}
    listParticles=$(<${workingDir}/${tag}/targetParticlesScript.sh)
    ./${workingDir}/${tag}/targetParticlesScript.sh

    emprove_optimizer generate_overview --directory ${workingDir}
    targetMap1=$(emprove_optimizer getTarget --overviewFile ${workingDir}/overview.txt --map1)
    targetMap2=$(emprove_optimizer getTarget --overviewFile ${workingDir}/overview.txt --map2)
    #targetStar=$(emprove_optimizer getTarget --overviewFile ${workingDir}/overview.txt --particles)
    #sigma=$(emprove_optimizer getTarget --overviewFile ${workingDir}/overview.txt --sigma)

'''

    run_script_cmd += "}\n\n\n"
    run_script_cmd += 'sigma='+data.get("sigma", 1)+'\n'
    run_script_cmd += 'targetStar="'+data.get("particles", "None")+'"\n'
    run_script_cmd += 'targetMap1="'+data.get("map", "None")+'"\n'
    run_script_cmd += 'targetMap2="'+data.get("map2", "None")+'"\n'
    run_script_cmd += 'targetMask="'+data.get("mask", "None")+'"\n'
    run_script_cmd += 'previousLocresFile=None\n'
    run_script_cmd += '#previousLocresFile="'+data.get("session_name")+'/bestRanked_locres_values.csv"\n'
    run_script_cmd += 'workingDir="'+data.get("session_name")+'"\n'
    maskingCrop_string = " --masked_crop " if data.get("maskingCrop", "False") == "True" else " "
    run_script_cmd += 'masked_crop_otpion="'+ maskingCrop_string +'"\n'
    run_script_cmd += 'apix='+data.get("apix", "None")+'\n\n'
    run_script_cmd += 'for ite in {1..'+data.get("maxSelections", 8)+'}; do\n'
    run_script_cmd += '     process_iteration "$ite" "$workingDir" "$masked_crop_otpion" \n'
    run_script_cmd += 'done\n\n\n'


    #runScriptName=str(data.get("session_name"))+"_run.sh"
    runScriptName=os.path.join(data.get("dir", "./"),data.get("session_name"),str(data.get("session_name"))+"_run.sh")

    with open(runScriptName, 'w') as f:
        f.write(run_script_cmd)
    os.chmod(runScriptName, 0o755)







#######################################################
## CREATE NEW SESSION, what the user has to call first
#######################################################
emprove_new_select_session = command.add_parser (
    "new_select_session", description="new_select_session", help='Create a new selection session with experiment parameters'
)
emprove_new_select_session.add_argument("--name", required=True, type=str, help="name of the new session, it creates a dir with that name, and a toml file in that directory with that name")
emprove_new_select_session.add_argument("--particles", required=True, type=str, help="star file with list of particles")
emprove_new_select_session.add_argument("--map", required=True, type=str, help="reference map, can be first half map if map2 is insterted")
emprove_new_select_session.add_argument("--map2", required=False, type=str, help="second half maps")
emprove_new_select_session.add_argument("--mask", required=True, type=str, help="mask for the region")
emprove_new_select_session.add_argument("--angpix", required=False, type=str, help="pixel spacing in Ansgtrom, if not given it is inferred from the map")
emprove_new_select_session.add_argument("--sigma", required=False, type=float, default=1.0, help="sigma used for the SCI score")
emprove_new_select_session.add_argument("--CC", action='store_true', help="use CC score instead of SCI score")
emprove_new_select_session.add_argument("--maskingCrop", action='store_true', help="in the sake of speed, cropping the file according to mask for computing locres")
emprove_new_select_session.add_argument("--maxSelections", required=False, default=8, type=int, help="Max number of selections (by default 8 selections)")
emprove_new_select_session.add_argument("--numRecs", required=False, default=10, type=int, help="num Samplin gReconstructions. Number of reconstructions for each selections. More reconstruction, more precise is the selection")
emprove_new_select_session.add_argument("--randomSeed", action='store_true', help="Same random sequence at each call. Get non-random initialization. Recommended to use only it only for debug or for fully reproducible reconstruction sampling (i.e. it is not completely random, ).")
#emprove_new_select_session.add_argument("--signalSubtraction", action='store_false', help="workWithSignalSubtraction")
#emprove_new_select_session.add_argument("--randomSeed", required=False, default=8, type=int, help="Same random sequene at each call. Get non-random initialization. Recommended to use only it only for debug or for fully reproducible reconstruction sampling (i.e. it is not completely random, ).")
emprove_new_select_session.add_argument("--mpi", required=False, default=5, type=str, help="number of mpi processes")
emprove_new_select_session.add_argument("--numViews", required=False, default=350, type=str, help="number of euler views for ranking")


#emprove_new_select_session.add_argument("--comparisonType", required=False, default="halfmaps", type=str, help="type of selection, select from [halfmaps,singlemap,]")
#emprove_new_select_session.add_argument("--typeSession", required=False, default="halfmaps", type=str, help="type of selection, select from [halfmaps,singlemap,]")


def new_select_session(args):
    # Create the directory if it doesn't exist
    if not os.path.exists(args.name):
        os.makedirs(args.name)

    # Path for the settings file
    settings_file_path = os.path.join(args.name, 'session_settings.toml')
    if not os.path.isfile(args.map):
        print("map file ", args.map, "not valid")

    if not os.path.isfile(args.mask):
        print("mask file ", args.mask, "not valid")



    # Create the settings file if it doesn't exist
    #if not os.path.isfile(settings_file_path):
    with open(settings_file_path, 'w') as file:
            # Write a comment and the session_name variable
            #file.write("# name of the session\n")
            #toml.dump({'session_name': args.name}, file)
            file.write("###########################n\n")
            file.write("\n# name of the session\n")
            file.write(f'session_name = "{args.name}"\n')
            file.write("\n# particles stack\n")
            file.write(f'particles = "{args.particles}"\n')            
            file.write("\n# Pixel size in angstrom for the map\n")
            file.write(f'apix = "{args.angpix}"\n')
            file.write("\n# Path of the mask file\n")
            file.write(f'mask= "{args.mask}"\n')
            file.write("\n# maps. If only one map it is given, emprove uses only this map, but this will break the assumption of independency for the two half maps, and ab-initio angular assignment might be required.\n")
            file.write(f'map = "{args.map}"\n')
            file.write(f'map2 = "{args.map2}"\n')
            file.write("\n# sigma used for the SCI score\n")
            file.write(f'sigma = "{args.sigma}"\n')
            file.write("\n# minimum sigma allowed (advanced option for optimization purposes, suggested to leave as it is)\n")
            file.write(f'minimum_sigma_allowed = 0.6\n')
            file.write("\n# sigma  (advanced option for optimization purposes, suggested to leave as it is)\n")
            file.write(f'sigma_decreasing_step = 0.05\n')
            file.write("\n# maskingCrop, if true crop the file according to mask for computing locres, in the sake of speed\n")
            file.write(f'maskingCrop = "{args.maskingCrop}"\n')
            file.write("\n# maxSelections, max number of selection iterations\n")
            file.write(f'maxSelections = "{args.maxSelections}"\n')
            file.write("\n# numRecs, Number of reconstructions for each selections. More reconstruction, more precise is the selection\n")
            file.write(f'numRecs = "{args.numRecs}"\n')
            file.write("\n# num of mpi processes\n")
            file.write(f'mpi = "{args.mpi}"\n')
            file.write("\n# num of eulerian views\n")
            file.write(f'numViews = "{args.numViews}"\n')



    # Check if the --map2 argument was provided
    #if args.map2 is not None:
    #    print(f"The --map2 argument was provided and its value is: {args.map2}")
    #else:
    #    print("The --map2 argument was not provided.")
    generate_run_script(settings_file_path)






#######################################################
#######################################################
## GENERATE run script for random selection section
def generate_random_run_script(fileSettings):
    with open(fileSettings, 'r') as file:
        data = toml.load(fileSettings)
    run_script_cmd  = '#!/bin/bash\n\n'
    run_script_cmd += '''
ensure_directory() {
    local dir="$1"
    [ ! -d "$dir" ] && mkdir -p "$dir"
}

process_random() {
    local name=$1
    local workingDir=$2
    local particles=$3
    local listSplits=$4
    local maskedCropOption=$5

    local tag=${name}
    ensure_directory ${workingDir}
    ensure_directory "${workingDir}/${tag}"
    emprove_session_manager produce_reconstructions_script --mode random --i ${particles} --mask ${targetMask}  --outDir ${workingDir}/${tag} --manualParticleSubsets $listSplits --scriptName script_reconstructions.sh ${maskedCropOption} --resultFilename randomRanked_locres_values.csv
    ./${workingDir}/${tag}/script_reconstructions.sh
'''

    run_script_cmd += "}\n\n\n"
    run_script_cmd += 'targetStar="'+data.get("particles", "None")+'"\n'
    run_script_cmd += 'targetMask="'+data.get("mask", "None")+'"\n'
    run_script_cmd += 'particles="'+data.get("particles", "None")+'"\n'
    run_script_cmd += 'listSplits="'+data.get("splits", "")+'"\n'
    run_script_cmd += 'workingDir="'+  os.path.join(data.get("dir", "./"),  data.get("session_name", "randomSelection"))  +'"\n'
    run_script_cmd += 'masked_crop_otpion=" --masked_crop "\n'
    run_script_cmd += 'process_random "$ite" "$workingDir" "$particles" "$listSplits" "$masked_crop_otpion" \n'
    run_script_cmd += './${workingDir}/script_reconstructions.sh'
    run_script_cmd += '\n\n\n'

    runScriptName=os.path.join(data.get("dir", "./"),data.get("session_name"),str(data.get("session_name"))+"_run.sh")
    with open(runScriptName, 'w') as f:
        f.write(run_script_cmd)
    os.chmod(runScriptName, 0o755)


#######################################################
## CREATE random selection section
#######################################################
emprove_random_selection_session = command.add_parser (
    "random_selection_session", description="random_selection_session", help='Create a new session with random selection'
)
emprove_random_selection_session.add_argument("--name", required=False, default="random_selection", type=str, help="name of the random_selection section, it creates a dir with that name, in not given it will assign the default name 'random_selection'")
emprove_random_selection_session.add_argument("--dir", required=False, default="./", type=str, help="working directory, if not defined will go with the current directory")
emprove_random_selection_session.add_argument("--particles", required=True, type=str, help="filename with the particles to be created")
emprove_random_selection_session.add_argument("--mask", required=True, type=str, help="filename with the mask")
emprove_random_selection_session.add_argument("--subsets", required=False, type=int, default=10, help="number of subsets to reconstruct (default=10)")
def random_selection_session(args):
    # Create the directory if it doesn't exist
    if not os.path.exists(args.name):
        os.makedirs(args.name)

    # Path for the settings file
    if not os.path.isfile(args.particles):
        print("particles star file ", args.particles, "not valid")

    if not os.path.isfile(args.mask):
        print("mask file ", args.mask, "not valid")

    params=starHandler.readStar(args.particles)
    #print('size=',params.shape[0])
    #print('splits=',args.subsets)
    interval_values = np.linspace(0, params.shape[0], args.subsets)
    interval_values_integer = map(int, interval_values[1:])
    interval_values_string = ",".join(map(str, interval_values_integer))
    #print('values=',interval_values_string)
    settings_file_path=os.path.join(args.name,"session_random_settings.txt")

    with open(settings_file_path, 'w') as file:
            # Write a comment and the session_name variable
            #file.write("# name of the session\n")
            #toml.dump({'session_name': args.name}, file)
            file.write("###########################\n")
            file.write("\n# name of the session\n")
            file.write(f'session_name = "{args.name}"\n')
            file.write("\n# particles stack\n")
            file.write(f'particles = "{args.particles}"\n')            
            file.write("\n# directory\n")
            file.write(f'dir = "{args.dir}"\n')
            file.write("\n# Path of the mask file\n")
            file.write(f'mask= "{args.mask}"\n')
            file.write("\n# number of particles to reconstruct:\n")
            file.write(f'splits= "{interval_values_string}"\n')
    generate_random_run_script(str(settings_file_path))








#######################################################
#######################################################
## GENERATE run script for random selection section
def generate_classification_run_script(fileSettings):
    with open(fileSettings, 'r') as file:
        data = toml.load(fileSettings)
    run_script_cmd  = '#!/bin/bash\n\n'
    run_script_cmd += '''
ensure_directory() {
    local dir="$1"
    [ ! -d "$dir" ] && mkdir -p "$dir"
}

rec_subset() {
        fileIn=$1
        fileOut_basename=$2
        subset=$3
        if [ -f ${fileOut_basename}_recH${subset}.mrc ]; then
            echo "DOING NOTHING: Reconstructed file ${fileOut_basename}_recH${subset}.mrc exists"
        else
            echo "DOING Reconstruction for file ${fileOut_basename}_recH${subset}.mrc"
            relion_reconstruct --i ${fileIn} --o  ${fileOut_basename}_recH${subset}.mrc --subset ${subset} --ctf &
            sleep 40
        fi
}

process_classification() {
    local name=$1
    local workingDir=$2
    local scoringDir=$3
    local classesDir=$4
    local particles=$5
    local angpix=$6
    local sigma=$7
    local mask=$8
    local listEqualizedMaps=$9
    local listClassNames=${10}
    local numMPI=${11}

    ensure_directory ${workingDir}

    
    # Convert space-separated strings to arrays
    IFS=' ' read -r -a mapsArray <<< "$listEqualizedMaps"
    IFS=' ' read -r -a classNamesArray <<< "$listClassNames"

    results_classification_array=""
    for i in "${!mapsArray[@]}"; do
        local map=${mapsArray[$i]}
        local className=${classNamesArray[$i]}
        echo "Processing map: $map with class name: $className"
        emprove scoreParticles --i ${particles} --mask ${mask} --map ${map} --apix ${angpix} --sigma  ${sigma} --o "${scoringDir}/${className}".star --mpi "$numMPI"
        results_classification_array="${results_classification_array} ${scoringDir}/${className}.star"
    done

    emprove_utils scores_to_csv --i $results_classification_array --o ${classesDir}/${name}_classified.star --csv ${scoringDir}/${name}.csv
    
    for i in "${!mapsArray[@]}"; do
        class_number=$((i + 1))
        emprove_utils extract_particles_from_label_value --i  ${classesDir}/${name}_classified.star --o ${classesDir}/class_${class_number}.star --label _rlnClassNumber --value ${class_number}
        rec_subset  ${classesDir}/class_${class_number}.star ${classesDir}/class_${class_number} 1
        rec_subset  ${classesDir}/class_${class_number}.star ${classesDir}/class_${class_number} 2
    done
    
'''
    #emprove_session_manager produce_reconstructions_script --mode random --i ${particles} --mask ${targetMask}  --outDir ${workingDir}/${tag} --manualParticleSubsets $listSplits --scriptName script_reconstructions.sh ${maskedCropOption} --resultFilename randomRanked_locres_values.csv
    #./${workingDir}/${tag}/script_reconstructions.sh

    run_script_cmd += "}\n\n\n"
    run_script_cmd += 'name="'+data.get("session_name", "classification")+'"\n'
    run_script_cmd += 'targetStar="'+data.get("particles", "None")+'"\n'
    run_script_cmd += 'targetMask="'+data.get("mask", "None")+'"\n'
    run_script_cmd += 'particles="'+data.get("particles", "None")+'"\n'
    run_script_cmd += 'listEqualizedMaps="'+data.get("equalized_maps", "")+'"\n'
    run_script_cmd += 'listClassNames="'+data.get("scored_classes", "")+'"\n'
    run_script_cmd += 'workingDir="'+  data.get("dir", "./")  +'"\n'
    run_script_cmd += 'scoring_dir="'+ os.path.join(data.get("dir", "./"),str(data.get("scoring_dir", "scoring_dir"))) +'"\n'
    run_script_cmd += 'classes_dir="'+ os.path.join(data.get("dir", "./"),str(data.get("classes_dir", "classes_dir"))) +'"\n'
    run_script_cmd += 'angpix="'+data.get("angpix", "1.0")+'"\n'
    run_script_cmd += 'sigma="'+data.get("sigma", "1.0")+'"\n'
    run_script_cmd += 'numMPI="'+data.get("numMPI", "8")+'"\n'
    run_script_cmd += "emprove_utils equalize_images  --i "+data.get("input_maps", "None")+ " --o_suffix _" + data.get("equalized_suffix", "equalized")+ " --dir " + os.path.join(data.get("dir", "./"), data.get("equalized_suffix", "equalized")) +"\n"
    run_script_cmd += 'process_classification "$name" "$workingDir" "$scoring_dir" "$classes_dir" "$particles" "$angpix" "$sigma" "$targetMask" "$listEqualizedMaps" "$listClassNames" "$numMPI"\n'
    #run_script_cmd += './${workingDir}/script_reconstructions.sh'
    run_script_cmd += '\n\n\n'

    runScriptName=os.path.join(data.get("dir", "./"),str(data.get("session_name"))+"_run.sh")
    with open(runScriptName, 'w') as f:
        f.write(run_script_cmd)
    os.chmod(runScriptName, 0o755)


#######################################################
## CREATE classification selection section
#######################################################
emprove_classification_session = command.add_parser (
    "classification_session", description="classification_session", help='Create a new session with classification'
)
emprove_classification_session.add_argument("--name", required=False, default="classification", type=str, help="name of the classification_selection section, it creates a dir with that name, in not given it will assign the default name 'random_selection'")
emprove_classification_session.add_argument("--dir", required=False, default="./", type=str, help="working directory, if not defined will go with the current directory")
emprove_classification_session.add_argument("--particles", required=True, type=str, help="filename with the particles to be created")
emprove_classification_session.add_argument("--mask", required=True, type=str, help="filename with the mask")
emprove_classification_session.add_argument("--maps", required=True, nargs='+', type=str, help="files with the input MRC maps")
emprove_classification_session.add_argument("--scoring_dir", required=False, default="scored_classes", type=str, help="directory with the scores to be created")
emprove_classification_session.add_argument("--classes_dir", required=False, default="final_classes", type=str, help="directory with the final classes to be created")
emprove_classification_session.add_argument("--o", required=False, type=str, default="classification.star", help="filename with the updated classes")
emprove_classification_session.add_argument("--angpix", required=True, type=float, default="1", help="pixel spacing")
emprove_classification_session.add_argument("--sigma", required=False, type=float, default="1", help="sigma for the SCI SCORE")
emprove_classification_session.add_argument("--mpi", required=False, type=int, default="8", help="number of MPI values")


def classification_session(args):

    # Path for the settings file
    if not os.path.isfile(args.particles):
        print("particles star file ", args.particles, "not valid")
        return

    if not os.path.isfile(args.mask):
        print("mask file ", args.mask, "not valid")
        return

    for map_file in args.maps:
        if not os.path.isfile(map_file):
            print("Map file", map_file, "not valid")
            return

    # Create the directory if it doesn't exist
    if not os.path.exists(args.dir):
        os.makedirs(args.dir)
    workingDir=os.path.join(args.dir,args.name)
    if not os.path.exists(workingDir):
        os.makedirs(workingDir)
    scores_dir=os.path.join(workingDir,args.scoring_dir)
    if not os.path.exists(scores_dir):
        os.makedirs(scores_dir)
    classes_dir=os.path.join(workingDir,args.classes_dir)
    if not os.path.exists(classes_dir):
        os.makedirs(classes_dir)

    params=starHandler.readStar(args.particles)
    #print('size=',params.shape[0])
    #print('splits=',args.subsets)
    settings_file_path=os.path.join(args.name,"session_classification_settings.txt")
    maps_string = ' '.join(args.maps)
    equalizedDir=os.path.join(workingDir,'equalized')
    modified_maps = [f"{equalizedDir}/{os.path.basename(file).split('.')[0]}_{'equalized'}.mrc" for file in args.maps]
    modified_maps_string = ' '.join(modified_maps)
    scored_classes_names = [f"{os.path.basename(file).split('.')[0]}_{'scoredClass'}" for file in args.maps]
    scored_classes_names_string = ' '.join(scored_classes_names)


    with open(settings_file_path, 'w') as file:
            # Write a comment and the session_name variable
            #file.write("# name of the session\n")
            #toml.dump({'session_name': args.name}, file)
            file.write("###########################\n")
            file.write("\n# name of the session\n")
            file.write(f'session_name = "{args.name}"\n')
            file.write("\n# particles stack\n")
            file.write(f'particles = "{args.particles}"\n')            
            file.write("\n# working directory\n")
            file.write(f'dir = "{workingDir}"\n')
            file.write("\n# Path of the mask file\n")
            file.write(f'mask= "{args.mask}"\n')
            file.write("\n# maps to analyze:\n")
            file.write(f'input_maps= "{maps_string}"\n')
            file.write(f'equalized_maps= "{modified_maps_string}"\n')
            file.write(f'scored_classes= "{scored_classes_names_string}"\n')
            file.write(f'equalized_suffix= "equalized"\n')
            file.write("\n# scoring dir:\n")
            file.write(f'scoring_dir= "{args.scoring_dir}"\n')
            file.write(f'classes_dir= "{args.classes_dir}"\n')
            file.write("\n# output file with updated classification:\n")
            file.write(f'output_classificated_file= "{args.o}"\n')
            file.write("\n# angpix, pixel spacing:\n")
            file.write(f'angpix= "{args.angpix}"\n')
            file.write("\n# sigma for the SCI score:\n")
            file.write(f'sigma= "{args.sigma}"\n')
            file.write("\n# number of MPI processing:\n")
            file.write(f'numMPI= "{args.mpi}"\n')


    generate_classification_run_script(str(settings_file_path))







def main(command_line=None):
    args = emprove_parser.parse_args(command_line)
    if args.command == "produce_reconstructions_script":
        produce_reconstructions_script(args)
    elif args.command == "random_selection_session":
        random_selection_session(args)
    elif args.command == "classification_session":
        classification_session(args)
    elif args.command == "new_select_session":
        new_select_session(args)
    else:
        emprove_parser.print_help()







if __name__ == "__main__":
    main()
