# EMprove

emprove is a program to select a subset of high quality particle from a larger subset.

## Simple particle selection
Simple particle selection allows to simply rank and select a stack of particle.
Pro: it is simple to run
Cons: 

## Create a session for iterative particle selection using SCI scoring function for selection, and script files to run
go in the directory where you want the new project to be create and execute the following command:
```
emprove_session_manager new_session \
            --name J182_emprove \
            --particles J185_particlesStack.star \
            --map maps/J182_003_volume_map_half_A.mrc \
            --map2 maps/J182_003_volume_map_half_B.mrc \
            --mask maps/J185_maskLocal.mrc \
            --angpix 0.8400
```
it creates a directory with the given name ("J182_emprove" in the example), and a script file to run, in this case J182_emprove_run.sh


## Create a session for random particle selection, and script files to run
go in the directory where you want the session to be created and execute the following command:
```
emprove_session_manager random_selection_session \
                        --name J188_fullMask_random \
                        --particles J188_fullMask_emprove_target.star \
                        --mask maps/J188_mask0577_dilated.mrc \
                        --dir J188_fullMask_emprove
./J188_fullMask_emprove/J188_fullMask_random_run.sh
```
