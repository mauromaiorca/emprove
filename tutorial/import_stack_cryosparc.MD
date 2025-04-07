# Import stack from cryoSPARC

EMPROVE, currently uses STAR file-based stack. So, if you need to convert a cryoSPARC stack for use in EMPROVE, you can do as following.
Firstly, locate where your particle stack currently is. For example you need to import a full stack which is result from NU refinement job (e.g. J44), the file you are looking for is the:
J44_006_particles.cs, which contains the locations of the particles in the cryoSPARC directory. You can use the pyem program csparc2star.py

(1) to convert the "particles.cs" file to star file. You can do it by using the command csparc2star.py from pyem. The syntax is something like:
csparc2star.py J44_006_particles.cs J44_006_particles.star
(2) point the star file to the mrcs filename instead of the mrc.You can do it by opening the star file with any decent text editor (gedit on linux works fine) and add an s to all the mrc filenames in the _rlnImageName parameters (something like "edit -> replace all" in the text editor). Or using regular expressions just do:

```
sed -i 's/_particles\.mrc /_particles.mrcs /g' J44_006_particles.star
```

(3) create the stack with relion using the command:
relion_stack_create --i J44_006_particles.star --o full_stack --one_by_one




