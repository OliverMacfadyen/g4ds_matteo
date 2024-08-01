# G4DS 11

## DEVELOPMENT PHASE

> Porting g4ds from GEANT4-10.00 to GEANT4.11.02

Currently working on conpiling problems

* `DSLScintilation.cc` was removed. It was designed for the DS50 veto, but never used.

* **NOTE** `G4Cherenkov.cc` renamed to G4Cherenkov.cc.ori. Old G4Cherenkov class no more compatible. We need to check if old Geant4.10 issues have been solved with Geant4.11. 

* **NOTE** `G4NeutronHPCaptureFS.cc` was modified to account for the most accurate gamma cascade model from neutron capture on Gd. The G4NeutronHPCaptureFS class was substituted with a more generic hadronic class. Need to find another solution. 

* **NOTE** `G4ParticleDefinition.cc` was modified to allow the setting of the ion lifetime to 0. This was needed as g4ds uses  particle global time, the timne elapsed from the generatio of the primary. Scintillation photons induced by radioactive decays have global times equal to isotope decay time + photon propagation. If the decay time is of the order of GYear, the photon time is at the 10^25-10^26 ns, loosing the precision at the ns level. The new approach of G4ParticleDefinition is radically changed, but maybe could work without the custom modification. This has to be checked when g4ds11 will compile.  

* **NOTE** Need to remove any reference to ROOT. See dsparameters (ITO optical model) and ArDM detector construction which uses TVector2. 


## Compiling g4ds 

within g4ds create a folder 

`mkdir build`

enter into the folder 

`cd build`

set the cmake environment 

`cmake ../`

and then compile it 

`make -j`



## Visaulization 

Use the "Qt" flag after the mac name to activate the visaulization (example: ./g4ds vis.mac Qt)

Please refer to the vis.mac macro for a list of useful instructions to setup the geometry navigator 


## Useful git commands

To clone the sources from the last release on the repository:

**git clone git@gitlab.in2p3.fr:darkside/g4ds10.git**

This command will download the last version on the server, will create a folder call g4ds10
and will track modification of g4ds10 on the server. It create a local version of the git server code.

To know the status of your local version (uncommitted files, unaided files, ….)

**git status**

To update your local version from the server one (**be careful with unmarked files**): 

**git pull**


To add a file or a folder in the git repository 
(you cannot commit a file which isn’t added on the repo):

**git add path/to/my/file.ext/or/my/folder**

To commit all your modifications on your local version (**be very careful !!!**):

**git commit -a -m 'your lovely commit message'**

To commit a specific file: 

**git commit -m 'your lovely commit message' path/to/my/file.ext**

To propagate your commit on the server: 

**git push**

To restore a file (erased or modified):

**git checkout path/to/my/file.ext**

**After a commit, you always need to push your modifications !!!**
