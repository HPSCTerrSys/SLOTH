[[_TOC_]]

# Welcome to SLOTH
This repository does holed smaller and bigger **helper** scripts for simulations based on [TSMP](https://www.terrsysmp.org/).  
The overall idea is to support the analysis progress of [TSMP](https://www.terrsysmp.org/)-simulations by providing  easy accessible functions and methods helping the user to focus on the real analysis task.  
Further, example-scripts should provide ideas and hits of how to tackle different analysis steps if those are not easily mapped within a modular function or method.  
All function, methods, and examples collected within the **SLOTH**-repository does work out-of-the-box by using a provided default data-set under `/p/project/cslts/shared_data` to enable an direct usage.

In general it is to be noted, that the **SLOTH**-repository is not a full collection of analysis-scripts, but a living repository, aimed to grow with upcoming tasks and providing the found solutions in a prepared and documented way for everyone who is facing a similar task at a later time.  
Further, the **SLOTH**-repository is not aimed to act as a single solution for analysis tasks, and therefore make use of already existing repositories living in IBG-3 [GitLab](https://icg4geo.icg.kfa-juelich.de/) to avoid redundancies.  

All questions regarding **SLOTH** and its content could be addressed to:  
Name: Niklas WAGNER  
E-mail: n.wagner@fz-juelich.de

## Getting Started

### Prepare the SLOTH repository
As this repository make use of submodules (git-repositories inside of git-repositories are called submodules) a little extra treatment is needed to clone this repo. Basically there are two options:

**Option 1**
Clone the repository as usual 
``` bash
git clone https://icg4geo.icg.kfa-juelich.de/SoftwareTools/SLOTH.git
```
and initialize and update the submodules afterwards
``` bash
cd SLOTH
git submodule init 
git submodule update
```
>**Note:** you have to type the GitLab password for each submodule!

**Option 2**
Combine individual steps of Option 1:
``` bash
git clone --recurse-submodules https://icg4geo.icg.kfa-juelich.de/SoftwareTools/SLOTH.git
```
> **Note:** you have to type the GitLab password for each submodule!

Now the repository is ready to test.

### Use the SLOTH repository  
All scripts inside of the **SLOTH**-repository are developed and tested with a default tool-chain based on `Stage2020` and provided under  
`/p/project/cslts/local/juwels/env_ini.JUWELS.stage2020.GCC`.  
So to use **SLOTH** you first have to source this environment file:  
```
source /p/project/cslts/local/juwels/env_ini.JUWELS.stage2020.GCC
```  
Afterwards all scripts inside of **SLOTH** are working out-of-the-box. Test this e.g. by moving to `examples/` and run one of the example-scripts:  
```
cd SLOTH/examples/  
python ex_SanityCheck_Season.py
``` 
To use **SLOTH** outside the repository, as for example in you own workflow, you have to extend your local `PYTHONPATH`, to tell python where to find **SLOTH**. You can do this by:  
```
cd SLOTH  
export PYTHONPATH=$PYTHONPATH:${pwd}
```
Afterwards you can simply import **SLOTH** inside any of your python scripts by:
```
import sloth
```
You find this also within the example-scripts.

## How to contribute
You can contribute to **SLOTH** with your own functions, methods, and classes which does calculate a specific quantity or fulfill a specific task.  
To do so, simply create a new branch and upload an example script to the `examples/` directory within the new branch. This way everything stays clean and the master-branch (which should be a working branch) is not messed up. Within this new branch you can develop whatever you want, as you are not interfering with the master branch.  
If you are done with developing / uploading contact the maintainer of **SLOTH**, who will decide then if your development is ready to go to the master-branch and therefore to be part of **SLOTH**.  

To do so:  
```
# clone SLOTH as described within the 'Getting Started' section
cd SLOTH
git checkout -b YOURDEVELOPBRANCHNAME
cd examplse/
# Start developing
# [...]
# If ready push to GitLab (ALWAYS PULL BEFORE!)
git pull origin YOURDEVELOPBRANCHNAME
git push origin YOURDEVELOPBRANCHNAME
```
