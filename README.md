# Welcome to SLOTH

Check also our **[documentation](https://hpscterrsys.github.io/SLOTH/README.html)!**  

This repository does holed smaller and bigger **helper** scripts for simulations 
based on [TSMP](https://www.terrsysmp.org/). The overall idea is to support the 
analysis progress of [TSMP](https://www.terrsysmp.org/)-simulations by providing  
easy accessible functions and methods helping the user to focus on the real 
analysis task. Further, example-scripts should provide ideas and hints of how to 
tackle different analysis steps if those are not easily mapped within a modular 
function or method.  

In general it is to be noted, that the **SLOTH**-repository is not a full 
collection of analysis-scripts, but a living repository, aimed to grow with 
upcoming tasks and providing the found solutions in a prepared and documented 
way for everyone who is facing a similar task at a later time.    
Further, the **SLOTH**-repository is not aimed to act as a single solution for 
analysis tasks, but as supporting lib.


## Getting Started

### Prepare the SLOTH repository
As this repository could make use of submodules (repositories inside of 
repositories are called submodules) a little extra treatment is needed to clone 
this repo. Basically there are two options:

**Option 1**
Clone the repository as usual 
``` bash
git clone https://github.com/HPSCTerrSys/SLOTH.git
```
and initialize and update the submodules afterwards
``` bash
cd SLOTH
git submodule init 
git submodule update
```
>**Note:** you have to type the credentials for each submodule!

**Option 2**
Combine all steps of option 1 in one command:
``` bash
git clone --recurse-submodules https://github.com/HPSCTerrSys/SLOTH.git
```
> **Note:** you have to type the credentials for each submodule!


### Use the SLOTH repository  
All scripts inside of the **SLOTH**-repository are developed and tested on 
[JURECA-DC](https://www.fz-juelich.de/en/ias/jsc/systems/supercomputers/jureca) 
with a default tool-chain based on `Stage2020`, which is provided under  
`/p/project/cslts/local/jureca/env_ini.JURECA.stage2020.GCC`.  
So to use **SLOTH** you first have to source this environment file:  
```
source /p/project/cslts/local/jureca/env_ini.JURECA.stage2020.GCC
```
  
To use **SLOTH** within other projects, you have to extend your local 
`PYTHONPATH`, to tell python where to find **SLOTH**. You can do this by:  
```
cd SLOTH  
export PYTHONPATH=$PYTHONPATH:$(pwd)
```
Afterwards you can simply import **SLOTH** inside any of your python scripts 
e.g. by:
```
import sloth.IO
```
You find this also within the example-scripts.

## How to contribute
You can contribute to **SLOTH** with your own functions, methods, and classes.  
To do so, simply create a new branch and upload an example script to the 
`examples/` directory within the new branch. This way everything stays clean 
and the master-branch is not messed up. Within this new branch you can develop 
whatever you want, as you are not interfering with the master branch.   
If you are done with developing open a PR (Pull Request), and someone from the 
maintaining team wil check your request and decide if your development is ready 
to go to the master-branch and therefore to be part of **SLOTH**.  

To do so:  
```
# clone SLOTH as described within the 'Getting Started' section
cd SLOTH
git checkout -b YOURDEVELOPBRANCHNAME
# Start developing
# [...]
# If ready push to GitLab (ALWAYS PULL BEFORE!)
git pull origin YOURDEVELOPBRANCHNAME
git push origin YOURDEVELOPBRANCHNAME
# And open a PR from the GitHUb website.
```
