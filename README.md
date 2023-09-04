# Welcome to SLOTH

Check also our **[documentation](https://hpscterrsys.github.io/SLOTH/README.html)!**  

This repository does holed smaller and bigger **helper** scripts for simulations based on [TSMP](https://www.terrsysmp.org/). The overall idea is to support the analysis progress of [TSMP](https://www.terrsysmp.org/)-simulations by providing easy accessible functions and methods helping the user to focus on the real analysis task. Further, example-scripts should provide ideas and hints of how to tackle different analysis steps if those are not easily mapped within a modular function or method.  

In general it is to be noted, that the **SLOTH**-repository is not a full 
collection of analysis-scripts, but a living repository, aimed to grow with 
upcoming tasks and providing the found solutions in a prepared and ocumented way for everyone who is facing a similar task at a later time.    
Further, the **SLOTH**-repository is not aimed to act as a single solution for analysis tasks, but as supporting lib.


## Getting Started

### Prepare the SLOTH repository
As this repository could make use of submodules (repositories inside of 
repositories are called submodules) a little extra treatment is needed to clone this repo. Basically there are two options:

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
[JURECA-DC](https://www.fz-juelich.de/en/ias/jsc/systems/supercomputers/jureca) with a default tool-chain, which is provided under `/p/project/cslts/local/jureca/`.  
So to use **SLOTH** you first have to source this environment file:  
```
source /p/project/cslts/local/jureca/ONEOFTHEDEFAULTENVFILES
```
  
To use **SLOTH** within other projects, you have to extend your local 
[PYTHONPATH](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONPATH), to tell python where to find **SLOTH**. You can do this by:  
```
cd SLOTH  
export PYTHONPATH=$PYTHONPATH:$(pwd)
```
Afterwards you can simply import **SLOTH** components from inside any of your python scripts 
e.g. by:
```
import sloth.IO
```
You find this also within the example-scripts.

## How to contribute
You can contribute to **SLOTH** with your own functions, methods, and classes. How to do so is explained withtin the [How to contribute](./content/maintenance.md#how-to-contribute) section. 
