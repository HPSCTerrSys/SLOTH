# examples/
## example scripts showing what is possible with **SLOTH**
This directory (`examples/`) is the heart of the **SLOTH** repository.  
Each script located in this directory does show a real usecase and introcudes how to use individual methods of **SLOTH** to solve them.  
In all scripts lots of comments are used serving as a kind of tutorial / guide.  
Next to many scripts a equal names plot is placed, created by the scritp itself and showing what to epect from the script -- kind of helping to brows through the availabale scritps.  

All scripts are running out of the box on JSC machines (JUWELS and JURECA) as a sample data-set is provided under `/p/scratch/cslts/shared_data/tmp_TestData` which is used by all the example scripts here. This ensures that everyone can clone this repo and directly start testing stuff, without the need to find all dependencies and correct datasets beforhand.  

Even if this directry does contain example scripts showing how to use methods of **SLOTH** the idea of developing **SLOTH** is the other way around. Everything taks with the potential to contribute to **SLOTH** which is coming up in the future, should be solved in a example scripte here. Than, the core-methods which turns out to be essential and resuable by other shoudl be extracted and added to the core-methods of **SLOTH** under the `sloth/` directory.  
This is because many times a 'ready to use' function is hard to imaging beforhand and essential elements does show up with some time onyl.
