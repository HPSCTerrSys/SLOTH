# sloth/
## core-methods

Within this directory (`sloth/`) the core-methods of **SLOTH** are located.  
To structure those methods, there are different python files collecting methods for different topics. So for example general / non-specific methods as `createNetCDF()` are stored in `toolBox.py`, where as more specific and more complex methods as the class `GRDCdataset()` are stord in seperate python files (`GRDCdataset.py`).  
This structure does also show up in the way individual functions are accessable within **SLOTH**. That for example `createNetCDF()` could be accessed by `sloth.toolBox.createNetCDF()`. This should become more clear inspecting the example scripts under `examples/`.  

Currently there is no extra documentation of all the methods inside **SLOTH**, but all methods are assigned with a so called `docstring` (or at least it is aimed that all methods does have a docstring at some point in the future). This docstring does provide the user with short information of:  
-) What does this method do  
-) What does this method need as input  
-) What does this methos return  

> Important:  
It is to be noted, that the structure of this directory might change in the future!
