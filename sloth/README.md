# sloth/
## core-methods

Within this directory (`sloth/`) the core-methods of **SLOTH** are located.  
To structure those methods, there are different python files collecting methods for different topics. So for example input/output (IO) related methods as `createNetCDF()` are stored in `IO.py`, where as more complex methods as individual classes like `GRDCdataset()` are stord in standalone python files (`GRDCdataset.py`).  
This structure does also show up in the way individual functions are accessable within **SLOTH**. That for example `createNetCDF()` could be accessed by `sloth.IO.createNetCDF()`, similar to the path this method is located `sloth/IO.py -> createNetCDF()`. This should become more clear inspecting the example scripts under `examples/`.  

Currently there is no extra documentation of all the methods inside **SLOTH**, but all methods are assigned with a so called `docstring` (or at least it is aimed that all methods does have a docstring at some point in the future). This docstring does provide the user with short information of:  
-) What does this method do  
-) What does this method need as input  
-) What does this methos return  

> Important:  
It is to be noted, that the structure of this directory might change in the future!
