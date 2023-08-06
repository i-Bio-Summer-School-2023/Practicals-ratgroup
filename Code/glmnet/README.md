# Glmnet recomplied for Matlab 2022a in Windows 10 64bit
This package was orginally downloaded from https://github.com/junyangq/glmnet-matlab. I then recomplied this orignal package for MATLAB 2022a in Windows 10.  
The recomplied package runs smoothly. In particular, it solves the issue that the orginal package makes connected workers in parallel computing to crash and abort. This recomplied package can be plug-n-play for the above condition, but has not been tested in other conditions. 

# Process to recomplie glmnet
1. Install [Microsoft Visual Studio Community 2019](https://docs.microsoft.com/en-us/visualstudio/releases/2019/release-notes) (select Desktop development with C++)
2. Install [Intel® oneAPI Base & HPC Toolkit](https://www.intel.com/content/www/us/en/developer/tools/oneapi/fortran-compiler.html#gs.umuow8) (You don't need to install everything, but make sure to select Fortran complier when installing the HPC Toolkit)  
		This step is slightly different from previous descriptions, as Intel recently retired the previous package Intel® Parallel Studio XE, which contains the Fortran complier that we need. However, Intel made Fortran complier available for free in the new Intel® oneAPI HPC Toolkit. I downloaded and installed both the Base and HPC Toolkit, but what you really need is just the Fortran complier from the HPC Toolkit. 
3. In MATLAB, navigate to the glmnet folder, run the following code
```
mex -v COMPFLAGS='$COMPFLAGS /real_size:64 /integer_size:64' glmnetMex.F GLMnet.f
```
Reference:
https://github.com/lachioma/glmnet_matlab

# Glmnet in Matlab (original description)
Lasso and elastic-net regularized generalized linear models. See https://web.stanford.edu/~hastie/glmnet_matlab/

We updated the package with compiled Mex-files for the newer verions of Matlab (2020a/b). Previous version is stored in the [matlab/R2013_older](https://github.com/junyangq/glmnet-matlab/tree/matlab/R2013_older) branch.

The new version was tested with Matlab 2020b on Mac OS 11 and with Matlab 2020a on CentOS Linux 7 (Core).
