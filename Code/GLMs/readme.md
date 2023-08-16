# Generalized Linear Models (GLMs)

This folder contains MATLAB functions for estimating Generalized Linear Models (GLMs) and evaluating their significance in neural recordings. 
It includes two main functions:  `SetGLMsParams` for defining options such as
which variables to use, how to bin them and over which subset of the data these GLMs should be estimated. `GLMsAnalysis` then takes these options as an input,
together with a structure containing the variables to use as predictors (called `Nav` here) and an array of responses for which you want to fit GLMs (e.g `Spk.spikeTrain`).
See the Load function section to see how the `Nav` and `Spk` structures should be constructed.

## SetGLMsParams

The `SetGLMsParams` function is used to define options required for estimating GLMs using the `GLMAnalysis` function.
It takes as an input a structure containing explanatory variables (`Nav`) and neural response data (`Srep`).
The output is a structure `glmsparams` containing parameters for running the GLM analysis.

Parameters included in `glmsparams`:
- `subset`: Define conditions and their values for data subsetting.
- `cellidx`: Logical array to select a subset of cells for analysis.
- `sampleRate`: Sampling rate of the data.
- `scalingFactor`: Scaling factor applied to response data.
- `variablenames`: Names of predictor variables (up to 2).
- `smthNbins`, `binedges`: Smoothing window and bin edges for each predictor.
- `occ_th`: Occupancy threshold for predictor inclusion.
- `nspk_th`: Minimal spikes required for cell inclusion.
- `kfold`: Number of folds for cross-validation.
- `intercept`: Flag for including a constant term in the GLM model.
- `alpha`: Regularization type in glmnet (lasso, ridge, elastic net).
- `maxit`, `thresh`: Convergence parameters for glmnet.
- `pval_th`: P-value threshold for predictor significance.

## GLMAnalysis

The `GLMAnalysis` function estimates tuning curves using a Poisson GLM model. it performs k-fold cross-validation, evaluating model significance.
Model comparison is performed by comparing likelihood values of held-out data under the model.
This function only works up to two explanatory variables. It takes behavioral data (`Nav`), neural response data (`Srep`), and `glmsparams` as inputs and
outputs a structure `GLMs` containing various results about the optimal GLM models.

Results included in `GLMs`:
- `glmsparams`: Input parameters for the GLM analysis.
- `bestmodel`: Best model indicator for each cell.
- `LLH`, `LLH_cst`: Log likelihood values for different models and constant mean model.
- `tuning`: Array of structures for each variable containing the following fields:
  - `bincenters`: Bin centers for the variable.
  - `map`, `mapcv`: Tuning curves for selected cells.
  - `map_SE`: Standard error estimates for tuning curves.
  - `pval`: P-values for variable significance.
Results returned in `tuning` corresponds to the best model as estimated by a likelihood ratio test on held-out data.

## Usage

1. **Load behavioral data and spike data:**

      datapath = 'path/to/your/data'

    `loadparams = SetLoadParams(datapath);`

    `Nav = LoaddataNav(loadparams);`

    `Spk = LoaddataSpk(loadparams, Nav.sampleTimes);`

    `Srep = Spk.spikeTrain;`

2. **Define parameters using SetGLMsParams:**

    `glmsparams = SetGLMsParams(Nav, Srep);`

4. **Modify paramters in `glmsparams` if needed. For instance:**

   `glmsparams.subset.Condition = [1 3 5];`
   
    `glmsparams.subset.Condition_op = 'ismember';`
   
    `glmsparams.subset.Spd = 2.5;`
   
    `glmsparams.subset.Spd = '>=';`

    `glmsparams.variablename = {'Xpos', 'Spd'};`
   
    `glmsparams.binedges = {0:2:100, [0:5:50 inf];`

6. **Estimate GLMs using GLMsAnalysis**

   `GLMs = GLMsAnalysis(Nav, Srep, glmsparams);`

For more detailed information and usage examples of each function, please refer to the function documentation, `Tutorial2.2` and `Tutorial_handson`.

Developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
