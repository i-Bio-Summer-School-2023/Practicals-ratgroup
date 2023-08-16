# Mapping neural responses along one or two variables

This folder contains MATLAB functions for estimating tuning curves along up to 2 dimensions.It includes two main functions: `SetMapsParams` for defining options such as
which variables to use, how to bin them and over which subset of the data these tuning curves should be estimated. `MapsAnalysis` then takes these options as an input,
together with a structure containing your explanatory variables (called Nav here) and an array of responses for which you want to compute tuning curves (e.g Spk.spikeTrain).
See the Load function section for how the `Nav` and `Spk` structures should be constructed.

## SetMapsParams

The `SetMapsParams` function is used to define parameters required for computing tuning curves. 
It takes behavioral data (`Nav`) and neural response data (`Srep`) as inputs and outputs a structure `mapsparams` containing a set of options to estimate tuning curves.

Parameters included in `mapsparams`:
- `subset`: Define conditions and their values for data subsetting.
- `cellidx`: Logical array to select a subset of cells in `Srep` for analysis.
- `sampleRate`: Sampling rate of the data (in Hz).
- `scalingFactor`: Scaling factor applied to response data (typically 1 / samplingRate).
- `Xvariablename`, `XsmthNbins`, `Xbinedges`: Parameters for X-axis mapping.
- `Yvariablename`, `YsmthNbins`, `Ybinedges`: Parameters for Y-axis mapping.
- `occ_th`: Occupancy threshold for place field estimation.
- `nspk_th`: Minimum number of spikes required for cell inclusion.
- `nShuffle`, `kfold`: Parameters for shuffling and cross-validation.

## MapsAnalysis

The `MapsAnalysis` function estimates tuning curves up to two dimensions. It also estimate their significance using either shuffling or model comparison on cross-validated predictions. 
It takes behavioral data (`Nav`), neural response data (`Srep`), and `mapsparams` as inputs and outputs a structure `Maps` containing various results of the tuning curve analysis.

Results included in `Maps`:
- `map`, `map_cv`, `map_SE`: tuning curves and their standard errors.
- `Xbincenters`, `Ybincenters`: Bin centers along X and Y axes.
- `occmap`: Occupancy map.
- `SI`, `SparsityIndex`, `SelectivityIndex`, `DirectionalityIndex`: metrics to quantify selectivity.
- `SI_pval`, `SparsityIndex_pval`, `SelectivityIndex_pval`, `DirectionalityIndex_pval`: P-values for those metrics.
- `EV`, `EV_cst`, `LLH`, `LLH_cst`: Cross-validated metrics and likelihoods.
- `LLH_pval`: P-values for model comparison.

## Usage

1. **Load behavioral data and spike data:**

    datapath = 'pat/to/your/data'

    `loadparams = SetLoadParams(datapath);`

    `Nav = LoaddataNav(loadparams);`

    `Spk = LoaddataSpk(loadparams, Nav.sampleTimes);`

    `Srep = Spk.spikeTrain;`

3. **Define parameters using SetMapsParams:**
 
    `mapsparams = SetMapsParams(Nav, Srep);`

5. **Modify paramters in `mapsparams` if needed. For instance:**

    `mapsparams.subset.Condition = [1 3 5];`
   
    `mapsparams.subset.Condition_op = 'ismember';`
   
    `mapsparams.subset.Spd = 2.5;`
   
    `mapsparams.subset.Spd = '>=';`

    `mapsparams.Xvariablename = 'Xpos';`
   
    `mapsparams.Xbinedges = 0:2:100;`

    `mapsparams.Yvariablename = [];`

7. **Perform place field analysis using MapsAnalysis**

    `Maps = MapsAnalysis(Nav, Srep, mapsparams);`

For more detailed information and usage examples of each function, please refer to the function documentation, `Tutorial2.1` and `Tutorial_handson`.

Developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
