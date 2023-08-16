# Cross-Correlation Analysis Toolbox

This folder contains MATLAB functions for estimating cross-correlations between pairs of neural responses and extracting noise/trial-by-trial correlations using a shuffling procedure. 
It includes two main functions: `SetCrossCorrParams` defines options for pair-wise correlations such as over which subset of the data to compute correlations or variables
that should be accounted for to subtract the influence of such variables from the overall correlation and estimate the trial-by-trial correlation. `CrossCorrAnalysis` then takes these options as an input,
together with a structure containing the variables to account for (stored in `Nav`) and an array of responses from which pair-wise correlations will be computed (e.g `Spk.spikeTrain`).
See the Load function section to see how the `Nav` and `Spk` structures should be constructed.

## SetCrossCorrParams

The `SetCrossCorrParams` function defines parameters needed for estimating cross-correlations between pairs of neural responses. 
It takes explanatory variables (stored in `Nav`) and response data (e.g. `Spk.spikeTrain`) as inputs and outputs a structure `crossparams` containing analysis parameters.

Parameters included in `crossparams`:

- `subset`: Conditions and values used to subset the data in `Nav` and over which correlations will be estimated.
- `cellidx`: Subset of cells for which correlations will be computed.
- `sampleRate`: Sampling rate of the data.
- `lag`: Range of the time lags over which cross-correlations will be computed (in seconds).
- `nShuffle`: Number of shuffle controls for establishing a distribution of correlation values and estimate correlations due to the variables indicated in `crossparams.variablenames`.
- `variablenames`: Cell array of variable names in `Nav` to account for when computing the signal correlations.
- `binedges`: Cell array of bin edges for discretizing variables indicated in `crossparams.variablenames`.
- `nspk_th`: Minimal number of spikes to consider a cell.
- `timewin`: Spike count window in seconds.

## CrossCorrAnalysis

The `CrossCorrAnalysis` function estimates cross-correlations between pairs of responses.
It uses a shuffling procedure within bins of the variables defined in `crossparams.variablenames` to estimate signal correlations and then computes trial-by-trial (noise) correlations
by subtracting the signal correlations from the overall correlations.

Fields included in the `Cross` structure:

- `crossparams`: Parameters used for cross-correlation analysis.
- `lagbins`: Time bins for the cross-correlation lags.
- `ccAll`: Pair-wise cross-correlation of original signals.
- `ccNoise`: Pair-wise noise or trial-by-trial correlations (ccAll - ccSig).
- `ccSig`: Pair-wise cross-correlation expected from shared selectivity to variables indicated in `crossparams.variablenames`.
- `ccSigSD`: Standard deviation of signal cross-correlation estimated by shuffling.
- `ccNoiseSD`: Standard deviation of trial-by-trial correlations.
- `pval`: P-value matrix for the peak of the trial-by-trial cross-correlation.
- `bestcc`: Maximum cross-correlation value.
- `bestlag`: Lag corresponding to the maximum cross-correlation value.

## Usage Example

1. **Load behavioral data and spike data:**
   
    `datapath = 'pat/to/your/data'`

    `loadparams = SetLoadParams(datapath);`

    `Nav = LoaddataNav(loadparams);`

    `Spk = LoaddataSpk(loadparams, Nav.sampleTimes);`

    `Srep = Spk.spikeTrain;`

2. **Define parameters using SetCrossCorrParams:** 

   `crossparams = SetCrossCorrParams(Nav, Srep);`

3. **Modify paramters in `crossparams` if needed. For instance:**

   `crossparams.subset.Condition = [1 3 5];`
   
    `crossparams.subset.Condition_op = 'ismember';`
   
    `crossparams.subset.Spd = 2.5;`
   
    `crossparams.subset.Spd = '>=';`

    `crossparams.variablename = {'Xpos', 'Spd'};`
   
    `crossparams.binedges = {0:2:100, [0:5:50 inf];`

4. **Perform cross-correlation analysis using CrossCorrAnalysis**

   `Cross = CrossCorrAnalysis(Nav, Srep, crossparams);`


For more detailed information and usage examples of each function, please refer to the function documentation, `Tutorial3.1` and `Tutorial_handson`.

Developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
