# Cell Assembly Detection

This folder contains MATLAB functions for detecting cell assembly patterns from neural spike response data using independent component analysis (ICA). 
It includes two main functions: `SetPattParams` defines options for the detection of cell assemblies such as over which subset of the data to do the detection or the variables 
that should be accounted for if one wants to remove the influence of such variables from the correlated discharge of cells before running the ICA. `PatternAnalysis` then takes these options as an input,
together with a structure containing variables to account for (called `Nav` here) and an array of responses from which assemblies will be detected (e.g `Spk.spikeTrain`).
See the Load function section to see how the `Nav` and `Spk` structures should be constructed.

## SetPattParams

The `SetPattParams` function defines options required for detecting cell assembly patterns from spike response data. 
It takes experimental condition and variable data (stored in `Nav`) and response data (`Srep`) as inputs. 
It outputs a structure `pattparams` containing parameters for the pattern analysis.

Parameters included in `pattparams`:

- `subset`: A structure defining conditions over the fields of Nav for which patterns will be detected.
- `cellidx`: A logical array indicating a subset of cells for pattern detection.
- `nspk_th`: Minimal number of spikes over the training set to consider a cell.
- `Marcenko`: If true, select principal components (PCs) based on the Marcenko-Pastur law.
- `nShuffle`: Number of shuffle controls to perform for eigenvalue selection if `pattparams.Marcenko` is false.
- `variablenames`: Cell array of variable names in Nav used for estimating the signal correlations.
- `binedges`: Cell array of bin edges for discretizing variables in `variablenames`.
- `NoiseCov`: If true, remove signal covariance from the covariance matrix before running the ICA.
- `pvalshf_th`: P-value threshold for selecting PCs when using a shuffling approach (i.e. when 'pattparams.Marcenko = false').
- `strength_th`: Pattern activation threshold to convert strength into "spikes".
- `sampleRate`: Sampling rate of the data (in Hz).
- `timewin`: Size of the spike count window over which to detect correlated firing (in seconds).

## PatternAnalysis

The `PatternAnalysis` function identifies cell assembly patterns in spike response data using independent component analysis (ICA). 
As inputs, it takes a set of data defining experimental conditions and explanatoy variables (stored in `Nav`), spike train responses (`Srep`), and a set of parameters (`pattparams`).
Although the detection of the cell assembly may be done over a subset of the data, `PatternAnalysis` will return the activation strength of those assemblies over the entire dataset.
`PatternAnalysis` computes various pattern analysis results stored in an output structure.

Fields included in the `Patt` structure:

- `pattparams`: Structure of parameters used for pattern detection.
- `weights`: Weight matrix representing the contribution of each neuron to each assembly pattern.
- `cellAssemblies`: Logical matrix indicating which neurons belong to each detected cell assembly.
- `Sparsity`: Sparsity measure of each assembly pattern (between 0 and 1).
- `strength`: Expression strength of each assembly pattern over time.
- `activation`: Events of assembly pattern activation.

## Usage Example

1. **Load behavioral data and spike train responses**

   `datapath = 'path/to/your/data';`

    `loadparams = SetLoadParams(datapath);`

    `Nav = LoaddataNav(loadparams);`

    `Spk = LoaddataSpk(loadparams, Nav.sampleTimes);`

    `Srep = Spk.spikeTrain;`

3. **Define pattern detection parameters using SetPattParams**

    `pattparams = SetPattParams(Nav, Srep);`

    adjust the subset over which assemblies will be detected

    `pattparams.subset.Condition = [1 3 5];`

    `pattparams.subset.Condition_op = 'ismember';`

    `pattparams.subset.Spd = 2.5;`

    `pattparams.subset.Spd_op = '>=';`

4. **Perform pattern analysis using PatternAnalysis**

   `Patt = PatternAnalysis(Nav, Srep, pattparams);`

For more detailed information and usage examples of each function, please refer to the function documentation, `Tutorial3.2` and `Tutorial_handson`.

Developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
