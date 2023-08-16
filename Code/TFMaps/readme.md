# Time-Frequency Analysis Toolbox

This folder contains MATLAB funcitons to perform time-frequency analysis of continuous data and compute power and coherence maps along a set of explanaotry variables. 
It includes two main functions: `SetTFMapsParams` for defining options such as
which channels to use, which explanatory variables to use, how to bin them and over which subset of the data these maps should be computed. `SetTFMapsParams` then takes these options as an input,
together with a structure containing your explanatory variables (called `Nav` here) and an array of responses for which you want to compute tuning curves in the spectral domain (e.g `Lfp.Lfp_raw`).
See the Load function section for how the `Nav` and `Lfp` structures should be constructed.

## SetTFMapsParams

The `SetTFMapsParams` function defines a set of parameters needed for time-frequency analysis of continuous data. 
It takes explanatory variables (`Nav`), response data (`Lrep`), and an optional input `sampleTimes_Lrep` (sampling times for signals in `Lrep`). 
It outputs a structure `TFmapsparams` containing parameters for the time-frequency analysis.

Parameters included in `TFmapsparams`:

- `subset`: Conditions and values used to subset the data in `Nav`.
- `freqrange`: Frequency range for time-frequency analysis.
- `chidx`: Subset of signals in `Lrep` considered for analysis.
- `sampleRate`: Sampling rate of the independent variables in `Nav`.
- `sampleRate_raw`: Sampling rate of the signals in `Lrep`.
- `Spectrogram`, `Coherence`: Flags to compute wavelet transforms and coherences.
- `scalingFactor`: Scaling factor on the response data for map computations.
- `Xvariablename`, `Xbinedges`, `XsmthNbins`: Parameters for mapping along X.
- `Yvariablename`, `Ybinedges`, `YsmthNbins`: Parameters for mapping along Y.
- `occ_th`: Occupancy threshold for map computations.
- `nShuffle`, `kfold`: Shuffle and cross-validation controls.
- `parallel`: Parallel computing option.

## TFMapsAnalysis

The `TFMapsAnalysis` function performs time-frequency analysis and computes power and coherence maps based on the parameters defined in `TFmapsparams`. 
It takes as inputs a set of explanatory variables and experimental conditions (stored in `Nav`), response data (a ntimes x nchannels array, e.g. Lfp.Lfp_raw), a set of options 
to run the time-frequency analysis (stored in `TFmapsparams`), and `sampleTimes_Lrep`, the timestamps of samples in the response data array.
This function returns tuning curves along the explanatory variables indicated in `TFmapsparams` for the power and the coherence across frequencies.

Fields included in the `TFMaps` structure:

- `TFmapsparams`: Parameters used for time-frequency analysis.
- `wtMaps`: Cell array of structures containing power maps across frequencies.
- `wcMaps`: Cell array of structures containing coherence maps across frequencies.
- `fqbins`: Frequency bins used for the analysis.

## Usage Example

1. **Load behavioral data and Lfp signals**

    `datapath = 'path/to/your/data';`
   
    `loadparams = SetLoadParams(datapath);`

    `Nav = LoaddataNav(loadparams);`

    `Lfp = LoaddataLfp(loadparams, Nav.sampleTimes);`

    `Lrep = Lfp.Lfp_raw;`

3. **Define time-frequency analysis parameters using SetTFMapsParams**

    `TFmapsparams = SetTFMapsParams(Nav, Lrep, Lfp.sampleTimes);`

     adjust the paramters if needed. For instance:

     `TFmapsparams.freqrange = [20 150];`

     `TFmapsparams.subset.Condition = [1 3 5];`
   
      `TFmapsparams.subset.Condition_op = 'ismember';`
   
      `TFmapsparams.subset.Spd = 2.5;`
   
      `TFmapsparams.subset.Spd = '>=';`

      `TFmapsparams.Xvariablename = 'Xpos';`
   
      `TFmapsparams.Xbinedges = 0:2:100;`

      `TFmapsparams.Xvariablename = [];`

5. **Perform time-frequency analysis using TFMapsAnalysis**
   
    `TFMaps = TFMapsAnalysis(Nav, Lrep, TFmapsparams, Lfp.sampleTimes);`


For more detailed information and usage examples of each function, please refer to the function documentation, `Tutorial4.1` and `Tutorial_handson`.

Developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
