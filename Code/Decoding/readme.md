# Bayesian decoding from spike trains 

These functions perform decoding of up to two variables from neural spike train data using a Bayesian approach. 
It includes two main functions: `SetDecParams` defines options for the decoder. `DecodingAnalysis` then takes these options as an input,
together with a structure containing the variables to decode (called `Nav` here) and an array of responses from which variables will be decoded (e.g `Spk.spikeTrain`).
See the Load function section to see how the `Nav` and `Spk` structures should be constructed.

## SetDecParams

The `SetDecParams` function allows you to define parameters required for decoding positions from spike trains. 
It takes as inputs a structure (`Nav`) containing the variables to decode and the response data (`Srep`). 
The output is a structure `decparams` containing options for the decoding analysis.

Parameters included in `decparams`:
- `subset`: conditions and values used to define the subset over which the decoder will be trained.
- `cellidx`: Logical array to select a subset of cells for decoding.
- `sampleRate`: Sampling rate of the data.
- `scalingFactor`: Scaling factor applied to the response data.
- `Xvariablename`, `XsmthNbins`, `Xbinedges`: Parameters for decoding the X-axis position.
- `Yvariablename`, `YsmthNbins`, `Ybinedges`: Parameters for decoding the Y-axis position.
- `occ_th`: Occupancy threshold for decoding inclusion.
- `nspk_th`: Minimal spikes required for decoding.
- `kfold`: Number of folds for cross-validation.
- `dectimewin`: Size of the decoding window.

## DecodingAnalysis

The `DecodingAnalysis` function performs decoding of up to two variables (X and Y positions) from a set of spike trains using 2D tuning curves. 
Given a set of variables (stored in `Nav`), spike train responses (e.g. `Spk.spikeTrain`), and decoding parameters (`decparams`), this function computes various decoding results and analyses.
Although the decoder may be trained on a subset of the data, decoding is performed over the entoire length of the responses. 
To decode over the trainin set, `DecodingAnalysis` uses a k-fold cross-validation approach.

Results included in the `Dec` structure:
- `decparams`: Input parameters for decoding.
- `Xbincenters`, `nXbins`, `Ybincenters`, `nYbins`: Binning information.
- `mapXY`, `mapXY_cv`: Place field maps and cross-validated maps.
- `occmap`, `X`, `Y`: Occupancy map and discretized position values.
- Decoded positions and errors (`XDecMax`, `YDecMax`, `XErrMax`, `YErrMax`, etc.).
- Confusion matrices (`Xdecmat`, `Ydecmat`) estimated over the training set.

## Usage

1. **Load behavioral data and spike data:**

    `datapath = 'path/to/your/data'`

    `loadparams = SetLoadParams(datapath);`

    `Nav = LoaddataNav(loadparams);`

    `Spk = LoaddataSpk(loadparams, Nav.sampleTimes);`

    `Srep = Spk.spikeTrain;`

2. **Define parameters using SetDecParams:**

    `decparams = SetDecParams(Nav, Srep);`

4. **Modify paramters in `decparams` if needed. For instance:**

    `decparams.subset.Condition = [1 3 5];`
   
    `decparams.subset.Condition_op = 'ismember';`
   
    `decparams.subset.Spd = 2.5;`
   
    `decparams.subset.Spd = '>=';`

    `decparams.Xvariablename = 'Xpos';`
   
    `decparams.Xbinedges = 0:2:100;`

    `decparams.Yvariablename = [];`

6. **Perform decoding using DecodingAnalysis**

     `Dec = DecodingAnalysis(Nav, Srep, decparams);`

For more detailed information and usage examples of each function, please refer to the function documentation, `Tutorial3.1` and `Tutorial_handson`.

Developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
