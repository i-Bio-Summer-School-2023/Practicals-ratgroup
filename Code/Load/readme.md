# Load functions

This folder contains MATLAB functions for loading and preprocessing neural and behavioral data using the functions `SetLoadParams`, `LoaddataNav`, `LoaddataSpk`, and `LoaddataLfp`. 
These functions should be customized to your specific datasets. They should however keep the main logic which is that:
* `LoaddataNav` should return a MATLAB structure (called Nav here) that contains all the explanatory
variables along which you want to investigate the selectivity of neural responses (for instance the position of the animal or a specific stimulus);
* `LoaddataSpk` should return a second MATLAB structure (called Spk here) that contains the spiking data
(both the spike counts resampled at the sampling rate of the explanatory variables and the original spike times);
* `LoaddataLfp` should return a third MATLAB (called Lfp here) that contain signals that need to be sampled at a higher
sampling rate than the explanatory variables stored in Nav (typically Lfp signals for which you want to investigate
power at frequencies higher than the sampling rate of the explanatory variables).

These returned structures (Nav, Spk and Lfp) should have one field called sampleTimes containing the timestamps of the samples stored in the other fields. The sampleTimes for variables stored in Nav and for the spike trains stored in Spk should be the same.
The name of the other fields can be customized according to your preference or your specific experiment. All signals returned in those three structures should be aligned with one another.

## SetLoadParams

The `SetLoadParams` function defines a set of parameters required to load the data. It takes the path to the directory where the data is stored as an input and outputs
a structure `loadparams` containing various parameters needed for data `LoaddataNav`, `LoaddataSpk` and `LoaddataLfp`.

## LoaddataNav

The `LoaddataNav` function should load all the explanotory variables that you are interested in and return it into
a MATLAB structure (called Nav here). It should contain at least a field called `sampleTimes` containing the timestamps of the samples, 
and some other signals in different fields named according to your preference.
It takes the `loadparams` structure as input and outputs a structure `Nav` containing resampled data, usually behavioral data such as positions, speed, experimental conditions, and more.

## LoaddataSpk

The `LoaddataSpk` function loads spiking data into a MATLAB structure, resampled according to a set of time stamps, usually defined as the same as for the data stored  in Nav.sampleTimes.
It takes the `loadparams` structure and these `sampleTimes` as input and outputs a structure `Spk` containing resampled spike train data, spike times and associated information about neurons.

## LoaddataLfp

The `LoaddataLfp` function loads signals that need to be sampled at a higher sampling rate than the explanatory variables.
These signals are returned into a MATLAB structure which should contain at least a field called `sampleTimes` with the timestamps of the samples. 
It takes the `loadparams` structure and a set of `sampleTimes` as input (for alignement) and outputs a structure `Lfp` containing the resampled data and related information.

Example usage:

`datadirpath = '/path/to/data_directory'; %Set data directory path`

`loadparams = SetLoadParams(datadirpath); % Define load parameters`

`loadparams.session = 20230802; %Change loadparams here if necessary.`

`Nav = LoaddataNav(loadparams); % Load behavioral data`

`Spk = LoaddataSpk(loadparams, Nav.sampleTimes); % Load spiking data`

`Lfp = LoaddataLfp(loadparams, Nav.sampleTimes); % Load LFP data`

Then perform analyses and explore the data using the created structures (Nav, Spk, Lfp).

For more usage examples and more information about each function, please refer to the function documentation, comments in the code and Tutorial1.

These functions were developed by J. Fournier in August 2023 for the Summer school "Advanced computational analysis for behavioral and neurophysiological recordings."
