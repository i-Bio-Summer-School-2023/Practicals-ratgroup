%Read me for the load functions:
% provides MATLAB functions to streamline loading and preprocessing of
% behavioral and neural data for the dataset example. These functions
% should be customized to your specific datasets. The main logic is that
% one MATLAB structure (called Nav here) should contain all the explanatory
% variables along which you want to investigate the selectivity of neural 
% responses; a second MATLAB structure (called Spk here) should contain the 
% spiking data (both the spike counts resampled at the sampling rate of the
% explanatory variables and the original spike times); a third MATLAB 
% structure shoud contain signals that need to be sampled at a higher
% sampling rate than the explanatory variables, typically Lfp signals for
% which you want to investigate power at frequencies higher than the
% sampling rate of the explanatory variables.
% 
% *SetLoadParams
% Defines parameters for loading and preprocessing data, used by other loading functions.
%
% *LoaddataNav
% Loads behavioral data into a MATLAB structure using user-defined 
% parameters. The resulting structure contains key behavioral metrics for analysis.
% 
% *LoaddataSpk
% Loads spiking data into a MATLAB structure (Spk) with resampled 
% spike trains and neuron information.
% 
% *LoaddataLfp
% Loads local field potential (LFP) data into a MATLAB structure (Lfp) 
% with resampled LFP signals, theta oscillations, and related metrics.
%
% USAGE:
% %Set data directory path
% datadirpath = '/path/to/data_directory';
% 
% % Define load parameters
% loadparams = SetLoadParams(datadirpath);
%
% %Change loadparams if necessary. For instance
% %loadparams.session = xxxx;
%
% % Load behavioral data
% Nav = LoaddataNav(loadparams);
% 
% % Load spiking data
% Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
% 
% % Load LFP data
% Lfp = LoaddataLfp(loadparams, Nav.sampleTimes);
% 
% % Perform analyses and explore the data using the created structures (Nav, Spk, Lfp)

