function crossparams = DefineCrossSpkParams(Nav,Spk)
% DefineCrossSpkParams - Define parameters for computing noise correlations using a shuffling procedure across position and speed bins.
%
% Usage:
%   crossparams = DefineCrossSpkParams(Nav, Spk)
%
% Inputs:
%   Nav: Structure containing navigation data (timestamps, positions, speeds, etc.).
%   Spk: Structure containing spike train data for each neuron (timestamps, spike counts, etc.).
%
% Outputs:
%   crossparams: Structure containing cross-spike correlation analysis parameters.
%
% Cross-Spike Correlation Analysis Parameters (within the output structure crossparams):
%   condition: Experimental conditions for correlation estimation.
%   dir: Lap types (directions) for correlation estimation.
%   laptype: Types of laps for correlation estimation.
%   spdthreshold: Minimum speed threshold for correlation computation.
%   cellidx: Subset of cells for which correlations will be computed.
%   sampleRate: Sampling rate of the data.
%   scalingFactor: Scaling factor on the response data.
%   lag: Range of lags over which cross-correlations will be computed (in seconds).
%   Xrange: Range of positions over which correlations will be estimated (in cm).
%   Xbinsize: Size of the position bins (in cm).
%   Xbinedges: Edges of position bins used for discretization.
%   Spdrange: Range of speeds over which correlations will be estimated (in cm/s).
%   Spdbinsize: Size of the speed bins (in cm/s).
%   Spdbinedges: Edges of speed bins used for discretization.
%   nspk_th: Minimal number of spikes to consider a cell.
%   nShuffle: Number of shuffle controls for randomization.
%   timewin: Spike count window in seconds.
%
% Written by J. Fournier in 08/2023 for the iBio Summer school.

%Experimental condition over which correlations will be estimated
crossparams.condition = [1 3 5];

%Directions over which correlations will be estimated
crossparams.dir = [-1 1];

%lap types over which correlations will be estimated
crossparams.laptype = [-1 0 1];

%Minimum speed threshold over which correlations will be computed
crossparams.spdthreshold = 2.5;

%Subset of cells for which correlations will be computed
crossparams.cellidx = true(1, size(Spk.spikeTrain, 2));

%Sampling rate of the data
crossparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
crossparams.scalingFactor = 1 / crossparams.sampleRate;

%minimal and maximal lag to compute cross-correlation (in seconds).
crossparams.lag = 0.1;

%Range of positions over which correlations will be estimated (in cm).
crossparams.Xrange = [0 100];%

%Size of the position bins (in cm).
crossparams.Xbinsize = 4;

%Edges of position bins used to discretize positions
crossparams.Xbinedges = crossparams.Xrange(1):crossparams.Xbinsize:crossparams.Xrange(2);

%Range of speeds over which correlations will be estimated (in cm).
crossparams.Spdrange = [2.5 52.5];%(check here)

%Size of the speed bins (in cm / s).
crossparams.Spdbinsize = 5;

%Edges of bins used to discretize speeds
crossparams.Spdbinedges = crossparams.Spdrange(1):crossparams.Spdbinsize:crossparams.Spdrange(2);

%Minimal number of spikes to consider a cell
crossparams.nspk_th = 0;

%Number of shuffle controls to perform for randomization
crossparams.nShuffle = 20;%100; takes a long time without paralellization

%Spike count window in seconds.
crossparams.timewin = 0.1;
end