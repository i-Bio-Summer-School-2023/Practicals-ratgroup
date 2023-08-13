function crossparams = SetCrossSpkParams(Nav,Spk)
% SetCrossSpkParams - Define parameters for computing noise correlations using a shuffling procedure across position and speed bins.
%
% Usage:
%   crossparams = SetCrossSpkParams(Nav, Spk)
%
% Inputs:
%   Nav: Structure containing navigation data (timestamps, positions, speeds, etc.).
%   Spk: Structure containing spike train data for each neuron (timestamps, spike counts, etc.).
%
% Outputs:
%   crossparams: Structure containing cross-spike correlation analysis parameters.
%
% Cross-Spike Correlation Analysis Parameters (within the output structure crossparams):
%   subset: a structure where field names correspond to the name of the
%           fields in Nav that we want to apply the condition on. Fields of subset
%           define the value of these fields of Nav that will be used to subset the
%           data. For instance, if we want to subset data corresponding to
%           Nav.Condition = 1 or 3 and Nav.Spd >= 5, crossparams.subset should be
%           defined as:
%           crossparams.subset.Condition = [1 3];
%           crossparams.subset.Condition_op = 'ismember';
%           crossparams.subset.Spd = 5;
%           crossparams.subset.Spd_op = '>=';
%   cellidx: Subset of cells for which correlations will be computed.
%   sampleRate: Sampling rate of the data.
%   scalingFactor: Scaling factor on the response data.
%   lag: Range of lags over which cross-correlations will be computed (in seconds).
%   - nShuffle: number of shuffle control to perform to establish a
%   distribution of correlation values expected from shared selectivity to
%   the variables indicated in crossparams.variablenames.
%   - variablenames: cell array of names of variables in Nav used to build
%   the shuffle controls: shuffling across time will be performed within 
%   bins of these variables. If variablenames is empty, shuffling will be
%   performed by circularly shifting responses by a random amount > 1
%   second.
%   - binedges: cell array of bin edges to discretize the variables indicated
%   in crossparams.variablenames.
%   nspk_th: Minimal number of spikes to consider a cell.
%   nShuffle: Number of shuffle controls for randomization.
%   timewin: Spike count window in seconds.
%
% Written by J. Fournier in 08/2023 for the iBio Summer school.

%Conditions over the fields of Nav for which pair-wise correlations will be 
%computed. crossparams.subset should be a structure where fields have names
%of the fields of Nav to which the condition should apply to.
crossparams.subset = [];

%For instance, for the example data set, we can define the following fields
crossparams.subset.Condition = [1 3 5];
crossparams.subset.Condition_op = 'ismember';

crossparams.subset.XDir = [-1 1];
crossparams.subset.XDir_op = 'ismember';

crossparams.subset.laptype = [-1 0 1];
crossparams.subset.laptype_op = 'ismember';

crossparams.subset.Spd =  2.5;
crossparams.subset.Spd_op = '>=';

crossparams.subset.Xpos =  0;
crossparams.subset.Xpos_op = '>=';

crossparams.subset.Xpos =  100;
crossparams.subset.Xpos_op = '<=';

%Subset of cells for which correlations will be computed
crossparams.cellidx = true(1, size(Spk.spikeTrain, 2));

%Sampling rate of the data
crossparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
crossparams.scalingFactor = 1 / crossparams.sampleRate;

%minimal and maximal lag to compute cross-correlation (in seconds).
crossparams.lag = 0.1;

%Names of varaibles in Nav that will be used to build the shuffle
%distribution of correlation values. Time points will be shuffled wihtin bins of
%these variables. If empty, the shuffling procedure will simply circularly 
%shift the spike counts by a random amount > 1 second.
crossparams.variablenames{1} = 'Xpos';
crossparams.variablenames{2} = 'Spd';
crossparams.variablenames{3} = 'XDir';

%Bin edges to discretize variables indicated in crossparams.variablenames.
crossparams.binedges{1} = 0 : 4: 100;
crossparams.binedges{2} = [0 : 5 : 50 inf];
crossparams.binedges{3} = [-2 0 2];

%Minimal number of spikes to consider a cell
crossparams.nspk_th = 0;

%Number of shuffle controls to perform for randomization
crossparams.nShuffle = 20;%100; takes a long time without paralellization

%Spike count window in seconds.
crossparams.timewin = 0.1;
end