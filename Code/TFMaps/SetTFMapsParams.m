function TFmapsparams = SetTFMapsParams(Nav, Lrep, sampleTimes_Lrep)
% SetTFMapsParams - Define a set of parameters needed for time-frequency
% analysis of continuous data.
%
% INPUTS:
%   - Nav: Structure containing explanatory variables (e.g. Nav.Xpos) and 
%   timestamps of the samples (in Nav.sampleTimes).
%   - Lrep: response array where each column represents the signal from a 
%   specific channel. Typically Lfp data.
%   - sampleTimes_Lrep (optional): Sampling times for signals in Lrep. If not 
%   provided, assumes the same sampling as in Nav.sampleTimes.
%
% OUTPUTS:
%   TFmapsparams: Structure containing parameters for time-frequency 
%                 analysis.
%
%   TFmapsparams contains the following fields:
%   - subset: a structure where field names correspond to the name of the
%           fields in Nav that we want to apply the condition on. Fields of
%           subse define the value of these fields of Nav that will be used
%           to subset the data. For instance, if we want to subset data 
%           corresponding to Nav.Condition = 1 or 3 and Nav.Spd >= 5, 
%           TFmapsparams.subset should be defined as:
%           TFmapsparams.subset.Condition = [1 3];
%           TFmapsparams.subset.Condition_op = 'ismember';
%           TFmapsparams.subset.Spd = 5;
%           TFmapsparams.subset.Spd_op = '>=';
%   - freqrange: Frequency range for time-frequency analysis.
%   - chidx: Subset of signals in Lrep for which time-frequency analysis 
%     will be performed.
%   - sampleRate: Sampling rate of the independent variables in Nav with 
%     which frequency power will be correlated to.
%   - sampleRate_raw: Sampling rate of the Lrep signals.
%   - Spectrogram: If true, compute the wavelet transforms.
%   - Coherence: If true, compute the wavelet coherences.
%   - scalingFactor: Scaling factor on the response data to compute maps 
%     (default is 1 for original scale).
%   - Xvariablename: Name of the independent variable used to map the 
%     response along X (default is 'Xpos').
%   - Xbinedges: Edges of position bins used to discretize X.
%   - XsmthNbins: Size of the Gaussian window for smoothing place fields 
%     along X (in bins).
%   - Yvariablename: Name of the independent variable used to map the 
%     response along Y (default is 'XDir').
%   - Ybinedges: Edges of Y bins used to discretize Y.
%   - YsmthNbins: Size of the Gaussian window for smoothing place fields 
%     along Y (in bins).
%   - occ_th: Occupancy threshold above which positions are included in the
%     place field estimate (in seconds).
%   - nShuffle: Number of shuffle controls to perform for randomization
%     when computing maps
%   - kfold: Number of folds to consider for cross-validation of maps.
%   - parallel: If true, the wavelet transforms will run in parallel. 
%     Otherwise, they'll run serial. Parallel is not always the best 
%     option: if the data are too large, the overhead for sending the data 
%     to the workers will be time consuming.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Lfp = LoaddataLfp(loadparams, Nav.sampleTimes);
%    Lrep = Lfp.Lfp_raw;
%    TFmapsparams = SetTFMapsParams(Nav, Lrep, Lfp.sampleTimes)
%
% Written by J Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological recordings"
%

%%
%Conditions over the fields of Nav over which time frequency analysis will
%be performed
%TFmapsparams.subset should be a structure where fields have names of the 
%fields of Nav to which the condition should apply to.
TFmapsparams.subset = [];

%For instance, for the example data set, we can define the following fields
TFmapsparams.subset.Condition = [1 3 5];
TFmapsparams.subset.Condition_op = 'ismember';

TFmapsparams.subset.XDir = [-1 1];
TFmapsparams.subset.XDir_op = 'ismember';

TFmapsparams.subset.laptype = [-1 0 1];
TFmapsparams.subset.laptype_op = 'ismember';

TFmapsparams.subset.Spd =  2.5;
TFmapsparams.subset.Spd_op = '>=';

TFmapsparams.subset.Xpos =  0;
TFmapsparams.subset.Xpos_op = '>=';

TFmapsparams.subset.Xpos =  100;
TFmapsparams.subset.Xpos_op = '<=';

%Frequency range for time-frequency analysis
TFmapsparams.freqrange = [1 200];

%Subset of signals in Lrep for which time-freuqency analysis will be
%performed
TFmapsparams.chidx = 1:size(Lrep,2);

%Sampling rate of the independent variables in Nav with which frequency
%power will be correlated to.
TFmapsparams.sampleRate = 1 / mean(diff(Nav.sampleTimes), 'omitnan');

%if sampleTimes_Lrep is not provided, we assume the signals in Lrep are
%sampled just as signals in Nav.
if nargin < 3
    sampleTimes_Lrep = Nav.sampleTimes; 
end

%Sampling rate of the signal to analyse in the frequency domain
TFmapsparams.sampleRate_raw = 1 / mean(diff(sampleTimes_Lrep), 'omitnan');

%If true, we'll compute the wavelet transforms
TFmapsparams.Spectrogram = true;

%If true, we'll compute the wavelet coherences
TFmapsparams.Coherence = true;

%The following parameters are for computing 2D maps from the time frequency
%decomposition of signals in Lrep.
%Scaling factor on the response data (default is 1 so that
%Lfp power is returned on the original scale)
TFmapsparams.scalingFactor = 1;

%Name of the independent variable used to map the response along X. Default
%is Xpos
TFmapsparams.Xvariablename = 'Xpos';

%Edges of position bins used to discretize X
TFmapsparams.Xbinedges = 0: 4: 100;

%Size of the gaussian window for smoothing place fields along X (in bins).
TFmapsparams.XsmthNbins = 1;

%Name of the independent variable used to map the response along Y. Default
%is XDir.
TFmapsparams.Yvariablename = 'XDir';

%Size of the gaussian window for smoothing place fields along Y (in bins).
TFmapsparams.YsmthNbins = 0;

%Edges of Y bins used to discretize Y
TFmapsparams.Ybinedges = [-2 0 2];

%Occupancy threshold above which positions are included in the place field
%estimate (in seconds)
TFmapsparams.occ_th = 0;

%Number of shuffle controls to perform for randomization
TFmapsparams.nShuffle = 100;

%Number of folds to consider for cross-validation
TFmapsparams.kfold = 10;

%If true, the wavelet transforms will run in parallel. Otherwise they'll 
%run serial. Parallel is not always the best options: if the data are too large, the
%overhead for sending the data to the workers will be time consuming. 
TFmapsparams.parallel = false;

end