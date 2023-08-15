function decparams = SetDecParams(Nav,Srep)
% SetDecParams - Define parameters for decoding positions from spike trains
% using a Bayesian decoding approach
%
%   decparams = SetDecParams(Nav, Srep) defines a set of parameters
%   required for decoding positions from spike trains.
%
% INPUTS:
% - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables to decode
% - Srep: Array of spike count responses (ntimes x ncells) used for 
%   decoding positions.
%
% OUTPUT:
% - decparams: A structure containing parameters for decoding positions.
%
% Fields of decparams are as follows:
%   subset: A structure where field names correspond to the names of the
%           fields in Nav on which conditions will be applied. Fields of
%           subset define the values of these fields of Nav to subset the data.
%           For instance, to subset data corresponding to Nav.Condition = 1
%           or 3 and Nav.Spd >= 5, decparams.subset should be defined as:
%           decparams.subset.Condition = [1 3];
%           decparams.subset.Condition_op = 'ismember';
%           decparams.subset.Spd = 5;
%           decparams.subset.Spd_op = '>=';
%   cellidx: A logical array indicating a subset of cells used for decoding.
%   sampleRate: The sampling rate of the data (in Hz).
%   scalingFactor: A scaling factor applied to the response data, typically
%                  set to 1 / samplingRate to convert spiking data to spikes per second.
%   Xvariablename: The name of the variable used to decode the X-axis position.
%                  Default is 'Xpos'.
%   XsmthNbins: The number of bins used for smoothing tuning curves along the
%               X-axis.
%   Xbinedges: The edges of bins used to discretize the X variable.
%   Yvariablename: The name of the variable used to decode the Y-axis position.
%                  Default is 'XDir'.
%   YsmthNbins: The number of bins used for smoothing tuning curves along the
%               Y-axis.
%   Ybinedges: The edges of bins used to discretize the Y variable.
%   occ_th: An occupancy threshold above which positions are included in the
%           decoding estimate (typically in seconds if scalingFactor is 1 / samplingRate).
%   nspk_th: The minimal number of spikes over the training set to consider a cell for decoding.
%   kfold: The number of folds considered for cross-validation on the training set.
%   dectimewin: The size of the decoding window in seconds.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    decparams = SetDecParams(Nav, Srep);
%
% Written by J Fournier in 08 / 2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%

%Conditions over the fields of Nav over which the decoder will be trained
%decparams.subset should be a structure where fields have names of the 
%fields of Nav to which the condition should apply to.
decparams.subset = [];

%For instance, for the example data set, we define the following fields
decparams.subset.Condition = [1 3 5];
decparams.subset.Condition_op = 'ismember';

decparams.subset.XDir = [-1 1];
decparams.subset.XDir_op = 'ismember';

decparams.subset.laptype = [-1 0 1];
decparams.subset.laptype_op = 'ismember';

decparams.subset.Spd =  2.5;
decparams.subset.Spd_op = '>=';

%Subset of cells used for decoding.
decparams.cellidx = true(1, size(Srep, 2));

%Sampling rate of the data
decparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
decparams.scalingFactor = 1 / decparams.sampleRate;

%Name of the X variable to decode. Default is Xpos
decparams.Xvariablename = 'Xpos';

%Size of the gaussian window for smoothing tuning curves along X (in bins).
decparams.XsmthNbins = 1;

%Edges of bins used to discretize the X variable
decparams.Xbinedges = 0 : 4: 100;

%Name of the Y variable to decode. Default is XDir.
decparams.Yvariablename = 'XDir';

%Size of the gaussian window for smoothing tuning curves along Y (in bins).
decparams.YsmthNbins = 0;

%Edges of Y bins used to discretize Y
decparams.Ybinedges = [-2 0 2];

%Occupancy threshold above which positions are included in the place field
%estimate
decparams.occ_th = 0;

%Minimal number of spikes over the train set to consider a cell for 
%decoding
decparams.nspk_th = 0;

%Number of folds to consider for cross-validation on the train set
decparams.kfold = 10;

%Size of the decoding window in seconds
decparams.dectimewin = .30;

end