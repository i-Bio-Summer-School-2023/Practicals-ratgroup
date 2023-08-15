function mapsparams = SetMapsParams(Nav,Srep)
% SetMapsParams - Define parameters for computing place fields.
%
%   mapsparams = SetMapsParams(Nav, Srep) defines a set of parameters
%   required for computing place fields.
%
% INPUTS:
% - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables.
% - Srep: Array of responses (ntimes x ncells) from which place fields will be
%   estimated.
%
% OUTPUT:
% - mapsparams: A structure containing parameters for running MapsAnalyses.
%
% Fields of mapsparams are as follows:
%
%   subset: A structure where field names correspond to the names of the
%           fields in Nav on which conditions will be applied. Fields of
%           subset define the values of these fields of Nav to subset the data.
%           For instance, to subset data corresponding to Nav.Condition = 1
%           or 3 and Nav.Spd >= 5, mapsparams.subset should be defined as:
%           mapsparams.subset.Condition = [1 3];
%           mapsparams.subset.Condition_op = 'ismember';
%           mapsparams.subset.Spd = 5;
%           mapsparams.subset.Spd_op = '>=';
%
%   cellidx: A logical array indicating a subset of cells for which place
%            fields will be computed.
%   sampleRate: The sampling rate of the data (in Hz).
%   scalingFactor: A scaling factor applied to the response data, typically
%                  set to 1 / samplingRate to convert spiking data to spikes per second.
%   Xvariablename: The name of the independent variable used to map the
%                  response along the X-axis.
%   XsmthNbins: The number of bins used for smoothing place fields along the
%               X-axis.
%   Xbinedges: The edges of position bins used to discretize the X-axis.
%   Yvariablename: The name of the independent variable used to map the
%                  response along the Y-axis.
%   YsmthNbins: The number of bins used for smoothing place fields along the
%               Y-axis.
%   Ybinedges: The edges of position bins used to discretize the Y-axis.
%   occ_th: An occupancy threshold above which positions are
%           included in the place field estimate (typically in seconds if 
%           scalingFactor is 1 / samplingRate)
%   nspk_th: The minimum number of spikes required to consider a cell.
%   nShuffle: The number of shuffle controls performed for randomization.
%   kfold: The number of folds considered for cross-validation.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    mapsparams = SetMapsParams(Nav, Srep);
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%

%Conditions over the fields of Nav for which place fields will be estimated
%mapsparams.subset should be a structure where fields have names of the 
%fields of Nav to which the condition should apply to.
mapsparams.subset = [];

%For instance, for the example data set, we define the following fields
mapsparams.subset.Condition = [1 3 5];
mapsparams.subset.Condition_op = 'ismember';

mapsparams.subset.XDir = [-1 1];
mapsparams.subset.XDir_op = 'ismember';

mapsparams.subset.laptype = [-1 0 1];
mapsparams.subset.laptype_op = 'ismember';

mapsparams.subset.Spd =  2.5;
mapsparams.subset.Spd_op = '>=';

%Subset of cells for which place fields will be computed
mapsparams.cellidx = true(1, size(Srep, 2));

%Sampling rate of the data
mapsparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
mapsparams.scalingFactor = 1 / mapsparams.sampleRate;

%Name of the independent variable used to map the response along X. Default
%is Xpos
mapsparams.Xvariablename = 'Xpos';

%Edges of position bins used to discretize X
mapsparams.Xbinedges = 0: 4: 100;

%Size of the gaussian window for smoothing place fields along X (in bins).
mapsparams.XsmthNbins = 1;

%Name of the independent variable used to map the response along Y. Default
%is XDir.
mapsparams.Yvariablename = 'XDir';

%Size of the gaussian window for smoothing place fields along Y (in bins).
mapsparams.YsmthNbins = 0;

%Edges of Y bins used to discretize Y
mapsparams.Ybinedges = [-2 0 2];

%Occupancy threshold above which positions are included in the place field
%estimate (in seconds)
mapsparams.occ_th = 0;

%Minimal number of spikes to consider a cell
mapsparams.nspk_th = 0;

%Number of shuffle controls to perform for randomization
mapsparams.nShuffle = 100;

%Number of folds to consider for cross-validation
mapsparams.kfold = 10;

end