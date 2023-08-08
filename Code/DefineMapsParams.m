function mapsparams = DefineMapsParams(Nav,Spk)
% mapsparams = DefineMapsParams(Nav,Spk)
%
%Define a set of parameters needed to compute place fields.
%
% INPUTS:
% - Nav: a structure with at least a field called sampleTimes containing
% sample times of the data from which maps will be computed.
% - Spk: a structure with at least a field called spikeTrain containing the
% time sereis of responses that will be mapped.
%
% OUTPUT:
%  - mapsparams: a structure whose fields contain parameters to run
%  MapsAnalysis1D and MapsAnalysis2D.
%
% Fields of mapsparams are the following:
%
% condition: experimental conditions over which place fields will be 
% estimated.
%
% dir: An array representing the directions over which place fields will 
% be estimated.
% 
% laptype: An array defining lap types over which place fields will be 
% estimated.
% 
% spdthreshold: The minimum speed threshold (in units of velocity) over 
% which place fields will be computed.
% 
% cellidx: A logical array indicating a subset of cells for which place 
% fields will be computed.
% 
% sampleRate: The sampling rate of the data (in Hz).
% 
% scalingFactor: A scaling factor applied to the response data, typically 
% set to 1 / samplingRate to convert spiking data to spikes per second.
% 
% Xvariablename: The name of the independent variable used to map the 
% response along the X-axis.
% 
% Xrange: The range of X values over which place fields will be estimated.
% 
% Xbinsize: The size of the X-axis bins.
% 
% Xsmthbinsize: The size of the Gaussian window for smoothing along the 
% X-axis.
% 
% XsmthNbins: The number of bins used for smoothing place fields along the 
% X-axis.
% 
% Xbinedges: The edges of position bins used to discretize the X-axis.
% 
% Yvariablename: The name of the independent variable used to map the 
% response along the Y-axis.
% 
% Yrange: The range of Y values over which place fields will be estimated.
% 
% Ybinsize: The size of the Y-axis bins.
% 
% Ysmthbinsize: The size of the Gaussian window for smoothing place fields
% along the Y-axis.
% 
% YsmthNbins: The number of bins used for smoothing place fields along the 
% Y-axis.
% 
% Ybinedges: The edges of position bins used to discretize the Y-axis.
% 
% occ_th: An occupancy threshold in seconds above which positions are 
% included in the place field estimate.
% 
% nspk_th: The minimal number of spikes required to consider a cell.
% 
% nShuffle: The number of shuffle controls performed for randomization.
% 
% kfold: The number of folds considered for cross-validation.
%
% USAGE:
% mapsparams = DefineMapsParams(Nav,Spk)
%
% written by J.Fournier 08/2023 for the iBio Summer school

%Experimental condition over which place fields will be estimated
mapsparams.condition = [1 3 5];

%Directions over which place fields will be estimated
mapsparams.dir = [-1 1];

%lap types over which place fields will be estimated
mapsparams.laptype = [-1 0 1];

%Minimum speed threshold over which place fields will be computed
mapsparams.spdthreshold = 2.5;

%Subset of cells for which place fields will be computed
mapsparams.cellidx = true(1, size(Spk.spikeTrain, 2));

%Sampling rate of the data
mapsparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
mapsparams.scalingFactor = 1 / mapsparams.sampleRate;

%Name of the independent variable used to map the response along X. Default
%is Xpos
mapsparams.Xvariablename = 'Xpos';

%Range of X over which place fields will be estimated.
mapsparams.Xrange = [0 100];

%Size of the X bins.
mapsparams.Xbinsize = 4;%(check here)

%Size of the gaussian window for smoothing along X.
mapsparams.Xsmthbinsize = 2;

%Size of the gaussian window for smoothing place fields along X (in bins).
mapsparams.XsmthNbins = mapsparams.Xsmthbinsize / mapsparams.Xbinsize;%(check here)

%Edges of position bins used to discretize X
mapsparams.Xbinedges = mapsparams.Xrange(1):mapsparams.Xbinsize:mapsparams.Xrange(2);%(check here)

%Name of the independent variable used to map the response along Y. Default
%is XDir.
mapsparams.Yvariablename = 'XDir';

%Range of Y over which place fields will be estimated.
mapsparams.Yrange = [-2 2];%(check here)

%Size of the Y bins.
mapsparams.Ybinsize = 2;%(check here)

%Size of the gaussian window for smoothing place fields along Y
mapsparams.Ysmthbinsize = 0;%(check here)

%Size of the gaussian window for smoothing place fields along Y (in bins).
mapsparams.YsmthNbins = mapsparams.Ysmthbinsize / mapsparams.Ybinsize;%(check here)

%Edges of Y bins used to discretize Y
mapsparams.Ybinedges = mapsparams.Yrange(1):mapsparams.Ybinsize:mapsparams.Yrange(2);%(check here)

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