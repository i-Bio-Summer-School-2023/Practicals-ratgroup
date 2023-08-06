function mapsparams = DefineMapsParams(Nav,Spk)
%Define a set of parameters needed to compute place fields

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
mapsparams.Xvariablename = 'Xpos';%(check here)

%Range of X over which place fields will be estimated.
mapsparams.Xrange = [0 100];%(check here)

%Size of the X bins.
mapsparams.Xbinsize = 4;%(check here)

%Size of the gaussian window for smoothing along X.
mapsparams.Xsmthbinsize = 2;%(check here)

%Size of the gaussian window for smoothing place fields along X (in bins).
mapsparams.XsmthNbins = mapsparams.Xsmthbinsize / mapsparams.Xbinsize;%(check here)

%Edges of position bins used to discretize X
mapsparams.Xbinedges = mapsparams.Xrange(1):mapsparams.Xbinsize:mapsparams.Xrange(2);%(check here)

%Name of the independent variable used to map the response along Y. Default
%is XDir.
mapsparams.Yvariablename = 'XDir';%(check here)

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