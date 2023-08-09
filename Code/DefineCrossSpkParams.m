function crossparams = DefineCrossParams(Nav,Spk)
%Define a set of parameters needed to compute noise correlations, using a 
%shuffling procedure across position and speed bins.

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