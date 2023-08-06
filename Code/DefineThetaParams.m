function thetaparams = DefineThetaParams(Nav,Spk)
%Define a set of parameters needed for theta related analyses

%Experimental condition over which place fields will be estimated
thetaparams.condition = [1 3 5];

%Directions over which place fields will be estimated
thetaparams.dir = [-1 1];

%lap types over which place fields will be estimated
thetaparams.laptype = [-1 0 1];

%Minimum speed threshold over which place fields will be computed
thetaparams.spdthreshold = 2.5;

%Subset of cells for which place fields will be computed
thetaparams.cellidx = true(1, size(Spk.spikeTrain, 2));

%Sampling rate of the data
thetaparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Range of positions over which place fields will be estimated (in cm).
thetaparams.Xrange = [0 100];

%Size of the position bins (in cm).
thetaparams.Xbinsize = 2;

%Size of the gaussian window for smoothing place fields (in cm).
thetaparams.Xsmthbinsize = 2;

%Size of the gaussian window for smoothing place fields (in bins).
thetaparams.XsmthNbins = thetaparams.Xsmthbinsize / thetaparams.Xbinsize;

%Edges of position bins used to discretize positions
thetaparams.Xbinedges = thetaparams.Xrange(1):thetaparams.Xbinsize:thetaparams.Xrange(2);

%Range of theta phases over which tuning curves will be estimated (in
%degrees).
thetaparams.Phsrange = [0 360];

%Size of the theta phase bins (in degrees).
thetaparams.Phsbinsize = 40;

%Size of the gaussian window for smoothing theta phase (in degrees).
thetaparams.Phssmthbinsize = 40;

%Size of the gaussian window for smoothing phase tuning curves (in bins).
thetaparams.PhssmthNbins = thetaparams.Phssmthbinsize / thetaparams.Phsbinsize;

%Edges of bins used to discretize theta phases
thetaparams.Phsbinedges = thetaparams.Phsrange(1):thetaparams.Phsbinsize:thetaparams.Phsrange(2);

%Occupancy threshold above which positions are included in the place field
%estimate
thetaparams.occ_th = 0;

%Minimal number of spikes to consider a cell
thetaparams.nspk_th = 0;

%Number of shuffle controls to perform for randomization
thetaparams.nShuffle = 100;
end