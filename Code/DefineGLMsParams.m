function glmsparams = DefineGLMsParams(Nav,Spk)
%Define a set of parameters needed to estimate GLMs

%Experimental condition over which GLMs will be estimated
glmsparams.condition = [1 3 5];

%lap type over which GLMs will be estimated
glmsparams.dir = 1;%[-1 1];

%lap types over which place fields will be estimated
glmsparams.laptype = [-1 0 1];

%Minimum speed threshold over which GLMs will be computed
glmsparams.spdthreshold = 2.5;

%Subset of cells for which GLMs will be computed
glmsparams.cellidx = true(1, size(Spk.spikeTrain, 2));

%Sampling rate of the data
glmsparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
glmsparams.scalingFactor = 1 / glmsparams.sampleRate;

%Range of spatial positions over which GLMs will be estimated (in cm).
glmsparams.Xrange = [0 100];

%Size of the position bins (in cm).
glmsparams.Xbinsize = 4;%2;

%Size of the gaussian window for smoothing position tuning curves (in cm).
glmsparams.Xsmthbinsize = 2;

%Size of the gaussian window for smoothing position tuning curves 
%(in bins).
glmsparams.XsmthNbins = glmsparams.Xsmthbinsize / glmsparams.Xbinsize;

%Edges of position bins used to discretize positions
glmsparams.Xbinedges = [glmsparams.Xrange(1):glmsparams.Xbinsize:glmsparams.Xrange(2)];

%Range of speeds over which GLMs will be estimated (in cm/s).
glmsparams.Spdrange = [2.5 52.5];

%Size of the speed bins (in cm/s).
glmsparams.Spdbinsize = 5;

%Size of the gaussian window for smoothing speed tuning curves (in cm/s).
glmsparams.Spdsmthbinsize = 5;

%Size of the gaussian window for smoothing speed tuning curves (in bins).
glmsparams.SpdsmthNbins = glmsparams.Spdsmthbinsize / glmsparams.Spdbinsize;

%Edges of speed bins used to discretize speed
glmsparams.Spdbinedges = glmsparams.Spdrange(1):glmsparams.Spdbinsize:glmsparams.Spdrange(2);

%Occupancy threshold above which predictors are included in the GLM
%estimate.
glmsparams.occ_th = 0;

%Minimal number of spikes to consider a cell for GLM estimation
glmsparams.nspk_th = 10;

%Number of folds to consider for cross-validation
glmsparams.kfold = 10;

%Whether or not the GLM model should include a constant term
glmsparams.intercept = true;

%Type of regularization in glmnet: 1 for lasso, 0 for ridge, elsticnet for
%everything in between
glmsparams.alpha = 1;

%Maximal number of passes over the data before reaching convergence
glmsparams.maxit = 10^3;

%Convergence threshold for coordinate descent, expressed as a fraction of
%the null deviance.
glmsparams.thresh = 10^-3;

%P-value threshold for considering a predictor as significantly
%contributing to the cell response.
glmsparams.pval_th = 0.05;
end