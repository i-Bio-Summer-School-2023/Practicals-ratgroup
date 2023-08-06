function decparams = DefineDecParams(Nav,Spk)
%Define a set of parameters needed to decode positions from spike trains

%Experimental condition over which the decoder will be trained
decparams.condition = [1 3 5];

%lap type over which the decoder will be trained
decparams.dir = [-1 1];

%lap types over which the decoder will be trained
decparams.laptype = [-1 0 1];

%Minimum speed threshold over which the decoder will be trained
decparams.spdthreshold = 2.5;

%Subset of cells used for decoding. By default we'll use only pyramidal
%cells since interneurons with high firing rates can strongly bias the
%decoder.
decparams.cellidx = Spk.PyrCell;

%Sampling rate of the data
decparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
decparams.scalingFactor = 1 / decparams.sampleRate;

%Name of the X variable to decode. Default is Xpos
decparams.Xvariablename = 'Xpos';%(check here)

%Range of positions (in cm) over which decoding will be performed, ie over
%which place fields will be computed
decparams.Xrange = [0 100];

%Size of the position bins (in cm).
decparams.Xbinsize = 4;

%Size of the gaussian window for smoothing place fields (in cm).
decparams.Xsmthbinsize = 2;

%Size of the gaussian window for smoothing place fields (in bins).
decparams.XsmthNbins = decparams.Xsmthbinsize / decparams.Xbinsize;

%Edges of position bins used to discretize positions
decparams.Xbinedges = decparams.Xrange(1):decparams.Xbinsize:decparams.Xrange(2);

%Name of the Y variable to decode. Default is XDir.
decparams.Yvariablename = 'XDir';%(check here)

%Range of Y over which place fields will be estimated.
decparams.Yrange = [-2 2];%(check here)

%Size of the Y bins.
decparams.Ybinsize = 2;%(check here)

%Size of the gaussian window for smoothing place fields along Y
decparams.Ysmthbinsize = 0;%(check here)

%Size of the gaussian window for smoothing place fields along Y (in bins).
decparams.YsmthNbins = decparams.Ysmthbinsize / decparams.Ybinsize;%(check here)

%Edges of Y bins used to discretize Y
decparams.Ybinedges = decparams.Yrange(1):decparams.Ybinsize:decparams.Yrange(2);%(check here)

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