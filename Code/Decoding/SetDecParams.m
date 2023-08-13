function decparams = DefineDecParams(Nav,Spk)
%Define a set of parameters needed to decode positions from spike trains

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