function Phsmapsparams = SetPhsMapsParams(Nav, Srep)
%Define a set of parameters needed for phase coupling between Lrep and
%Srep.

%Conditions over the fields of Nav over which phase coupling analysis will
%be performed.
%Phsmapsparams.subset should be a structure where fields have names of the 
%fields of Nav to which the condition should apply to.
Phsmapsparams.subset = [];

%For instance, for the example data set, we can define the following fields
Phsmapsparams.subset.Condition = [1 3 5];
Phsmapsparams.subset.Condition_op = 'ismember';

Phsmapsparams.subset.XDir = [-1 1];
Phsmapsparams.subset.XDir_op = 'ismember';

Phsmapsparams.subset.laptype = [-1 0 1];
Phsmapsparams.subset.laptype_op = 'ismember';

Phsmapsparams.subset.Spd =  2.5;
Phsmapsparams.subset.Spd_op = '>=';

Phsmapsparams.subset.Xpos =  0;
Phsmapsparams.subset.Xpos_op = '>=';

Phsmapsparams.subset.Xpos =  100;
Phsmapsparams.subset.Xpos_op = '<=';

%Subset of cells for which phase coupling will be estimated
Phsmapsparams.cellidx = true(1, size(Srep, 2));

%Sampling rate of the data
Phsmapsparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Frequency range over which phase will be estimated
Phsmapsparams.freqrange = [6 10];

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
Phsmapsparams.scalingFactor = 1 / Phsmapsparams.sampleRate;

%Name of the independent variable used to map the response along X. Default
%is Xpos
Phsmapsparams.Xvariablename = 'Xpos';

%Size of the gaussian window for smoothing place fields (in bins).
Phsmapsparams.XsmthNbins = 1;

%Edges of position bins used to discretize positions
Phsmapsparams.Xbinedges = 0 : 4 : 100;

%Size of the gaussian window for smoothing phase tuning curves (in bins).
Phsmapsparams.PhssmthNbins = 1;

%Edges of bins used to discretize theta phases
Phsmapsparams.Phsbinedges = 0 : 40 : 360;

%Occupancy threshold above which positions are included in the map
%estimate.
Phsmapsparams.occ_th = 0;

%Minimal number of spikes to consider a cell
Phsmapsparams.nspk_th = 0;

%Number of shuffle controls to perform for randomization
Phsmapsparams.nShuffle = 100;
end