function pattparams = DefinePattParams(Nav,Spk)
%Define a set of parameters needed to detect cell assembly patterns

%Experimental condition over which patterns will be estimated
pattparams.condition = [1 3 5];

%lap type over which patterns will be estimated
pattparams.dir = [-1 1];

%Minimum speed threshold over which patterns will be estimated
pattparams.spdthreshold = 2.5;

%Subset of cells used for pattern detection. By default we'll use only 
%pyramidal cells since interneurons with high firing rates can bias the
%covariance matrix.
pattparams.cellidx = Spk.PyrCell;

%Minimal number of spikes over the train set to consider a cell for 
%pattern detection
pattparams.nspk_th = 0;

%Pattern activation threshold to convert activation strength into
%activation "spikes".
pattparams.strength_th = 5;

%Sampling rate of the data
pattparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Size of the spike count window in seconds
pattparams.timewin = .02;
end