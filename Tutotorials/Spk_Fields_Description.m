% Fields of Nav are the following:
% - sampleTimes: time stamps of the resampled spike trains
%
% - spikeTimes: a ntimes x 1 array of all spike times
%
% - spikeID: cluster IDs of spikes in spikeTimes
%
% - spikeTrain: a nTimes x nCells array of spike counts in bins centered
% around sample times of Spk.sampleTimes
%
% - shankID: 1 x nCells array of ID of the shank where each cluster was 
% recorded
%
% - PyrCell: 1 x nCells logical array. true if the cluster is a putative
% Pyramidal neuron
%
% - IntCell: 1 x nCells logical array. true if the cluster is a putative
% interneuron
%
% - ripple: ntimes x 1 array with ones wherever there is a ripple.
%
% - ripplepeak: ntimes x 1 array with ones for ripple peaks.
%
% - rippleTimes: timestamps of the detected ripple peaks (in seconds)