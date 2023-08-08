function Spk = LoaddataSpk(loadparams, sampleTimes)
% Spk = LoaddataSpk(loadparams, sampleTimes)
%
%Load spiking data into Spk structure, resampled according to timestamps
%provided in sampleTimes and other parameters defined in loadparams. See
%DefineLoadParams.m for a description of those parameters.
%
% INPUT:
% - loadparams:  a structure whose fields contain the parameters necessary
% to load the spiking data data and resample them.
% See DefineLoadParams.m for a description of these parameters.
%
% OUTPUT:
% - Spk: a matlab structure whose fields contain the different types of
% behavioral data resampled at the desired sampling rate (defined in
% loadparams.samplingRate).
%
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
%
% USAGE:
%  Spk = LoaddataSpk(loadparams, sampleTimes)
%
%
% written by J.Fournier 08/2023 for the iBio Summer school

%%
%loading spike times and cluster ID from the prepared .mat file
S = load([loadparams.Datafolder filesep loadparams.spkfilename]);

%Removing spikes that are before or after behavior started
extraspk = S.HippoSpikes(:,1) < sampleTimes(1) | S.HippoSpikes(:,1) > sampleTimes(end);
S.HippoSpikes(extraspk,:) = [];

%Saving spike times and cluster IDs.
Spk.spikeTimes = S.HippoSpikes(:,1);% - Nav.sampleTimes(1);
Spk.spikeID = S.HippoSpikes(:,2);

%convert spike times into an array of spike trains, sampled according to 
%sampleTimes.
clustList = unique(S.HippoSpikes(:,2));
ncells = max(clustList);
nTimeSamples = numel(sampleTimes);
sampleRate = 1 / mean(diff(sampleTimes));

Spk.spikeTrain = zeros(nTimeSamples, ncells);
binEdges = [sampleTimes ; max(sampleTimes) + 1/sampleRate];

for icell = clustList(:)'
    s = S.HippoSpikes(S.HippoSpikes(:,2) == icell,1);
    Spk.spikeTrain(:,icell) = histcounts(s, binEdges);
end
Spk.sampleTimes = sampleTimes;

%Saving some cluster info into the Spk structure
S = load([loadparams.Datafolder filesep loadparams.spkinfofilename]);
HippoClustidx = ismember(S.IndexType(:,3), loadparams.ShankList);
Spk.shankID = S.IndexType(HippoClustidx,3)';
Spk.PyrCell = (S.IndexType(HippoClustidx,6) == 1)';
Spk.IntCell = (S.IndexType(HippoClustidx,6) == 2)';

%%
%loading ripples timestamps and filling in Spk.ripple with 1 when there
%is a ripple and Spk.ripplepeak with 1 for the peak of each ripple
Spk.ripple = zeros(numel(sampleTimes),1);
Spk.ripplepeak = zeros(numel(sampleTimes),1);

ripevents = LoadEvents([loadparams.Datafolder filesep loadparams.ripevtfilename]);
ripstart = ripevents.timestamps(contains(ripevents.description,'start'));
ripstop = ripevents.timestamps(contains(ripevents.description,'stop'));
rippeak = ripevents.timestamps(contains(ripevents.description,'peak'));

for idx = 1:numel(ripstart)
    Spk.ripple(sampleTimes >= ripstart(idx) & sampleTimes <= ripstop(idx)) = 1;
    [~, rippeakidx] = min(abs(sampleTimes - rippeak(idx)));
    Spk.ripplepeak(rippeakidx) = 1;
end

Spk.rippleTimes = rippeak;

end