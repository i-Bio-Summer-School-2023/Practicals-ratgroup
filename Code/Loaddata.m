function [Nav, Spk, Lfp] = Loaddata(loadparams)
%Load behavior, spike times and LFP data into Nav, Spk and Lfp structures

%%
%loading position data
S = load([loadparams.Datafolder filesep loadparams.posfilename]);

%Converting positions from pixels to cm
Nav.X = loadparams.pix2cm * S.positions(:,2);
Nav.Y = loadparams.pix2cm * S.positions(:,3);
Nav.sampleTimes = S.positions(:,1);

%%
%loading subsessions indices
catevents = LoadEvents([loadparams.Datafolder filesep loadparams.catevtfilename]);
%Filling in the Nav.Condition vector with 1 for preprun, 2 for presleep, 3
%for run, 4 for postsleep, 5 for postrun
Nav.Condition = NaN(size(Nav.X));
for idx = 1:2:numel(catevents.timestamps)
    strmatch = regexp(catevents.description{idx},{'-prerun','-presleep','-run','-postsleep','-postrun'});
    condID = find(~cellfun(@isempty,strmatch));
    condidx = Nav.sampleTimes >= catevents.timestamps(idx) & Nav.sampleTimes <= catevents.timestamps(idx + 1);
    Nav.Condition(condidx) = condID;
end

%%
%Computing Speed and direction of mvt. 
%These need to be computed for each sessions seperately since recordings 
%are discontinous between sessions.
Nav.Spd = NaN(size(Nav.X));
Nav.smthSpd = NaN(size(Nav.X));
Nav.XDir = NaN(size(Nav.X));%+1 for L to R; -1 for R to L; 0 for unclear.
sampleRate = 1 / mean(diff(Nav.sampleTimes), 'omitnan');
smthwin = round(0.5 * sampleRate);%500 ms smoothing window for speed and direction
for icond = 1:5
    %Subsession indices
    idxcond = find(Nav.Condition == icond);
    
    %Running speed
    Xdiff = [diff(Nav.X(idxcond)) ; NaN];
    Ydiff = [diff(Nav.Y(idxcond)) ; NaN];
    Tdiff = [diff(Nav.sampleTimes(idxcond)) ; Nav.sampleTimes(idxcond(end)) - Nav.sampleTimes(idxcond(end - 1))];
    Nav.Spd(idxcond) = sqrt(Xdiff.^2 + Ydiff.^2)./Tdiff;%(Check here)
    Nav.smthSpd(idxcond) = smooth(sqrt(Xdiff.^2 + Ydiff.^2)./Tdiff, smthwin, 'moving');%(Check here)
    
    %Speed along X
    Xspd = smooth(Xdiff ./ Tdiff, smthwin, 'moving');%(Check here)
    
    %Main direction of travel along X
    XDirLtoR = smooth(double(Xspd > 0), smthwin) > 0;%(Check here)
    XDirRtoL = smooth(double(Xspd < 0), smthwin) > 0;%(Check here)
    Nav.XDir(idxcond) = double(XDirLtoR - XDirRtoL > 0) - double(XDirLtoR - XDirRtoL < 0);%(Check here)
end

%%
%Resampling behavioral data to the final resolution (loadparams.sampleRate) 
sampleTimes_orig = Nav.sampleTimes;
sampleTimes_new = (min(sampleTimes_orig):(1/loadparams.sampleRate):max(sampleTimes_orig))';
Nav.X = interp1(sampleTimes_orig, Nav.X, sampleTimes_new, 'linear'); 
Nav.Y = interp1(sampleTimes_orig, Nav.Y, sampleTimes_new, 'linear');
Nav.Spd = interp1(sampleTimes_orig, Nav.Spd, sampleTimes_new, 'linear');
Nav.smthSpd = interp1(sampleTimes_orig, Nav.smthSpd, sampleTimes_new, 'linear');
Nav.XDir = interp1(sampleTimes_orig, Nav.XDir, sampleTimes_new, 'nearest');
Nav.Condition = interp1(sampleTimes_orig, Nav.Condition, sampleTimes_new, 'nearest');
Nav.sampleTimes = sampleTimes_new;


%%
%loading information about types of lap (left to right, right to left). 
%Nav.laptype equals 1 or -1 for left to right and right to left trials
%respectively.
Nav.laptype = zeros(size(Nav.X));
Nav.uturn = zeros(size(Nav.X));
LapType = load([loadparams.Datafolder filesep loadparams.laptypefilename]);
for idx = 1:size(LapType.LtoRlaps, 1)
    Nav.laptype(Nav.sampleTimes >= LapType.LtoRlaps(idx,1) & Nav.sampleTimes <= LapType.LtoRlaps(idx,2)) = 1;
end
for idx = 1:size(LapType.RtoLlaps, 1)
    Nav.laptype(Nav.sampleTimes >= LapType.RtoLlaps(idx,1) & Nav.sampleTimes <= LapType.RtoLlaps(idx,2)) = -1;
end

%Nav.uturn = 1 when the rat makes a uturn before the end of the trial.
%Nav.laptype equals zero when when the rats makes a u-turn before the end.
for idx = 1:size(LapType.Uturnlaps, 1)
    Nav.uturn(Nav.sampleTimes >= LapType.Uturnlaps(idx,1) & Nav.sampleTimes <= LapType.Uturnlaps(idx,2)) = 1;
end

%%
%Expressing positions along X as percentage of the track length
Xtrackstart = min(Nav.X(Nav.laptype ~= 0));
Xtrackend = max(Nav.X(Nav.laptype ~= 0));
Nav.Xpos = 100 * (Nav.X - Xtrackstart) / (Xtrackend - Xtrackstart);
Nav.Xpos(Nav.Xpos < 0 | Nav.Xpos > 100) = NaN;

%%
%loading left/right reward times and filling in Nav.reward with 1 for right
%reward and -1 for left reward
Nav.reward = zeros(size(Nav.X));
rrwevents = LoadEvents([loadparams.Datafolder filesep loadparams.rrwevtfilename]);
lrwevents = LoadEvents([loadparams.Datafolder filesep loadparams.lrwevtfilename]);
for idx = 1:numel(rrwevents.timestamps)
    [~, rewidx] = min(abs(Nav.sampleTimes - rrwevents.timestamps(idx)));
    Nav.reward(rewidx) = 1;
end
for idx = 1:numel(lrwevents.timestamps)
    [~, rewidx] = min(abs(Nav.sampleTimes - lrwevents.timestamps(idx)));
    Nav.reward(rewidx) = -1;
end

%%
%loading air puffs timestamps and filling in Nav.airpuff with 1 when there
%is an air puff
Nav.airpuff = zeros(size(Nav.X));
pufevents = LoadEvents([loadparams.Datafolder filesep loadparams.pufevtfilename]);
for idx = 1:numel(pufevents.timestamps)
    [~, puffidx] = min(abs(Nav.sampleTimes - pufevents.timestamps(idx)));
    Nav.airpuff(puffidx) = 1;
end

%%
%loading eeg states timestamps (wake, rem, nrem) and filling in
%Lfp.eegstate with 1 for wake, 0 for drowsy, -1 for REM sleep and -2 for
%non-REM sleep.
Nav.state = NaN(size(Nav.X));
States = load([loadparams.Datafolder filesep loadparams.statefilename]);
for idx = 1:size(States.wake, 1)
    Nav.state(Nav.sampleTimes >= States.wake(idx,1) & Nav.sampleTimes <= States.wake(idx,2)) = 1;
end
for idx = 1:size(States.drowsy, 1)
    Nav.state(Nav.sampleTimes >= States.drowsy(idx,1) & Nav.sampleTimes <= States.drowsy(idx,2)) = 0;
end
for idx = 1:size(States.sws, 1)
    Nav.state(Nav.sampleTimes >= States.sws(idx,1) & Nav.sampleTimes <= States.sws(idx,2)) = -2;
end
for idx = 1:size(States.Rem, 1)
    Nav.state(Nav.sampleTimes >= States.Rem(idx,1) & Nav.sampleTimes <= States.Rem(idx,2)) = -1;
end

%%
%Loading accelerometer data
acc = matfile([loadparams.Datafolder filesep loadparams.accfilename]);

%Sample times for the accelerometer
accsampleTimes = acc.acc(:,1);

%Resampling the 3 accelerometer channels at the final sampling rate
sampleTimes_orig = accsampleTimes;
sampleTimes_new = Nav.sampleTimes;
Nav.acc = interp1(sampleTimes_orig, acc.acc(:,2:end), sampleTimes_new, 'linear');

%%
%loading spike times and cluster ID from the prepared .mat file
S = load([loadparams.Datafolder filesep loadparams.spkfilename]);

%Removing spikes that are before or after behavior started
extraspk = S.HippoSpikes(:,1) < Nav.sampleTimes(1) | S.HippoSpikes(:,1) > Nav.sampleTimes(end);
S.HippoSpikes(extraspk,:) = [];

%Saving spike times and cluster IDs.
Spk.spikeTimes = S.HippoSpikes(:,1);% - Nav.sampleTimes(1);(Check here)
Spk.spikeID = S.HippoSpikes(:,2);

%convert spike times into an array of spike trains, at the same resolution
%as the behavioral data. The behavioral data and the spike recording have 
%already been aligned together.
clustList = unique(S.HippoSpikes(:,2));
ncells = max(clustList);
nTimeSamples = numel(Nav.sampleTimes);

Spk.spikeTrain = zeros(nTimeSamples, ncells);
binEdges = [Nav.sampleTimes ; max(Nav.sampleTimes) + 1/loadparams.sampleRate];

for icell = clustList(:)'
    s = S.HippoSpikes(S.HippoSpikes(:,2) == icell,1);
    Spk.spikeTrain(:,icell) = histcounts(s, binEdges);
end
Spk.sampleTimes = Nav.sampleTimes;

%Saving some cluster info into the Spk structure
S = load([loadparams.Datafolder filesep loadparams.spkinfofilename]);
HippoClustidx = ismember(S.IndexType(:,3), loadparams.ShankList);
Spk.shankID = S.IndexType(HippoClustidx,3)';
Spk.PyrCell = (S.IndexType(HippoClustidx,6) == 1)';
Spk.IntCell = (S.IndexType(HippoClustidx,6) == 2)';

%%
%loading ripples timestamps and filling in Lfp.ripple with 1 when there
%is a ripple and Lfp.ripplepeak with 1 for the peak of each ripple
Lfp.sampleTimes = Nav.sampleTimes;
Lfp.ripple = zeros(size(Nav.X));
Lfp.ripplepeak = zeros(size(Nav.X));
ripevents = LoadEvents([loadparams.Datafolder filesep loadparams.ripevtfilename]);
ripstart = ripevents.timestamps(contains(ripevents.description,'start'));
ripstop = ripevents.timestamps(contains(ripevents.description,'stop'));
rippeak = ripevents.timestamps(contains(ripevents.description,'peak'));
for idx = 1:numel(ripstart)
    Lfp.ripple(Nav.sampleTimes >= ripstart(idx) & Nav.sampleTimes <= ripstop(idx)) = 1;
    [~, rippeakidx] = min(abs(Nav.sampleTimes - rippeak(idx)));
    Lfp.ripplepeak(rippeakidx) = 1;
end

%%
%Loading Lfp oscillations
lfp = matfile([loadparams.Datafolder filesep loadparams.lfpfilename]);
lfpsampleTimes = lfp.hippoLFPs(:,1);
lfpsampRate = 1 / nanmean(diff(lfpsampleTimes));

%Resampling raw LFPs at the final sampling rate
sampleTimes_orig = lfpsampleTimes;
idxStart = find(lfpsampleTimes > Nav.sampleTimes(1), 1, 'first');
idxEnd = find(lfpsampleTimes > Nav.sampleTimes(end), 1, 'first');
sampleTimes_new = (lfpsampleTimes(idxStart):(1/loadparams.sampleRate_rawLfp):lfpsampleTimes(idxEnd))';
Lfp.LfpHpc_raw = interp1(sampleTimes_orig, lfp.hippoLFPs(:,1 + loadparams.LfpChannel_Hpc), sampleTimes_new, 'linear');
Lfp.LfpBla_raw = interp1(sampleTimes_orig, lfp.blaLFPs(:,1 + loadparams.LfpChannel_Bla), sampleTimes_new, 'linear');
Lfp.sampleTimes_raw = sampleTimes_new;

%Filtering in the theta band
Lfp.Theta = filterellip(lfp.hippoLFPs(:,1 + loadparams.LfpChannel_Hpc), loadparams.ThetaBand(1), loadparams.ThetaBand(2), lfpsampRate);

%Computing Hilbert transform of the filtered signal to get the phasepower
%and instantaneous frequency of Theta oscillations
hilbertTrans = hilbert(Lfp.Theta);
phs = angle(hilbertTrans);
pow = abs(hilbertTrans);

Lfp.ThetaPhase = mod(phs / pi * 180, 360);
Lfp.ThetaPower = pow;
Lfp.ThetaFreq = [0 ; diff(unwrap(phs)) * lfpsampRate / (2*pi)];

%Resampling Theta at the final sampling rate
sampleTimes_orig = lfpsampleTimes;
sampleTimes_new = Nav.sampleTimes;
Lfp.Theta = interp1(sampleTimes_orig, Lfp.Theta, sampleTimes_new, 'linear');
Lfp.ThetaPhase = mod(interp1(sampleTimes_orig, unwrap(Lfp.ThetaPhase / 180 * pi), sampleTimes_new, 'linear') / pi * 180, 360); 
Lfp.ThetaPower = interp1(sampleTimes_orig, Lfp.ThetaPower, sampleTimes_new, 'linear');
Lfp.ThetaFreq = interp1(sampleTimes_orig, Lfp.ThetaFreq, sampleTimes_new, 'linear');

%%
%We then need some function to compute and plots behavior. For instance,
%compute speed profile, triggered average of speed on puff delivery, etc
end