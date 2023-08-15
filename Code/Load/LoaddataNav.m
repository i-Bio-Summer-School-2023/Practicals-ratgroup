function Nav = LoaddataNav(loadparams)
% Nav = LoaddataNav(loadparams);
%
% Loads the behavioral data into a MATLAB structure, using parameters
% defined in loadparams to find the raw data file and preprocess the data.
% Typically, the output structure should contain a field called sampleTimes
% containing the timestamps of each time sample and a set of fields
% containing data which will be used as independent/explanatory variables 
% along which neural responses are investigated.
%
% INPUT:
% - loadparams: a structure whose fields contain the parameters necessary
% to load the behavioral data data and preprocess them.
% See SetLoadParams.m for a description of these parameters.
%
% OUTPUT:
% - Nav: a MATLAB structure whose fields contain different types of
% behavioral data resampled at the desired sampling rate (defined in
% loadparams.samplingRate).
%
% Fields of Nav are the following:
% * sampleTimes: time stamps of the samples for behavioral data
% * X: positions of the animal on the X axis
% * Y: position of the animal on the Y axis
% * Xpos: position along the X axis in percentage of the track length.
% * XDir: direction of movement along the X axis (+1 for left to right, -1
% for right to left)
% * Spd: speed of movement
% * smthSpd: speed smoothed by a box car window
% * Condition: experimental condition corresponding to different 
%   subsessions (1:preprun, 2:presleep, 3:run, 4:postsleep, 5:postrun)
% * laptype: +1 if the animal went from the left to the right platform;
%   -1 if it went from the right one to the left one;
%   0 if it went back to the same platform before reaching the
%   end of the track
% * uturn: +1 on trial where the animal returned to the same platform; 0
% otherwise
% * trialID: trial number. Trial were defined as any continuous period of
% time where the animal was on the track
% * reward: +1 when a reward was delivered; 0 otherwise.
% * airpuff: +1 when an air puff was delivered; 0 otherwise.
% * state: +1 for awake; 0 for drowsy; -1 for REM sleep; -2 for slow wave
% sleep
% * acc: 3-axis accelerometer data (ntimes x 3 array)
%
% All fields of Nav have time along the first dimension.
%
% USAGE:
% datadirpath = <path to the directory containing your data>
% loadparams = SetLoadParams(datadirpath);
% Nav = LoaddataNav(loadparams);
%
% See also: SetLoadParams, LoadEvents, LoaddataSpk, LoaddataLfp
%
% Written by J.Fournier in 08/2023 for the Summer school "Advanced
% computational analysis for behavioral and neurophysiological recordings"

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

%Interpolating Nav.XDir = 0 values to nearest non-zero value for convenience %(Check here)
Nav.XDir = interp1(Nav.sampleTimes(Nav.XDir ~= 0), Nav.XDir(Nav.XDir ~= 0), Nav.sampleTimes, 'nearest');

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
Nav.Xpos(~ismember(Nav.Condition, [1 3 5])) = NaN;

%% Check here
%Defining some trials indices for whenever the animal is on the track and
%spend more than 0.5 second there.
%Looking for potential start and end of trials
trialduration_th = 0.5;

trialStart = find(sign(diff(~isnan(Nav.Xpos(1:end-1)))) > 0) + 1;
trialEnd = find(sign(diff(~isnan(Nav.Xpos))) < 0);

%if the animal is already on the track at the beginning, we modify
%trialStart accordingly
if trialEnd(1) < trialStart(1)
    trialStart = [1 ; trialStart];
end
%if the recording is stopped while the animal is on the track, we modify 
%trialEnd accordingly.
if numel(trialEnd) < numel(trialStart)
    trialEnd = [trialEnd ; numel(Nav.Xpos)];
end

%Initializing the vector of trialIDs
Nav.trialID = NaN(size(Nav.Xpos));

%Checking that the trials are valid (i.e. longer than 1 s)
trialnum = 0;
for k = 1:numel(trialStart)
    if (Nav.sampleTimes(trialEnd(k))- Nav.sampleTimes(trialStart(k))) > trialduration_th
        trialnum = trialnum + 1;
        Nav.trialID(trialStart(k):trialEnd(k)) = trialnum;
    end
end

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
end