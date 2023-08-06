function TFlfp = TimeFreqAnalysis(Nav, Spk, Lfp, TFparams)

%%
%Time indices over which some TF analysis will be performed
tidx = ismember(Nav.Condition, TFparams.condition) &...
       ismember(Nav.XDir, TFparams.dir) &...
       Nav.Spd > TFparams.spdthreshold &...
       ~isnan(Nav.Xpos);
        
%Sampling frequency of raw Lfp signals
rawFs = 1 / mean(diff(Lfp.sampleTimes_raw));

%Sampling frequency of behavioral data
navFs = 1 / mean(diff(Nav.sampleTimes));

%%
%Finding segments of data that are less than 300 seconds apart and merge
%them
brkwin = 2 * floor(300 * navFs /2) + 1;
tidxmerge = smooth(double(tidx), brkwin, 'moving') > 0;

startidx = find(diff(tidxmerge(brkwin:end)) > 0) + floor(brkwin/2);
endidx = find(diff(tidxmerge(1:end - brkwin + 1)) < 0) - floor(brkwin/2);
if startidx(1) > endidx(1)
    startidx = [1 ; startidx];
end
if startidx(end) > endidx(end)
    endidx = [endidx ; numel(Nav.sampleTimes)];
end

%Converting the indices from behavioral indices to raw Lfp indices
startidx_lfp = floor((startidx - 1) / navFs * rawFs + 1);
endidx_lfp = floor((endidx - 1) / navFs * rawFs + 1);
nseg = numel(startidx_lfp);

%%
%Computing wavelet transform of the signal and resample power according to 
%Nav.sampleTimes.

%List of frequencies to resample the wavelet transform
fq_new = logspace(log10(TFparams.freqrange(1)),log10(TFparams.freqrange(2)),96);

%Initializing the final array of time x xfrequency wavelet transform
wthpc = NaN(numel(Nav.sampleTimes), numel(fq_new));

%Calculating the wavelet transform over the identified segments of data
for k = 1:nseg
    %wavelet transform with Morse wavelet
    [wt,f,~,~] = cwt(Lfp.LfpHpc_raw(startidx_lfp(k):endidx_lfp(k)),'morse',rawFs, 'FrequencyLimits', [fq_new(1) fq_new(end)]);
    
    %cwt returns frequencies from high to low. Flipping it up/down to avoid
    %confusions later.
    wt = flipud(wt);
    f = flipud(f);

    %Power of oscillations
    wt = abs(wt);
    
    %Time along the first dimension.
    wt = wt';

    %Saving wavelet transform into wtresp at the resolution
    t_old = Lfp.sampleTimes_raw(startidx_lfp(k):endidx_lfp(k));
    t_new = Nav.sampleTimes(startidx(k):endidx(k));
    wthpc(startidx(k):endidx(k),:) = interp2(f, t_old, wt, fq_new, t_new);
end

%%
%To compute the spatial profile of power across frequencies, we can simply
%call PlaceFieldAnalyses2D.
%Getting default parameters to compute spatial profiles
TFparams.mapsparams = DefineMapsParams(Nav,Spk);

%Replacing default values for time indices selection.
TFparams.mapsparams.condition = TFparams.condition;
TFparams.mapsparams.dir = TFparams.dir;
TFparams.mapsparams.spdthreshold = TFparams.spdthreshold;

%Replacing default values related to subsetting and scaling the response
%variable
TFparams.mapsparams.cellidx = true(1, size(wthpc, 2));
TFparams.mapsparams.nspk_th = -inf;
TFparams.mapsparams.scalingFactor = 1;

%Running place field analysis on the wavelet transform
tfMaps = PlaceFieldAnalyses2D(Nav,wthpc,TFparams.mapsparams);

%%
%Computing wavelet coherence between bla and hpc and resampling it 
%according to Nav.sampleTimes.

%List of frequencies to resample the wavelet coherence
fq_new = logspace(log10(TFparams.freqrange(1)),log10(TFparams.freqrange(2)),96);

%Initializing the final array of time x xfrequency wavelet transform
wcoh = NaN(numel(Nav.sampleTimes), numel(fq_new));

%Calculating the wavelet coherence over the identified segments of data
for k = 1:nseg
    %wavelet coherence with Morlet wavelet
    [wc,~,f] = wcoherence(Lfp.LfpHpc_raw(startidx_lfp(k):endidx_lfp(k)),Lfp.LfpBla_raw(startidx_lfp(k):endidx_lfp(k)),rawFs, 'FrequencyLimits', [fq_new(1) fq_new(end)]);
    
    %wcoherence returns frequencies from high to low. Flipping it up/down to avoid
    %confusions later.
    wc = flipud(wc);
    f = flipud(f);

    %Power of oscillations
    wc = abs(wc);
    
    %Time along the first dimension.
    wc = wc';

    %flipping

    %Saving wavelet transform into wtresp at the resolution
    t_old = Lfp.sampleTimes_raw(startidx_lfp(k):endidx_lfp(k));
    t_new = Nav.sampleTimes(startidx(k):endidx(k));
    wcoh(startidx(k):endidx(k),:) = interp2(f, t_old, wc, fq_new, t_new);
end

%%
%To compute the spatial profile of power across frequencies, we can simply
%call PlaceFieldAnalyses2D.

%Getting default parameters to compute spatial profiles
TFparams.mapsparams = DefineMapsParams(Nav,Spk);

%Replacing default values for time indices selection.
TFparams.mapsparams.condition = TFparams.condition;
TFparams.mapsparams.dir = TFparams.dir;
TFparams.mapsparams.spdthreshold = TFparams.spdthreshold;

%Replacing default values related to subsetting and scaling the response
%variable
TFparams.mapsparams.cellidx = true(1, size(wcoh, 2));
TFparams.mapsparams.nspk_th = -inf;
TFparams.mapsparams.scalingFactor = 1;

%Running place field analysis on wavelet coherence
cohMaps = PlaceFieldAnalyses2D(Nav,wcoh,TFparams.mapsparams);

%%
%Returning the results into TFlfp
TFparams.tidx = tidx;
TFlfp.TFparams = TFparams;
TFlfp.tfMaps = tfMaps;
TFlfp.cohMaps = cohMaps;
TFlfp.fqbins = fq_new;









%%
%Everything starting from here is optional
%Dependence of Theta oscillations on running speed

%Filter Lfp signal in the theta band
rawFs = 1 / mean(diff(Lfp.sampleTimes_raw));
theta = filterellip(Lfp.LfpHpc_raw, 6, 10, rawFs);  

%Computing the power and the instantaneous frequency of theta oscillations.
h = hilbert(theta);
thetapow = abs(h);
thetaphs = angle(h);
thetafreq = [0 ; diff(unwrap(thetaphs)) * rawFs / (2*pi)];

%Resampling theta power and frequency at the sampling rate of behavioral
%data.
thetapow = interp1(Lfp.sampleTimes_raw, thetapow, Nav.sampleTimes, 'linear');
thetafreq = interp1(Lfp.sampleTimes_raw, thetafreq, Nav.sampleTimes, 'linear');

%Discretizing speeds
spdbinsize = 2;
spdrange = [0 30];
spdbinedges = spdrange(1):spdbinsize:spdrange(2);
spdbincenters = spdbinedges(1:end-1) + spdbinsize / 2;
nspdbins = numel(spdbincenters);
spd = Nav.Spd;
spd(spd > spdbinedges(end)) = spdbinedges(end);
spd_discrete = discretize(spd, spdbinedges);

%Computing the occupancy map across speed bins
flat =  ones(size(spd_discrete));
occmap = Compute1DMap(spd_discrete(tidx), flat(tidx), nspdbins);

%Computing the sum of theta power and theta frequency across speeds
ThetaPowSmap = Compute1DMap(spd_discrete(tidx), thetapow(tidx), nspdbins);
ThetaFqSpdmap = Compute1DMap(spd_discrete(tidx), thetafreq(tidx), nspdbins);

%Smoothing maps with a Gaussian of s.d. of spdbinsize
occmap = GaussianSmooth1D(occmap, spdbinsize);
ThetaPowSmap = GaussianSmooth1D(ThetaPowSmap, spdbinsize);
ThetaFqSpdmap = GaussianSmooth1D(ThetaFqSpdmap, spdbinsize);

%Computing the mean maps by dividing with the occupancy map
ThetaPowSmap = ThetaPowSmap ./ occmap;
ThetaFqSpdmap = ThetaFqSpdmap ./ occmap;


%%
%Time frequency analysis of Hippocampal LFP around ripple peaks
%Resample ripple peak vector to the sampling frequency of the raw Lfps
ripFs = 1 / mean(diff(Spk.sampleTimes));
riptimes = find(Spk.ripplepeak == 1);
riptimes = round((riptimes - 1) / ripFs * rawFs + 1);

%Extract snippets of Lfp around ripple peaks.
idxwin = -round(0.1 * rawFs):round(0.1 * rawFs);
[~, ~, lrip] = ComputeTriggeredAverage(Lfp.LfpHpc_raw, riptimes, idxwin);

%Computing the average spectrogram around ripple time
nrip = numel(riptimes);
sp_rip = 0;
for i = 1:nrip
    [wt,f,~,fb] = cwt(lrip(i,:),'morse', rawFs, 'FrequencyLimits', TFparams.freqrange);
   
    %Taking the power of the spectrogram.
    wt = abs(wt);

    %averaging spectrogram
    sp_rip = sp_rip + wt / nrip;

    % %averaging spectrogram after normalization across frequencies to
    % %discard broadband increase in power.
    % sp_ripnorm = sp_ripnorm + (wt./ mean(wt, 1)) / nrip;
end

%Filter Lfp signal in the ripple band
lfpfilt = filterellip(Lfp.LfpHpc_raw, 100, 200, rawFs);  

%Computing the power in the ripple band
lfppow = abs(hilbert(lfpfilt));

%Computing the triggered average of power in the ripple band around ripple
%times
idxwin = -round(0.05 * rawFs):round(0.05 * rawFs);
[rippow, ~, ~] = ComputeTriggeredAverage(lfppow, riptimes, idxwin);

%Add some pltting functions to see each of the steps

%%
%Spectrogram triggered on Theta trough

%Filtering raw Lfp in the theta band.
theta = filterellip(Lfp.LfpHpc_raw, 6, 9, rawFs);

%Computing the Hilbert transform of the filtered signal to get the phase,
%power and instantaneous frequency of Theta oscillations
hilbertTrans = hilbert(theta);
phs = angle(hilbertTrans);

%Converting phases from [-pi,pi[ radians to [0,360[ degrees
ThetaPhase = mod(phs / pi * 180, 360);

%Resample state vector according to Lfp.sampleTimes_raw
cond = interp1(Nav.sampleTimes, Nav.Condition, Lfp.sampleTimes_raw, 'nearest');

%Resample state vector according to Lfp.sampleTimes_raw
state = interp1(Nav.sampleTimes, Nav.state, Lfp.sampleTimes_raw, 'nearest');

%Time frequency analysis of Hippocampal LFP around ripple peaks
%Resample ripple peak vector to the sampling frequency of the raw Lfps
troughidx = find(diff(sign(ThetaPhase - 180)) > 0 & ~ismember(cond(1:end-1), [1 3 5]));

%Extract snippets of Lfp around ritheta troughs.
idxwin = -round(0.1 * rawFs):round(0.1 * rawFs);
[~, ~, lthrough] = ComputeTriggeredAverage(Lfp.LfpHpc_raw, troughidx, idxwin);

%Removing snippets with missing values
nanidx = sum(isnan(lthrough),2) > 0;
lthrough(nanidx,:) = [];
troughidx(nanidx) = [];

%Also computing the average LFP oscillations around those times
[ltheta, ~, ~] = ComputeTriggeredAverage(theta, troughidx, idxwin);

%Computing the average spectrogram around ripple time
ntrough = min(numel(troughtimes), 10000);
sp_trough = 0;
% sp_ripnorm = 0;
for i = 1:ntrough
    [wt,f,~,fb] = cwt(lthrough(i,:),'morse',rawFs);
   
    %Taking the power of the spectrogram.
    wt = abs(wt);

    %averaging spectrogram
    sp_trough = sp_trough + wt / ntrough;
end

sp_troughZ = (sp_trough - mean(sp_trough, 2))./std(sp_trough, 1, 2);

figure;
subplot(5,1,1);
plot(idxwin / rawFs,ltheta)
set(gca,'Xlim',[-0.1 0.1])
ylabel('Theta')
title('Theta triggered wavelet transform')
c = colorbar;
set(c,'Visible','off')
subplot(5,1,2:5);
imagesc(idxwin / rawFs,f,sp_troughZ)
set(gca,'Yscale','log','Ydir','normal')
set(gca,'Xlim',[-0.1 0.1],'Ylim',[6 200])
set(gca,'Ytick',[5 10 20 40 80 160])
xlabel('time from theta trough (s)')
ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'z-scored power';


%%
%Spike triggered spectrograms

%Selecting neurons for which we'll compute the triggered average
cellidx = find(Spk.PyrCell);
ncells = numel(cellidx);

for icell = 1:ncells
    % %spike times
    % s = Spk.spikeTimes(Spk.spikeID == cellidx(icell));
    % 
    % %Binning spike times according to Lfp.sampleTimes_raw.
    % st = histcounts(s, [Lfp.sampleTimes_raw ; Lfp.sampleTimes_raw(end) + 1/ rawFs])';

    %Converting spike times to indices for the raw Lfps
    sidx = interp1(Lfp.sampleTimes_raw, 1:numel(Lfp.sampleTimes_raw), Spk.spikeTimes(Spk.spikeID == cellidx(icell)), 'nearest', NaN);

    %Extract snippets of Lfp around spike times.
    idxwin = -round(0.1 * rawFs):round(0.1 * rawFs);
    [~, ~, lspk] = ComputeTriggeredAverage(Lfp.LfpHpc_raw, sidx, idxwin);

    %Removing snippets with missing values
    nanidx = sum(isnan(lspk),2) > 0;
    lspk(nanidx,:) = [];
    lspk(nanidx) = [];

    %Computing the average spectrogram around spike times
    nspk = min(numel(st), 1000);
    sp_spk = 0;
    % sp_ripnorm = 0;
    for i = 1:nspk
        [wt,f,~,~] = cwt(lspk(i,:),'morse',rawFs);

        %Taking the power of the spectrogram.
        wt = abs(wt);

        %averaging spectrogram
        sp_spk = sp_spk + wt / ntrough;
    end

    sp_spkZ = (sp_spk - mean(sp_spk, 2))./std(sp_spk, 1, 2);

    %saving those averages
    if icell == 1
        STspec = NaN([ncells size(sp_spk)]);
        STspec_norm = NaN([ncells size(sp_spk)]);
    end
    STspec(icell,:,:) = sp_spk;
    STspec_norm(icell,:,:) = sp_spkZ;
end

figure;
c = colorbar;
set(c,'Visible','off')
subplot(5,1,2:5);
imagesc(idxwin / rawFs,f,sp_spkZ)
set(gca,'Yscale','log','Ydir','normal')
set(gca,'Xlim',[-0.1 0.1],'Ylim',[6 200])
set(gca,'Ytick',[5 10 20 40 80 160])
xlabel('time from spike (s)')
ylabel('Frequency (Hz)')
c = colorbar;
c.Label.String = 'z-scored power';
title('Spike triggered wavelet transform')

%%
%Power spectrums during different states

%Computing the average spectrogram from the wavelet transform of the raw 
%Lfp sampled for each state.
stateList = [-2 -1 0 1];

%list of frequencies to interpolate the mean power spectrums.
f_new = 0.5:0.5:100;
ps_state = NaN(numel(f_new),numel(stateList));

%Threshold on the time window over which to compute wavelet transform since
%matlab online may crash for too long segments of data.
maxTimeWin = 100;

for istate = 1:numel(stateList)
    %timestamps of start and end of the eeg state
    tstart = Nav.sampleTimes(diff(Nav.state == stateList(istate)) > 0);
    tend = Nav.sampleTimes(diff(Nav.state == stateList(istate)) < 0);
    if tstart(1) > tend(1)
        tstart = [Nav.sampleTimes(1) ; tstart];
    end
    if tstart(end) > tend(end)
        tend = [tend ; Nav.sampleTimes(end)];
    end
    n = numel(tstart);
    
    %Averaging power spectrum across all segments corresponding to the
    %current state
    ps_ave = 0;
    for l = 1:n
        rawlfpidx = find(Lfp.sampleTimes_raw >= tstart(l) & Lfp.sampleTimes_raw <= min(tend(l), tstart(l) + maxTimeWin));
        [wt,f] = cwt(Lfp.LfpHpc_raw(rawlfpidx),'morse',rawFs);
        
        %power of the wavelet transform
        wt = abs(wt);

        %averaging the wavelet transform across time
        wt = mean(wt, 2, 'omitnan');

        %Interpolating to the requested set of frequencies, since the
        %resoltuion of the wavelet transform depends on the number of time
        %points
        wt = interp1(f, wt, f_new);

        %Averaging the power spectrum corresponding to that snippet of
        %data
        ps_ave = ps_ave + wt / n;
    end

    %Saving the average power spectrum into ps_state.
    ps_state(:,istate) = ps_ave;
end

sampleTime_state = (Lfp.sampleTimes_raw(1):Lfp.sampleTimes_raw(end))';
%Filtering Lfp signal in the delta band
delta = filterellip(Lfp.LfpHpc_raw, 1, 2, rawFs);
%Computing the power in the delta band
delta = medfilt1(abs(hilbert(delta)), 2 * rawFs, 'moving');
delta = interp1(Lfp.sampleTimes_raw, delta, sampleTime_state);

%Filtering Lfp signal in 
% the theta band
theta = filterellip(Lfp.LfpHpc_raw, 6, 9, rawFs);
%Computing the power in the delta band
theta = medfilt1(abs(hilbert(theta)), 10 * rawFs);
theta = interp1(Lfp.sampleTimes_raw, theta, sampleTime_state);

%Filtering Lfp signal in the low gamma band
lowgamma = filterellip(Lfp.LfpHpc_raw, 30, 60, rawFs);
%Computing the power in the delta band
lowgamma = smooth(abs(hilbert(lowgamma)), 10 * rawFs, 'moving');
lowgamma = interp1(Lfp.sampleTimes_raw, lowgamma, sampleTime_state);

%Resample state vector according to sampleTime_state
state = interp1(Nav.sampleTimes, Nav.state, sampleTime_state, 'nearest');

%Add some pltting functions to look at spectrograms for each state, power
%spectrums and thet/delta ratio vs accelerometer data
end