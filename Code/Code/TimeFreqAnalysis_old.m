function TFlfp = TimeFreqAnalysis(Nav, Spk, Lfp, TFparams)

%%
%Time indices over which some TF analysis will be performed
tidx = ismember(Nav.Condition, TFparams.condition) &...
       ismember(Nav.XDir, TFparams.dir) &...
       Nav.Spd > TFparams.spdthreshold &...
       ~isnan(Nav.Xpos);
        
%Sampling frequency
rawFs = 1 / mean(diff(Lfp.sampleTimes_raw));

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
% sp_ripnorm = 0;
for i = 1:nrip
    [p,f,~,fb] = cwt(lrip(i,:),'morse',rawFs);
   
    %Taking the power of the spectrogram.
    p = abs(p);

    %averaging spectrogram
    sp_rip = sp_rip + p / nrip;

    % %averaging spectrogram after normalization across frequencies to
    % %discard broadband increase in power.
    % sp_ripnorm = sp_ripnorm + (p./ mean(p, 1)) / nrip;
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
    [p,f,~,fb] = cwt(lthrough(i,:),'morse',rawFs);
   
    %Taking the power of the spectrogram.
    p = abs(p);

    %averaging spectrogram
    sp_trough = sp_trough + p / ntrough;
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
    %spike times
    s = Spk.spikeTimes(Spk.spikeID == cellidx(icell));
    
    %Binning spike times according to Lfp.sampleTimes_raw.
    st = histcounts(s, [Lfp.sampleTimes_raw ; Lfp.sampleTimes_raw(end) + 1/ rawFs])';

    %Converting spike times to indices for the raw Lfps
    sidx = find(st > 0 & ~ismember(cond, [1 3 5]));

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
        [p,f,~,~] = cwt(lspk(i,:),'morse',rawFs);

        %Taking the power of the spectrogram.
        p = abs(p);

        %averaging spectrogram
        sp_spk = sp_spk + p / ntrough;
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
    tstart = Nav.sampleTimes(sign(diff(Nav.state == stateList(istate))) > 0);
    tend = Nav.sampleTimes(sign(diff(Nav.state == stateList(istate))) < 0);
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
        [p,f] = cwt(Lfp.LfpHpc_raw(rawlfpidx),'morse',rawFs);
        
        %power of the wavelet transform
        p = abs(p);

        %averaging the wavelet transform across time
        p = mean(p, 2, 'omitnan');

        %Interpolating to the requested set of frequencies, since the
        %resoltuion of the wavelet transform depends on the number of time
        %points
        p = interp1(f, p, f_new);

        %Averaging the power spectrum corresponding to that snippet of
        %data
        ps_ave = ps_ave + p / n;
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