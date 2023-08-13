function Lfp = LoaddataLfp(loadparams, sampleTimes)
%Load Lfp data into a structure. Theta oscillations are resampled according
%to timestamps provided in sampleTimes. Raw Lfp signals are sampled at a
%frequency defined by loadparams.sampleRate_rawLfp.

%%
%Loading Lfp oscillations
lfp = matfile([loadparams.Datafolder filesep loadparams.lfpfilename]);
lfpsampleTimes = lfp.hippoLFPs(:,1);
lfpsampRate = 1 / nanmean(diff(lfpsampleTimes));

%Resampling raw LFPs at the final sampling rate
sampleTimes_orig = lfpsampleTimes;
idxStart = find(lfpsampleTimes > sampleTimes(1), 1, 'first');
idxEnd = find(lfpsampleTimes > sampleTimes(end), 1, 'first');
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

%Resampling Theta according to sampleTimes
sampleTimes_orig = lfpsampleTimes;
sampleTimes_new = sampleTimes;
Lfp.Theta = interp1(sampleTimes_orig, Lfp.Theta, sampleTimes_new, 'linear');
Lfp.ThetaPhase = mod(interp1(sampleTimes_orig, unwrap(Lfp.ThetaPhase / 180 * pi), sampleTimes_new, 'linear') / pi * 180, 360); 
Lfp.ThetaPower = interp1(sampleTimes_orig, Lfp.ThetaPower, sampleTimes_new, 'linear');
Lfp.ThetaFreq = interp1(sampleTimes_orig, Lfp.ThetaFreq, sampleTimes_new, 'linear');
Lfp.sampleTimes_theta = sampleTimes_new;

end