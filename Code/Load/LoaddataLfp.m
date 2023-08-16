function Lfp = LoaddataLfp(loadparams, sampleTimes)
% Lfp = LoaddataLfp(loadparams, sampleTimes)
%
% Load Lfp data into a MATLAB structure. Theta oscillations are resampled according
% to timestamps provided in sampleTimes. Raw Lfp signals are sampled at a
% frequency defined by loadparams.sampleRate_raw. Typically the output
% structure should contain raw Lfp signal sampled at a higher sampling rate
% than the behavioral data. Timestamps for these raw signals should be
% stored in a field called sampleTimes_raw.
%
% INPUT:
% - loadparams: a structure whose fields contain the parameters necessary
%   to load the Lfp data and resample it.
%   See SetLoadParams.m for a description of these parameters.
% - sampleTimes: time stamps of the resampled Lfp data
%
% OUTPUT:
% - Lfp: a MATLAB structure whose fields contain the Lfp data resampled at
%   the desired sampling rate (defined in loadparams.sampleRate) and processed.
%
% Fields of Lfp are the following:
% - Lfp_raw: an nTimes x 2 array containing raw Lfp signals, with Hpc Lfp in
%   the first column and Bla Lfp in the second column.
% - sampleTimes_raw: time stamps of the resampled raw Lfp signals
% - Theta: filtered Lfp signal in the theta band
% - ThetaPhase: phase of the theta oscillations (in degrees)
% - ThetaPower: power of the theta oscillations
% - ThetaFreq: instantaneous frequency of the theta oscillations
% - sampleTimes_theta: time stamps of the resampled theta-related Lfp data
%
% USAGE:
% datadirpath = <path to the directory containing your data>
% loadparams = SetLoadParams(datadirpath);
% Lfp = LoaddataLfp(loadparams, sampleTimes)
%
% See also: SetLoadParams, LoaddataNav, LoaddataSpk
%
% Written by J.Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"


%%
%Loading Lfp oscillations
lfp = matfile([loadparams.Datafolder filesep loadparams.lfpfilename]);
lfpsampleTimes = lfp.hippoLFPs(:,1);
lfpsampRate = 1 / nanmean(diff(lfpsampleTimes));

%Resampling raw LFPs at the final sampling rate
sampleTimes_orig = lfpsampleTimes;
idxStart = find(lfpsampleTimes > sampleTimes(1), 1, 'first');
idxEnd = find(lfpsampleTimes > sampleTimes(end), 1, 'first');
sampleTimes_new = (lfpsampleTimes(idxStart):(1/loadparams.sampleRate_raw):lfpsampleTimes(idxEnd))';

%loading raw Lfps into Lfp_raw, a ntimes x 2 array with Hpc and Bla Lfps in
%its columns.
Lfp.Lfp_raw = NaN(numel(sampleTimes_new), 2);
Lfp.Lfp_raw(:,1) = interp1(sampleTimes_orig, lfp.hippoLFPs(:,1 + loadparams.LfpChannel_Hpc), sampleTimes_new, 'linear');
Lfp.Lfp_raw(:,2) = interp1(sampleTimes_orig, lfp.blaLFPs(:,1 + loadparams.LfpChannel_Bla), sampleTimes_new, 'linear');
Lfp.sampleTimes = sampleTimes_new;

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