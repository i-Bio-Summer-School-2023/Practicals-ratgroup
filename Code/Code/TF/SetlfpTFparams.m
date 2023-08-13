function TFparams = DefinelfpTFparams(Lfp)
%Define a set of parameters needed to for time frequency analysis of Lfp
%data

%Experimental condition over which time frequency analyses will be
%conducted
TFparams.condition = [1 3 5];

%directions along X over which time frequency analyses will be conducted
TFparams.dir = [-1 1];

%Minimum speed threshold over which time frequency analyses will be 
%conducted
TFparams.spdthreshold = 2.5;

%Frequency range for TF analysis
TFparams.freqrange = [1 200];

%Sampling rate of the raw lfps
TFparams.sampleRate_raw = 1 / nanmean(diff(Lfp.sampleTimes_raw));

%Sampling rate of theta
TFparams.sampleRate_theta = 1 / nanmean(diff(Lfp.sampleTimes_theta));
end