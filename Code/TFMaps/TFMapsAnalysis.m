function TFMaps = TFMapsAnalysis(Nav, Lrep, TFmapsparams, sampleTimes_Lrep)
%TFMaps = TFMapsAnalysis(Nav, Lrep, TFmapsparams, sampleTimes_Lrep)
%
% TFMapsAnalysis - Perform time-frequency analysis and computes maps from
% the power across frequencies as a function of explanatory variables
% defined in TFmapsparams
%
% INPUTS:
%   - Nav: Structure containing explanatory variables (e.g. Nav.Xpos) and 
%    timestamps of the samples (in Nav.sampleTimes).
%   - Lrep: response array where each column represents the signals from 
%    different channels. Typically raw Lfp signals sampled at > 100 Hz.
%   - TFmapsparams: Structure containing time-frequency analysis
%   parameters. See SetTFMapsParams for more details.
%   - sampleTimes_Lrep: Sampling times for the signals in Lrep.
%
% Outputs:
%   TFMaps: Structure containing time-frequency analysis results.
%
%   TFMaps has the following fields:
%   - TFmapsparams: Parameters used for time-frequency analysis.
%   - wtMaps: 1 x nCh cell array of structures containing maps of power across 
%     frequencies as a function of the explanatory variables indicated in 
%     TFmapsparams.Xvariablename and TFmapsparams.Yvariablename. Each
%     structure is the result of MapsAnalysis. See MapsAnalysis for more
%     details.
%   - wcMaps: nCh x nCh cell array of structures containing maps of coherence across 
%     frequencies as a function of the explanatory variables indicated in 
%     TFmapsparams.Xvariablename and TFmapsparams.Yvariablename. Each
%     structure is the result of MapsAnalysis for the coherence bewteen a 
%     specific pair of channels in Lrep. See MapsAnalysis for more
%     details.
%   - fqbins: Frequency bins used for the analysis.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Lfp = LoaddataLfp(loadparams, Nav.sampleTimes);
%    Lrep = Lfp.Lfp_raw;
%    TFmapsparams = SetTFMapsParams(Nav, Lrep, Lfp.sampleTimes);
%    %change parameters in TFmapsparams here if needed. For instance:
%    %TFmapsparams.freqrange = [5 100];
%    TFMaps = TFMapsAnalysis(Nav, Lrep, TFmapsparams, Lfp.sampleTimes);
%
% Written by J Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological recordings"
%
%%

%If no Y variable are indicated in mapsparams, we'll just compute a 1D map
if ~isempty(TFmapsparams.Xvariablename)
    X = Nav.(TFmapsparams.Xvariablename);
else
    X = ones(size(Nav.sampleTimes));
    TFmapsparams.Xbinedges = 1;
    TFmapsparams.XsmthNbins = 0;
end

%If no Y variable are indicated in TFmapsparams, we'll just compute a 1D map
if ~isempty(TFmapsparams.Yvariablename)
    Y = Nav.(TFmapsparams.Yvariablename);
else
    Y = ones(size(Nav.sampleTimes));
    TFmapsparams.Ybinedges = 1;
    TFmapsparams.YsmthNbins = 0;
end

%Selecting time indices over which to perform the time frequency 
%decompodition, according to parameters defined in TFmapsparams.subset
tidx = true(size(X));
pnames = fieldnames(TFmapsparams.subset);
fnames = fieldnames(Nav);
for i = 1:numel(pnames)
    if ismember(pnames{i},fnames)
        fn = str2func(TFmapsparams.subset.([pnames{i} '_op']));
        tidx = tidx & fn(Nav.(pnames{i}), TFmapsparams.subset.(pnames{i}));
    elseif ~strcmp(pnames{i}(end-2:end),'_op')
        warning('some fields of TFmapsparams.subset are not matching fields of Nav')
    end
end

%Selecting time indices where X and Y are in the range of bins
tidx = tidx &...
       X >= TFmapsparams.Xbinedges(1) & X <= TFmapsparams.Xbinedges(end) &...
       Y >= TFmapsparams.Ybinedges(1) & Y <= TFmapsparams.Ybinedges(end) &...
       ~isnan(X) & ~isnan(Y);


%Selecting signals to include in the analysis
chidx = find(TFmapsparams.chidx(:)');

%Subsetting Lrep
Lsignal = Lrep(:,chidx);

%Numbre of channels to analyze
nch = size(Lsignal, 2);
        
%Sampling frequency of raw Lrep signals
rawFs = 1 / mean(diff(sampleTimes_Lrep));

%Sampling frequency of independent variables
navFs = 1 / mean(diff(Nav.sampleTimes));

%We'll run the code in serial mode if TFmapsparams.parallel is set to 
%false. This is because the overhead of transferring the data to the 
%workers is time consuming and it's only worth when the number of channels
%to process is large enough, relative to the size of your data. 
%Compare with or without parallelization on your own data and adjust this 
%parameter accordingly.
if ~TFmapsparams.parallel
  parforArg = 0;
else
  parforArg = Inf;
end

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

%Converting the indices from behavioral indices to raw Lrep indices
startidx_lfp = floor((startidx - 1) / navFs * rawFs + 1);
endidx_lfp = floor((endidx - 1) / navFs * rawFs + 1);
nseg = numel(startidx_lfp);

%%
%Computing wavelet transform of the signal and resample power according to 
%Nav.sampleTimes.

%List of frequencies to resample the wavelet transform
fq_new = logspace(log10(TFmapsparams.freqrange(1)),log10(TFmapsparams.freqrange(2)),96);

%Initializing the final array of time x frequency x channels of wavelet 
%transforms
wti = NaN(numel(Nav.sampleTimes), numel(fq_new), nch);

%Copy of Nav.sampleTimes for parallel loop
sampleTimes_Nav = Nav.sampleTimes;

%Calculating the wavelet transform over the identified segments of data
if TFmapsparams.Spectrogram
    for k = 1:nseg
        idxseg_lfp = startidx_lfp(k):endidx_lfp(k);
        t_old = sampleTimes_Lrep(idxseg_lfp);
        idxseg = startidx(k):endidx(k);
        t_new = sampleTimes_Nav(idxseg);
        wtk = NaN(numel(idxseg), numel(fq_new), nch);
        Lrepk = Lrep(idxseg_lfp,:);
        parfor (ich = 1:nch, parforArg)
            Lrepkch = Lrepk(:,ich);
            %wavelet transform with Morse wavelet
            [wt,f,~,~] = cwt(Lrepkch,'morse',rawFs, 'FrequencyLimits', [fq_new(1) fq_new(end)]);

            %cwt returns frequencies from high to low. Flipping it up/down to avoid
            %confusions later.
            wt = flipud(wt);
            f = flipud(f);

            %Power of oscillations
            wt = abs(wt);

            %Time along the first dimension.
            wt = wt';

            %Saving wavelet transform into wtich, interpolated at
            %sampleTimes_Nav
            wtk(:,:,ich) = interp2(f, t_old, wt, fq_new, t_new);
        end
        wti(idxseg,:,:) = wtk;
    end
end

%%
%Computing the spatial profile of power across frequencies using 
%MapsAnalyses.

%Defining parameters to compute spatial profiles
mapsparams = TFmapsparams;
mapsparams.cellidx = true(1, size(wti, 2));
mapsparams.nspk_th = -inf;

wtMaps = cell(1,nch);

if TFmapsparams.Spectrogram
    %Running place field analysis on the wavelet transform
    for ich = 1:nch
        wtMaps{ich} = MapsAnalysis(Nav, wti(:,:,ich), mapsparams);
    end
end

%%
%Computing wavelet coherence between signals of Lrep and resampling it 
%according to Nav.sampleTimes.

%List of frequencies to resample the wavelet coherence
fq_new = logspace(log10(TFmapsparams.freqrange(1)),log10(TFmapsparams.freqrange(2)),96);

%Initializing the final array as a cell array which will contain the 
% time x xfrequency wavelet coherence for each pair of channels of Lrep.
wci = cell(nch);

%Copy of Nav.sampleTimes for parallel loop
sampleTimes_Nav = Nav.sampleTimes;
if nch > 1 && TFmapsparams.Coherence
    %Calculating the wavelet coherence between pairs of channels over the 
    %identified segments of data
    for k = 1:nseg
        idxseg_lfp = startidx_lfp(k):endidx_lfp(k);
        t_old = sampleTimes_Lrep(idxseg_lfp);
        idxseg = startidx(k):endidx(k);
        t_new = sampleTimes_Nav(idxseg);
        wck = NaN(numel(idxseg), numel(fq_new), nch, nch);
        Lrepk = Lrep(idxseg_lfp,:);
        for ich1 = 1 : (nch -1)
            wckch1 = NaN(numel(idxseg), numel(fq_new), nch);
            Lrepkch1 = Lrepk(:,ich1);
            parfor (ich2 = (ich1 + 1) : nch, parforArg)
                Lrepkch2 = Lrepk(:,ich2);
                %wavelet coherence with Morlet wavelet
                [wc,~,f] = wcoherence(Lrepkch1,Lrepkch2,rawFs, 'FrequencyLimits', [fq_new(1) fq_new(end)]);

                %wcoherence returns frequencies from high to low. Flipping it up/down to avoid
                %confusions later.
                wc = flipud(wc);
                f = flipud(f);

                %Power of oscillations
                wc = abs(wc);

                %Time along the first dimension.
                wc = wc';

                %Saving wavelet transform into wtich, interpolated at
                %sampleTimes_Nav
                wckch1(:,:,ich2) = interp2(f, t_old, wc, fq_new, t_new);
            end
            wck(:,:,ich1,:) = wckch1;
        end
        for ich2 = (ich1 + 1) : nch
            wci{ich1,ich2}(idxseg,:) = wck(:,:,ich1,ich2);
        end
    end
end
%%
%Computing the spatial profile of coherence across frequencies, using 
%MapsAnalyses.

wcMaps = cell(nch);

if nch > 1 && TFmapsparams.Coherence
    %Defining parameters to compute spatial profiles
    mapsparams = TFmapsparams;
    mapsparams.cellidx = true(1, size(wci{1,2}, 2));
    mapsparams.nspk_th = -inf;

    %Running place field analysis on the wavelet coherence
    for ich1 = 1:(nch - 1)
        for ich2 = (ich1 + 1):nch
            wcMaps{ich1, ich2} = MapsAnalysis(Nav, wci{ich1,ich2}, mapsparams);
        end
    end
end

%%
%Returning the results into TFMaps
TFmapsparams.tidx = tidx;
TFMaps.TFmapsparams = TFmapsparams;

nch_orig = size(Lrep, 2);
TFMaps.wtMaps = cell(1, nch_orig);
TFMaps.wcMaps = cell(nch_orig);

for ich1 = 1:nch
    TFMaps.wtMaps(chidx(ich1)) = wtMaps(ich1);
    for ich2 = 1:nch
        TFMaps.wcMaps(chidx(ich1),chidx(ich2)) = wcMaps(ich1,ich2);
    end
end

TFMaps.fqbins = fq_new;
end