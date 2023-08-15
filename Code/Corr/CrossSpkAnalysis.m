function Cross = CrossSpkAnalysis(Nav, Srep, crossparams)
% CrossSpkAnalysis - Estimates cross-correlations between responses provided in columns of Srep. 
% Signal correlations are estimated by shuffling time points within bins of
% the variables defined in crossparams.variablenames. Trial-by-trial (or
% noise) correlations are then computed as the overall correlation - the
% signal correlation (averaged across all shuffles).
%
% INPUTS:
%   - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables
%   - Srep: Spike response matrix where each column represents the spike 
%   counts of a neuron.
%   - crossparams: Structure containing cross-spike correlation analysis 
%   parameters (output from DefineCrossSpkParams).
%
% OUTPUTS:
%   Cross: Structure containing cross-correlation analysis results.
%
%   Cross has the following fields:
%   - crossparams: Parameters used for cross-correlation analysis.
%   - lagbins: Time bins for the cross-correlation lag.
%   - ccAll: Pair-wise cross-correlation of original signals.
%   - ccNoise: Pair-wise noise / trial-by-trial correlations (ccAll - ccSig).
%   - ccSig: Pair-wise cross-correlation expected from shared selectivity 
%     to variables indicated in crossparams.variablenames. Signal 
%     correlations are estimated by shuffling time points wihtin bins of
%     the explanatory variables
%   - ccSigSD: Standard deviation of signal cross-correlation estimated by
%   shuffling within bins of the explanatory variables
%   - pval: P-value matrix for the maximum cross-correlation peak.
%   - bestcc: Maximum cross-correlation value.
%   - bestlag: Lag corresponding to the maximum cross-correlation value.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    crossparams = SetCrossSpkParams(Nav, Srep);
%    %change parameters in crossparams here if needed. For instance:
%    %crossparams.lag = 0.5;
%    Cross = CrossSpkAnalysis(Nav, Srep, crossparams)
%
% Written by J Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"

%%
%Selecting time indices over which correlations will be estimated, 
% according to parameters defined in crossparams.subset.
tidx = true(size(Nav.sampleTimes));
pnames = fieldnames(crossparams.subset);
fnames = fieldnames(Nav);
for i = 1:numel(pnames)
    if ismember(pnames{i},fnames)
        fn = str2func(crossparams.subset.([pnames{i} '_op']));
        tidx = tidx & fn(Nav.(pnames{i}), crossparams.subset.(pnames{i}));
    elseif ~strcmp(pnames{i}(end-2:end),'_op')
        warning('some fields of crossparams.subset are not matching fields of Nav')
    end
end

%Selecting cell indices for which pair-wise correlations will be estimated
if islogical(crossparams.cellidx)
    cellidx = find(crossparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > crossparams.nspk_th);
else
    cellidx = crossparams.cellidx(sum(Srep(tidx,crossparams.cellidx), 1, 'omitnan') > crossparams.nspk_th);
end

%Subsetting Srep across cells
spikeTrain = Srep(:,cellidx);

%Sampling rate
sampleRate = 1 / mean(diff(Nav.sampleTimes));

%number of cells selected for place field analysis
ncells = size(spikeTrain, 2);

%%
%Computing pair-wise cross-correlations.

%Smoothing spike trains over a time window. 
spkCount = zeros(size(spikeTrain));
decbinwin = 2 * floor(0.5 * crossparams.timewin * crossparams.sampleRate) + 1;
for icell = 1:size(spikeTrain,2)
    spkCount(:,icell) = smooth(spikeTrain(:,icell), decbinwin, 'moving') * decbinwin;
end

%List of time indices to do the triggered average
idxwin = -round(crossparams.lag * sampleRate):round(crossparams.lag * sampleRate);
lagbins = idxwin / sampleRate;

%Initializing the cross-correlation matrix between cell pairs
ccAll = NaN(ncells, ncells, numel(idxwin));

%Computing the cross-correlation between cell pairs.
%Since spike trains are discrete it is actually much more efficient to
%compute the cross-correlation from the spike-triggered average than using
%xcorr.
parfor icell1 = 1:ncells
    ccAlltemp = NaN(ncells, numel(idxwin));
    for icell2 = icell1+1:ncells
        %Spike indices that are in the right range of time indices
        st1 = find(tidx & spkCount(:,icell1) > 0);

        %Removing values from the spike counts that are not in the
        %right range of time indices
        sp1 = spkCount(:,icell1);
        sp2 = spkCount(:,icell2);
        sp1(~tidx) = NaN;
        sp2(~tidx) = NaN;

        %Getting the snipets of cell2's spikes around cell1's spikes
        [~, ~, l] = ComputeTriggeredAverage(sp2, st1, idxwin, spkCount(st1,icell1));

        %Unnormalized cross-correlation
        c12 = sum(l, 1, 'omitnan');

        %Auto-correlations at zero lag
        c1 = sum(sp1.^2, 'omitnan');
        c2= sum(sp2.^2, 'omitnan');

        %Normalized cross-correlation
        ccAlltemp(icell2,:) = c12 / sqrt(c1 * c2);
    end
    ccAll(icell1,:,:) = ccAlltemp;
end

%filling in the symetrical lower part of the correlation matrix
for l = 1:size(ccAll,3)
    cslice = ccAll(:,:,l);
    csliceT = transpose(cslice);
    cslice(tril(ones(size(cslice)))>0) = csliceT(tril(ones(size(cslice)))>0);
    ccAll(:,:,l) = cslice;
end
%%
%Computing pair-wise cross-correlations after shuffling spikes within
%bins of variables indicated in 

%The shuffling procedure consists in establishing a distribution of
%eigenvalues obtained after shuffling time points within bins of the
%variables provided in crossparams.variablenames.
%We first start by discretizing these variables.
nVars = numel(crossparams.variablenames);
if nVars > 0
    vars_discrete = cell(1,nVars);
    sz = cellfun(@numel,crossparams.binedges) - 1;
    for i = 1:nVars
        vars_discrete{i} = discretize(Nav.(crossparams.variablenames{i}), crossparams.binedges{i});
    end
    %linearizing the indices across all variables, so each bin in the
    %indexed space corresponds to a specific combination of bins in the
    %original variable space
    varlin_discrete = sub2ind(sz, vars_discrete{:});
    nbins = prod(sz);
else
    varlin_discrete = ones(ntimepts, 1);
    nbins = 1;
end

%Initializing the cross-correlation matrix for the shuffle controls
ccSigShf = NaN(ncells, ncells, numel(idxwin), crossparams.nShuffle);

%Initializing the random number generator for reproducibility purposes
s = RandStream('mt19937ar','Seed',0);

%Computing cross-correlation after shuffling spikes within position and
%speed bins.
for ishf = 1:crossparams.nShuffle
    %Shuffling spike trains
    spikeTrainShf = NaN(size(spikeTrain));
    for k = 1:nbins
        idx = find(varlin_discrete == k);
        for icell = 1:ncells
            idxshf = idx(randperm(s,numel(idx)));
            spikeTrainShf(idx,icell) = spikeTrain(idxshf,icell);
        end
    end

    %Smoothing spike trains over a time window.
    spkCountShf = zeros(size(spikeTrainShf));
    decbinwin = 2 * floor(0.5 * crossparams.timewin * crossparams.sampleRate) + 1;
    for icell = 1:size(spikeTrainShf,2)
        spkCountShf(:,icell) = smooth(spikeTrainShf(:,icell), decbinwin, 'moving') * decbinwin;
    end
    
    parfor icell1 = 1:ncells
        ccSigShftemp = NaN(ncells, numel(idxwin));
        for icell2 = icell1+1:ncells
            %Spike indices that are in the right range of time indices
            st1 = find(tidx & spkCountShf(:,icell1) > 0);

            %Removing values from the spike counts that are not in the
            %right range of time indices
            sp1 = spkCountShf(:,icell1);
            sp2 = spkCountShf(:,icell2);
            sp1(~tidx) = NaN;
            sp2(~tidx) = NaN;

            %Getting the snipets of cell2's spikes around cell1's spikes
            [~, ~, l] = ComputeTriggeredAverage(sp2, st1, idxwin, spkCountShf(st1,icell1));

            %Unnormalized cross-correlation
            c12 = sum(l, 1, 'omitnan');

            %Auto-correlations at zero lag
            c1 = sum(sp1.^2, 'omitnan');
            c2= sum(sp2.^2, 'omitnan');

            %Normalized cross-correlation
            ccSigShftemp(icell2,:) = c12 / sqrt(c1 * c2);
        end
        ccSigShf(icell1,:,:,ishf) = ccSigShftemp;
    end

    %filling in the symetrical lower part of the correlation matrix
    for l = 1:size(ccSigShf,3)
        cslice = ccSigShf(:,:,l,ishf);
        csliceT = transpose(cslice);
        cslice(tril(ones(size(cslice)))>0) = csliceT(tril(ones(size(cslice)))>0);
        ccSigShf(:,:,l,ishf) = cslice;
    end
end

%%
%Estimating the noise correlation as the actual correlation minus the
%average of the shuffle controls.
ccSig = mean(ccSigShf, 4, 'omitnan');
ccNoise = ccAll - ccSig;
ccSigSD = std(ccAll - ccSigShf, [], 4);

%Estimating the p-value of the maximum of the cross-correlation.
[m, I] = max(ccAll, [], 3);
mshf = squeeze(max(ccSigShf, [], 3));
pval = sum(abs(m) < abs(mshf), 3) / crossparams.nShuffle;
pval(isnan(m)) = NaN;
bestcc = m;
bestlag = lagbins(I);
bestlag(isnan(m)) = NaN;

%%
%Returning results in the output structure
crossparams.tidx = tidx;
Cross.crossparams = crossparams;
Cross.lagbins = lagbins;

ncells_orig = size(Srep, 2);

nlagbins = size(ccAll, 3);
Cross.ccAll = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.ccNoise = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.ccSig = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.ccSigSD = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.pval = NaN(ncells_orig, ncells_orig);
Cross.bestcc = NaN(ncells_orig, ncells_orig);
Cross.bestlag = NaN(ncells_orig, ncells_orig);

for i = 1:ncells
    for j = 1:ncells
        Cross.ccAll(cellidx(i),cellidx(j),:) = ccAll(i,j,:);
        Cross.ccNoise(cellidx(i),cellidx(j),:) = ccNoise(i,j,:);
        Cross.ccSig(cellidx(i),cellidx(j),:) = ccSig(i,j,:);
        Cross.ccSigSD(cellidx(i),cellidx(j),:) = ccSigSD(i,j,:);
        Cross.pval(cellidx(i),cellidx(j)) = pval(i,j);
        Cross.bestcc(cellidx(i),cellidx(j)) = bestcc(i,j);
        Cross.bestlag(cellidx(i),cellidx(j)) = bestlag(i,j);
    end
end
end