function Cross = CrossCorrelationAnalyses(Nav, Srep, crossparams)
%Estimates cross-correlations between responses provided in columns of
%Srep. Noise correlations are computed by shuffling across position and
%speed bins.

%%
%Selecting time and cell indices over which to compute pair-wise
%correlations
tidx = ismember(Nav.Condition, crossparams.condition) &...
       ismember(Nav.XDir, crossparams.dir) &...
       ismember(Nav.laptype, crossparams.laptype) &...
       Nav.Spd >= crossparams.spdthreshold &...
       ~isnan(Nav.Xpos);

if islogical(crossparams.cellidx)
    cellidx = find(crossparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > crossparams.nspk_th);
else
    cellidx = crossparams.cellidx(sum(Srep(tidx,crossparams.cellidx), 1, 'omitnan') > crossparams.nspk_th);
end

%Subsetting Srep across cells
spikeTrain = Srep(:,cellidx);

%
Xpos = Nav.Xpos;
Spd = Nav.Spd;

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
ccSig = NaN(ncells, ncells, numel(idxwin));

%Computing the cross-correlation between cell pairs.
%Since spike trains are discrete it is actually much more efficient to
%compute the cross-correlation from the spike-triggered average than using
%xcorr.
for icell1 = 1:ncells
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
        ccSig(icell1,icell2,:) = c12 / sqrt(c1 * c2);
        ccSig(icell2,icell1,:) = c12 / sqrt(c1 * c2);
    end
end

%%
%Computing pair-wise cross-correlations after shuffling spikes within
%position and speed bins.

%Discretizing position vector according to crossparams.Xbinedges
Xpos_discrete = discretize(Xpos, crossparams.Xbinedges);

%number of position bins
nXbins = numel(crossparams.Xbinedges) - 1;

%Discretizing speed vector according to crossparams.Spdbinedges
%Speeds are first clamped within the speed range
Spd(Spd > crossparams.Spdbinedges(end)) = crossparams.Spdbinedges(end);
Spd(Spd < crossparams.Spdbinedges(1)) = crossparams.Spdbinedges(1);

Spd_discrete = discretize(Spd, crossparams.Spdbinedges);

%number of speed bins
nSpdbins = numel(crossparams.Spdbinedges) - 1;

%Initializing the cross-correlation matrix for the shuffle controls
ccSigShf = NaN(ncells, ncells, numel(idxwin), crossparams.nShuffle);

%Initializing the random number generator for reproducibility purposes
s = RandStream('mt19937ar','Seed',0);

%Computing cross-correlation after shuffling spikes within position and
%speed bins.
for ishf = 1:crossparams.nShuffle
    %Shuffling spike trains
    spikeTrainShf = spikeTrain;
    for ipos = 1:nXbins
        for ispd = 1:nSpdbins
            idx = find(tidx & Xpos_discrete == ipos & Spd_discrete == ispd);
            for icell = 1:ncells
                idxshf = idx(randperm(s,numel(idx)));
                spikeTrainShf(idx,icell) = spikeTrainShf(idxshf,icell);
            end
        end
    end

    %Smoothing spike trains over a time window.
    spkCountShf = zeros(size(spikeTrainShf));
    decbinwin = 2 * floor(0.5 * crossparams.timewin * crossparams.sampleRate) + 1;
    for icell = 1:size(spikeTrainShf,2)
        spkCountShf(:,icell) = smooth(spikeTrainShf(:,icell), decbinwin, 'moving') * decbinwin;
    end

    for icell1 = 1:ncells
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
            ccSigShf(icell1,icell2,:,ishf) = c12 / sqrt(c1 * c2);
            ccSigShf(icell2,icell1,:,ishf) = c12 / sqrt(c1 * c2);
        end
    end
end

%%
%Estimating the noise correlation as the actual correlation minus the
%average of the shuffle controls.
ccShf = mean(ccSigShf, 4, 'omitnan');
ccNoise = ccSig - ccShf;
ccShfSD = std(ccSig - ccSigShf, [], 4);

%Estimating the p-value of the maximum of the cross-correlation.
[m, I] = max(ccSig, [], 3);
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

nlagbins = size(ccSig, 3);
Cross.ccSig = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.ccNoise = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.ccShf = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.ccShfSD = NaN(ncells_orig, ncells_orig, nlagbins);
Cross.pval = NaN(ncells_orig, ncells_orig);
Cross.bestcc = NaN(ncells_orig, ncells_orig);
Cross.bestlag = NaN(ncells_orig, ncells_orig);

for i = 1:ncells
    for j = 1:ncells
        Cross.ccSig(cellidx(i),cellidx(j),:) = ccSig(i,j,:);
        Cross.ccNoise(cellidx(i),cellidx(j),:) = ccNoise(i,j,:);
        Cross.ccShf(cellidx(i),cellidx(j),:) = ccShf(i,j,:);
        Cross.ccShfSD(cellidx(i),cellidx(j),:) = ccShfSD(i,j,:);
        Cross.pval(cellidx(i),cellidx(j)) = pval(i,j);
        Cross.bestcc(cellidx(i),cellidx(j)) = bestcc(i,j);
        Cross.bestlag(cellidx(i),cellidx(j)) = bestlag(i,j);
    end
end
end