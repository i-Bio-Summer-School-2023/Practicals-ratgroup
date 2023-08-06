function Dec = DecodingAnalysis(Nav, Srep, decparams)
%Decoding spatial positions from a set of neuron's spike trains. 

%%

if isempty(decparams.Xvariablename)
    decparams.Xvariablename = 'Xpos';
end

%Smoothing the spike train over the decoding window to get spike counts
spkCount = zeros(size(Srep));
decbinwin = 2 * floor(0.5 * decparams.dectimewin * decparams.sampleRate) + 1;
for icell = 1:size(Srep,2)
    spkCount(:,icell) = smooth(Srep(:,icell), decbinwin) * decbinwin;
end

%Time indices over which place fields used for decoding will be estimated.
traintidx = ismember(Nav.Condition, decparams.condition) &...
            ismember(Nav.XDir, decparams.dir) &...
            Nav.Spd > decparams.spdthreshold &...
            ~isnan(Nav.(decparams.Xvariablename));
               
%Selecting cell indices over which to compute place fields
if islogical(decparams.cellidx)
    cellidx = find(decparams.cellidx(:)' & sum(Srep(traintidx,:), 1, 'omitnan') > decparams.nspk_th);
else
    cellidx = decparams.cellidx(sum(Srep(traintidx,decparams.cellidx), 1, 'omitnan') > decparams.nspk_th);
end

%Subsetting spike trains across cells.
spikeTrain = Srep(:,cellidx);
spkCount = spkCount(:,cellidx);

X = Nav.X;

%number of cells selected for decoding
ncells = size(spkCount, 2);

%number of position bins
nbins = numel(decparams.Xbinedges) - 1;

%%
%Discretizing X vectors according to decparams.Xbinedges
X_discrete = discretize(X, decparams.Xbinedges);

%%
X_discrete_trainset = X_discrete(traintidx);
spkTrain_trainset = spikeTrain(traintidx,:);
%Computing occupancy map (same for all cells) on the train set
flat = decparams.scalingFactor * ones(size(X_discrete_trainset));
occmap = Compute1DMap(X_discrete_trainset, flat, nbins);

%Removing occupancy for positions below the occupancy threshold
occmap(occmap <= decparams.occ_th) = NaN;

%Smoothing the occupancy map with a gaussian window (decparams.smthNbins
%of sd)
occmap = GaussianSmooth1D(occmap, decparams.XsmthNbins);


%Computing and smoothing spike count map for each cell
scmap = NaN(ncells, nbins);
for icell = 1:ncells
    scmap(icell,:) = Compute1DMap(X_discrete_trainset, spkTrain_trainset(:,icell), nbins);
    scmap(icell,isnan(occmap)) = NaN;
    scmap(icell,:) = GaussianSmooth1D(scmap(icell,:), decparams.XsmthNbins);
end

%Calculating the place field maps by dividing scmap and occmap
mapX = scmap ./ occmap;

%%
%Computing decoded positions for the entire data set.
[DecMax, DecMean] = ComputeBayesMAP(mapX, spkCount, decparams.dectimewin);

%%
%Actually, place fields computed on the train set are over-fitting the data
%so decoding is valid only for data points not included in the train set.

%number of data points
ntimepts = size(spkCount, 1);

%Filling in decoded position for time bins that were not included in the
%train set.
DecMax_full = NaN(ntimepts,1);
DecMean_full = NaN(ntimepts,1);
DecMax_full(~traintidx) = DecMax(~traintidx);
DecMean_full(~traintidx) = DecMean(~traintidx);

%Doing the same thing now with cross-validated data on the train set.
%First defining a partition of the data for k-fold cross-validation. NB: we
%should normally be more careful about the fact that the spike count data
%are actually smoothed over time...
ntimepts_trainset = sum(traintidx);
cv = crossvalPartition(ntimepts_trainset, decparams.kfold);

%Computing the spike train predicted from the place field using k-fold 
%cross-validation
mapX_cv = NaN(ncells, nbins, decparams.kfold);
DecMax_cv = NaN(ntimepts_trainset,1);
DecMean_cv = NaN(ntimepts_trainset,1);
X_discrete_trainset = X_discrete(traintidx);
spkTrain_trainset = spikeTrain(traintidx,:);
spkCount_trainset = spkCount(traintidx,:);
for i = 1:decparams.kfold
    %Subsetting X and spiketrain according to the train set of the
    %current fold
    Xtraining = X_discrete_trainset(cv.trainsets{i});
    Spktraining = spkTrain_trainset(cv.trainsets{i},:);
    
    %Computing occupancy map for the current fold
    flat = decparams.scalingFactor * ones(size(Xtraining));
    occmap_cv = Compute1DMap(Xtraining, flat, nbins);
    occmap_cv(occmap_cv <= decparams.occ_th) = NaN;
    occmap_cv = GaussianSmooth1D(occmap_cv, decparams.XsmthNbins);
    
    %Computing the spike count map and place field of each cell for the
    %current fold
    for icell = 1:ncells
        scmap_cv = Compute1DMap(Xtraining, Spktraining(:,icell), nbins);
        scmap_cv(isnan(occmap_cv)) = NaN;
        scmap_cv = GaussianSmooth1D(scmap_cv, decparams.XsmthNbins);
        mapX_cv(icell,:,i) = scmap_cv ./ occmap_cv;
    end
end

%Now that we've got cross-validated place fields for the train set, we can
%compute decoded positions on the train set using the same k-fold
%partition.
for i = 1:decparams.kfold
    spkCountTest = spkCount_trainset(cv.testsets{i},:);
    [DecMax_cv(cv.testsets{i}), DecMean_cv(cv.testsets{i})] = ComputeBayesMAP(mapX_cv(:,:,i), spkCountTest, decparams.dectimewin);
end

%Filling in cross-validated decoded positions for the train set.
DecMax_full(traintidx) = DecMax_cv;
DecMean_full(traintidx) = DecMean_cv;

%%
%Populate the output structure with results to be saved
decparams.traintidx = traintidx;
Dec.decparams = decparams;

ncells_orig = size(Srep, 2);
Dec.Xbincenters = decparams.Xbinedges(1:end-1) + decparams.Xbinsize / 2;
Dec.nXbins = nbins;
Dec.mapX = NaN(ncells_orig, nbins);
Dec.mapX_cv = NaN(ncells_orig, nbins, decparams.kfold);
Dec.occmap = NaN(1, nbins);

Dec.mapX(cellidx,:) = mapX;
Dec.mapX_cv(cellidx,:,:) = mapX_cv;
Dec.occmap = occmap;

Dec.X = X_discrete;

Dec.XDecMax = DecMax_full;
Dec.XDecMean = DecMean_full;

%Calculating the distribution of decoded variable as a function of the 
%actual variable over the training set (just to be able to quickly check 
%the quality of decoding on that portion of the data).
Dec.Xdecmat = Compute2DMap(Dec.X(traintidx), Dec.XDecMax(traintidx), ones(size(Dec.X(traintidx))), Dec.nXbins, Dec.nXbins);
end