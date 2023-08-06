function Dec = DecodingAnalysis2D(Nav, Srep, decparams)
%Decoding spatial positions x direction of travel from a set of neuron's 
%spike trains. Same as DecodingAnalysis but based on 2D spatial maps of
%size nXbins x 2 (ie place fields corresponding to L to R and R to L
%directions).
%%

if isempty(decparams.Xvariablename)
    decparams.Xvariablename = 'X';
end
if isempty(decparams.Yvariablename)
    decparams.Yvariablename = 'XDir';
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
            ~isnan(Nav.(decparams.Xvariablename)) &...
            ~isnan(Nav.(decparams.Yvariablename));

%Selecting cell indices over which to compute place fields
if islogical(decparams.cellidx)
    cellidx = find(decparams.cellidx(:)' & sum(Srep(traintidx,:), 1, 'omitnan') > decparams.nspk_th);
else
    cellidx = decparams.cellidx(sum(Srep(traintidx,decparams.cellidx), 1, 'omitnan') > decparams.nspk_th);
end

%Subsetting spike trains across cells.
spikeTrain = Srep(:,cellidx);
spkCount = spkCount(:,cellidx);

X = Nav.(decparams.Xvariablename);
Y = Nav.(decparams.Yvariablename);

%number of cells selected for decoding
ncells = size(spikeTrain, 2);

%number of X bins
nXbins = numel(decparams.Xbinedges) - 1;

%number of Y bins
nYbins = numel(decparams.Ybinedges) - 1;

%%
%Discretizing X according to decparams.Xbinedges
X_discrete = discretize(X, decparams.Xbinedges);

%Discretizing Y according to decparams.Ybinedges.
Y_discrete = discretize(Y, decparams.Ybinedges);
Y_discrete = NaN(size(Y));
Y_discrete(Y == 1) = 1;
Y_discrete(Y == -1) = 2;

%%
%Subsetting X_discrete, Y_discrete and spikeTrain on the training set
X_discrete_trainset = X_discrete(traintidx);
Y_discrete_trainset = Y_discrete(traintidx);
spkTrain_trainset = spikeTrain(traintidx,:);

%Computing occupancy map (same for all cells) on the train set
flat = decparams.scalingFactor * ones(size(X_discrete_trainset));
occmap = Compute2DMap(Y_discrete_trainset, X_discrete_trainset, flat, nYbins, nXbins);

%Removing occupancy for position bins below the occupancy threshold
occmap(occmap <= decparams.occ_th) = NaN;

%Smoothing the occupancy map with a 2D gaussian window.
occmap = GaussianSmooth(occmap, [decparams.XsmthNbins decparams.YsmthNbins]);

%Computing and smoothing the spike count map for each cell
scmap = NaN(ncells, nXbins, nYbins);
for icell = 1:ncells
    scmap(icell,:,:) = Compute2DMap(Y_discrete_trainset, X_discrete_trainset, spkTrain_trainset(:,icell), nYbins, nXbins);
    scmap(icell,isnan(occmap)) = NaN;
    scmap(icell,:,:) = GaussianSmooth(squeeze(scmap(icell,:,:)), [decparams.XsmthNbins decparams.YsmthNbins]);
end

%Calculating the place field x direction maps by dividing scmap and occmap
occmap = permute(occmap, [3 1 2]);%permuting dimension for convenience
mapXY = scmap ./ occmap;

%%
%number of data points
ntimepts = size(spkCount, 1);

%Initializing decoded variables
XDecMax = NaN(ntimepts,1);
YDecMax = NaN(ntimepts,1);
XDecMean = NaN(ntimepts,1);
YDecMean = NaN(ntimepts,1);
%Computing decoded positions for data points that are not included in the
%train set.
[YDecMax(~traintidx), XDecMax(~traintidx), YDecMean(~traintidx), XDecMean(~traintidx)] = ...
    ComputeBayesMAP2D(mapXY, spkCount(~traintidx,:), decparams.dectimewin);

%%
%Doing the same thing now with cross-validated data on the train set.
%First defining a partition of the data for k-fold cross-validation. NB: we
%should normally be more careful about the fact that the spike count data
%are actually smoothed over time...
ntimepts_trainset = sum(traintidx);
cv = crossvalPartition(ntimepts_trainset, decparams.kfold);

%Computing the place field using k-fold cross-validation
mapXY_cv = NaN(ncells, nXbins, nYbins, decparams.kfold);
XDecMax_cv = NaN(ntimepts_trainset,1);
YDecMax_cv = NaN(ntimepts_trainset,1);
XDecMean_cv = NaN(ntimepts_trainset,1);
YDecMean_cv = NaN(ntimepts_trainset,1);
X_discrete_trainset = X_discrete(traintidx);
Y_discrete_trainset = Y_discrete(traintidx);
spkTrain_trainset = spikeTrain(traintidx,:);
spkCount_trainset = spkCount(traintidx,:);
for i = 1:decparams.kfold
    %Subsetting X and spiketrain according to the train set of the
    %current fold
    Xtraining = X_discrete_trainset(cv.trainsets{i});
    Ytraining = Y_discrete_trainset(cv.trainsets{i});
    Spktraining = spkTrain_trainset(cv.trainsets{i},:);
    
    %Computing occupancy map for the current fold
    flat = decparams.scalingFactor * ones(size(Xtraining));
    occmap_cv = Compute2DMap(Ytraining, Xtraining, flat, nYbins, nXbins);
    occmap_cv(occmap_cv <= decparams.occ_th) = NaN;
    occmap_cv = GaussianSmooth(occmap_cv, [decparams.XsmthNbins decparams.YsmthNbins]);
    
    %Computing the spike count map and place field of each cell for the
    %current fold
    for icell = 1:ncells
        scmap_cv = Compute2DMap(Ytraining, Xtraining, Spktraining(:,icell), nYbins, nXbins);
        scmap_cv(isnan(occmap_cv)) = NaN;
        scmap_cv = GaussianSmooth(scmap_cv, [decparams.XsmthNbins decparams.YsmthNbins]);
        mapXY_cv(icell,:,:,i) = scmap_cv ./ occmap_cv;
    end
end

%Now that we've got cross-validated place fields for the train set, we can
%compute decoded positions on the train set using the same k-fold
%partition.
for i = 1:decparams.kfold
    spkCountTest = spkCount_trainset(cv.testsets{i},:);
    
    [YDecMax_cv(cv.testsets{i}), XDecMax_cv(cv.testsets{i}),...
     YDecMean_cv(cv.testsets{i}), XDecMean_cv(cv.testsets{i})] = ComputeBayesMAP2D(mapXY_cv(:,:,:,i), spkCountTest, decparams.dectimewin);
end

%Filling in cross-validated decoded positions for the train set.
XDecMax(traintidx) = XDecMax_cv;
YDecMax(traintidx) = YDecMax_cv;
XDecMean(traintidx) = XDecMean_cv;
YDecMean(traintidx) = YDecMean_cv;

%%
%Populate the output structure with results to be saved
decparams.traintidx = traintidx;
Dec.decparams = decparams;

ncells_orig = size(Srep, 2);
Dec.Xbincenters = decparams.Xbinedges(1:end-1) + decparams.Xbinsize / 2;
Dec.nXbins = nXbins;
Dec.Ybincenters = decparams.Ybinedges(1:end-1) + decparams.Ybinsize / 2;
Dec.nYbins = nYbins;
Dec.mapXY = NaN(ncells_orig, nXbins, nYbins);
Dec.mapXY_cv = NaN(ncells_orig, nXbins, nYbins, decparams.kfold);
Dec.occmap = NaN(1, nXbins);

Dec.mapXY(cellidx,:,:) = mapXY;
Dec.mapXY_cv(cellidx,:,:,:) = mapXY_cv;
Dec.occmap = occmap;

Dec.X = X_discrete;
Dec.Y = Y_discrete;

Dec.XDecMax = XDecMax;
Dec.YDecMax = YDecMax;
Dec.XDecMean = XDecMean;
Dec.YDecMean = YDecMean;

Dec.XErrMax = (XDecMax - X_discrete);
Dec.YErrMax = YDecMax - Y_discrete;
Dec.XErrMean = (XDecMean - X_discrete);
Dec.YErrMean = YDecMean - Y_discrete;

%Calculating the distribution of decoded variables as a function of the 
%actual variables over the training set (just to be able to quickly check 
%the quality of decoding on that portion of the data).
Dec.Xdecmat = Compute2DMap(Dec.X(traintidx), Dec.XDecMax(traintidx), ones(size(Dec.X(traintidx))), Dec.nXbins, Dec.nXbins);
Dec.Ydecmat = Compute2DMap(Dec.Y(traintidx), Dec.YDecMax(traintidx), ones(size(Dec.Y(traintidx))), Dec.nYbins, Dec.nYbins);
end