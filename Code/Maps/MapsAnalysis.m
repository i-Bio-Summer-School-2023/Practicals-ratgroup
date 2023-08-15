function Maps = MapsAnalysis(Nav, Srep, mapsparams)
% MapsAnalysis - Estimates two-dimensional place fields and their significance.
%
%   Maps = MapsAnalysis(Nav, Srep, mapsparams)
%
%   This function estimates two-dimensional maps and their significance using either
%   shuffling or model comparison on cross-validated predictions.
%
%   INPUTS:
%   - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables
%   - Srep: Array of responses (ntimes x ncells).
%   - mapsparams: Struct containing parameters for place field estimation.
%
%   OUTPUT:
%   - Maps: Struct containing place field analysis results, including fields such as:
%     * map: Two-dimensional place fields (ncells x nYbins x nXbins array)
%     * map_cv: Two-dimensional place fields estimated using k-fold 
%       cross-validation (ncells x nYbins x nXbins x k-fold array)
%     * map_SE: Jackknife estimate of standard error for place fields, 
%       (ncells x nYbins x nXbins array)
%     * mapsparams: Structure of parameters used for analysis
%     * Xbincenters: Bin centers along X-axis
%     * Ybincenters: Bin centers along Y-axis
%     * occmap: Occupancy map, a nYbins x nXbins arrays (scaled by mapsparams.scalingFactor)
%     * SI: Spatial information for each cell (ncells x 1 array).
%     * SparsityIndex: Sparsity index for each cell.
%     * SelectivityIndex: Selectivity index for each cell.
%     * DirectionalityIndex: Directionality index for each cell.
%     * SI_pval: P-values for spatial information, based on shuffle controls
%     * SparsityIndex_pval: P-values for sparsity index, based on shuffle controls
%     * SelectivityIndex_pval: P-values for selectivity index, based on shuffle controls
%     * DirectionalityIndex_pval: P-values for directionality index, based on shuffle controls
%     * EV: Cross-validated percentage of explained variance from place field model
%     * EV_cst: Cross-validated percentage of explained variance from constant mean model
%     * LLH: Cross-validated Log likelihood from place field model
%     * LLH_cst: Cross-validated Log likelihood from constant mean model
%     * LLH_pval: P-values for model comparison from likelihood ratio test
%
%   USAGE:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    mapsparams = SetMapsParams(Nav,Srep);
%
%    % Change parameters of mapsparams here if needed. For instance
%    mapsparams.Yvariablename = [];%compute 1D maps along X variable.
%
%    Maps = MapsAnalysis(Nav, Spk.spikeTrain, mapsparams);
%
%   SEE ALSO:
%   ComputeMap, GaussianSmooth, computeEV, computeLLH_normal, crossvalPartition,
%   getSpatialinfo, getSparsity, getSelectivity, getDirectionality,
%   lratiotest (requires econometrics toolbox).
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological recordings"
%
%%

%If no Y variable are indicated in mapsparams, we'll just compute a 1D map
if ~isempty(mapsparams.Xvariablename)
    X = Nav.(mapsparams.Xvariablename);
else
    X = ones(size(Nav.sampleTimes));
    mapsparams.Xbinedges = 1;
    mapsparams.XsmthNbins = 0;
end

%If no Y variable are indicated in mapsparams, we'll just compute a 1D map
if ~isempty(mapsparams.Yvariablename)
    Y = Nav.(mapsparams.Yvariablename);
else
    Y = ones(size(Nav.sampleTimes));
    mapsparams.Ybinedges = 1;
    mapsparams.YsmthNbins = 0;
end

%Selecting time indices over which to compute maps, according to parameters
%defined in mapsparams.subset
tidx = true(size(X));
pnames = fieldnames(mapsparams.subset);
fnames = fieldnames(Nav);
for i = 1:numel(pnames)
    if ismember(pnames{i},fnames)
        fn = str2func(mapsparams.subset.([pnames{i} '_op']));
        tidx = tidx & fn(Nav.(pnames{i}), mapsparams.subset.(pnames{i}));
    elseif ~strcmp(pnames{i}(end-2:end),'_op')
        warning('some fields of mapsparams.subset are not matching fields of Nav')
    end
end

%Selecting time indices where X and Y are in the range of bins
tidx = tidx &...
       X >= mapsparams.Xbinedges(1) & X <= mapsparams.Xbinedges(end) &...
       Y >= mapsparams.Ybinedges(1) & Y <= mapsparams.Ybinedges(end) &...
       ~isnan(X) & ~isnan(Y);

%Selecting cell indices for which to compute maps
if islogical(mapsparams.cellidx)
    cellidx = find(mapsparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > mapsparams.nspk_th);
else
    cellidx = mapsparams.cellidx(sum(Srep(tidx,mapsparams.cellidx), 1, 'omitnan') > mapsparams.nspk_th);
end

%Subsetting spikeTrain, X and Y.
spikeTrain = Srep(tidx,cellidx);
X = X(tidx);
Y = Y(tidx);

%number of bins along X
nXbins = max(1,numel(mapsparams.Xbinedges) - 1);

%number of bins along Y
nYbins = max(1,numel(mapsparams.Ybinedges) - 1);

%If no Y variable are indicated, we'll just compute a 1D maps so
%bining parameters are changed accordingly
if isempty(mapsparams.Yvariablename)
    %number of bins along Y
    nYbins = 1;
end

%number of selected cells
ncells = size(spikeTrain, 2);

%number of data points
ntimepts = size(spikeTrain, 1);

%%
%Discretizing X position vectors according to mapsparams.Xbinedges
X_discrete = discretize(X, mapsparams.Xbinedges);

%Discretizing Y position vectors according to mapsparams.Ybinedges 
Y_discrete = discretize(Y, mapsparams.Ybinedges);

%%
%Computing occupancy map (same for all cells)
flat =  mapsparams.scalingFactor * ones(size(X_discrete));
occmap = ComputeMap(X_discrete, Y_discrete, flat, nXbins, nYbins);

%Removing occupancy for bins below the occupancy threshold
occmap(occmap <= mapsparams.occ_th) = NaN;

%Smoothing the occupancy map with a 2D gaussian window
occmap = GaussianSmooth(occmap, [mapsparams.YsmthNbins mapsparams.XsmthNbins]);

%Computing and smoothing spike count map for each cell
scmap = NaN(ncells, nYbins, nXbins);
parfor icell = 1:ncells
    scmapcell = ComputeMap(X_discrete, Y_discrete, spikeTrain(:,icell), nXbins, nYbins);
    scmapcell(isnan(occmap)) = NaN;
    scmapcell = GaussianSmooth(scmapcell, [mapsparams.YsmthNbins mapsparams.XsmthNbins]);
    scmap(icell,:,:) = scmapcell;
end

%Calculating the maps by dividing scmap and occmap
mapXY = scmap ./ permute(occmap, [3 1 2]);
occmap = squeeze(occmap);

%%
%Quantifying selectivity by computing the spatial information (SI),
%the sparsity index, the selectivity index and the directionality index.
SI = NaN(ncells, 1);
SparsityIndex = NaN(ncells, 1);
SelectivityIndex = NaN(ncells, 1);
DirectionalityIndex = NaN(ncells, 1);
for icell = 1:ncells
    SI(icell) = getSpatialinfo(mapXY(icell,:), occmap(:));
    SparsityIndex(icell) = getSparsity(mapXY(icell,:), occmap(:));
    SelectivityIndex(icell) = getSelectivity(mapXY(icell,:));
    if nYbins == 2
        DirectionalityIndex(icell) = getDirectionality(mapXY(icell,1,:), mapXY(icell,2,:));
    end
end

%%
%Computing shuffle controls by randomly shifting time bins of positions and
%calculate the selectivity metrics for each shuffle control
SI_Shf = NaN(ncells, mapsparams.nShuffle);
SparsityIndex_Shf = NaN(ncells, mapsparams.nShuffle);
SelectivityIndex_Shf = NaN(ncells, mapsparams.nShuffle);
DirectionalityIndex_Shf = NaN(ncells, mapsparams.nShuffle);
%Initializing the random number generator for reproducibility purposes
nShf = mapsparams.nShuffle;
%Calculating the place field for each shuffle permutation
parfor icell = 1:ncells
    s = RandStream('mt19937ar','Seed',icell);
    for iperm  = 1:nShf
        %Shifting X and Y by a random amount larger than 1 second
        tshift = randi(s,ntimepts - 2 * mapsparams.sampleRate) + 1 * mapsparams.sampleRate;
        X_discrete_shf = circshift(X_discrete, round(tshift));
        Y_discrete_shf = circshift(Y_discrete, round(tshift));
        
        %Computing maps after shuffling
        scmap_shf = ComputeMap(X_discrete_shf, Y_discrete_shf, spikeTrain(:,icell), nXbins, nYbins);
        scmap_shf(isnan(occmap)) = NaN;
        scmap_shf = GaussianSmooth(scmap_shf, [mapsparams.YsmthNbins mapsparams.XsmthNbins]);
        mapX_shf = scmap_shf ./ occmap;

        %saving only the spatial selectivity metrics for each permutation
        SI_Shf(icell,iperm) = getSpatialinfo(mapX_shf(:), occmap(:));
        SparsityIndex_Shf(icell,iperm) = getSparsity(mapX_shf(:), occmap(:));
        SelectivityIndex_Shf(icell,iperm) = getSelectivity(mapX_shf(:));
        if nYbins == 2
            DirectionalityIndex_Shf(icell,iperm) = getDirectionality(mapX_shf(:,1), mapX_shf(:,2));
        end
    end
end

%Computing p-values from the distribution of selectivity measures obtained
%from the shuffle controls
SI_pval = sum(SI_Shf > SI, 2) / mapsparams.nShuffle;
SparsityIndex_pval = sum(SparsityIndex_Shf > SparsityIndex, 2) / mapsparams.nShuffle;
SelectivityIndex_pval = sum(SelectivityIndex_Shf > SelectivityIndex, 2) / mapsparams.nShuffle;
DirectionalityIndex_pval = sum(DirectionalityIndex_Shf > DirectionalityIndex, 2) / mapsparams.nShuffle;

%%
%Defining a partition of the data for k-fold cross-validation
cv = crossvalPartition(ntimepts, mapsparams.kfold);

%Computing the spike train predicted from the place field using k-fold 
%cross-validation
mapXY_cv = NaN(ncells, nYbins, nXbins, mapsparams.kfold);
Ypred = NaN(ntimepts,ncells);
for i = 1:mapsparams.kfold
    %Subsetting X and spiketrain according to the train set of the
    %current fold
    Xtraining = X_discrete(cv.trainsets{i});
    Ytraining = Y_discrete(cv.trainsets{i});
    Spktraining = spikeTrain(cv.trainsets{i},:);
    
    %Computing occupancy map for the current fold
    flat = mapsparams.scalingFactor * ones(size(Xtraining));
    occmap_cv = ComputeMap(Xtraining, Ytraining, flat, nXbins, nYbins);
    occmap_cv(occmap_cv <= mapsparams.occ_th) = NaN;
    occmap_cv = GaussianSmooth(occmap_cv, [mapsparams.YsmthNbins mapsparams.XsmthNbins]);
    
    %Computing the spike count map and place field of each cell for the
    %current fold
    parfor icell = 1:ncells
        scmap_cv = ComputeMap(Xtraining, Ytraining, Spktraining(:,icell), nXbins, nYbins);
        scmap_cv(isnan(occmap_cv)) = NaN;
        scmap_cv = GaussianSmooth(scmap_cv, [mapsparams.YsmthNbins mapsparams.XsmthNbins]);
        mapXY_cv(icell,:,:,i) = scmap_cv ./ occmap_cv;
    end
    
    %Subsetting X and Y according to the test set of the current fold
    Xtest = X_discrete(cv.testsets{i});
    Ytest = Y_discrete(cv.testsets{i});
    
    %Computing the spike train predicted on the test set from the place 
    %computed from the train set
    for icell = 1:ncells
        XYlinidx = sub2ind([ncells, nYbins, nXbins, mapsparams.kfold], icell*ones(size(Xtest)), Ytest, Xtest, i*ones(size(Xtest)));
        Ypred(cv.testsets{i},icell) = mapXY_cv(XYlinidx) * mapsparams.scalingFactor;
    end
end


%Now computing the spike train predicted from the mean firing rate of the 
%cell using the same k-fold partition as above
Ypred_cst = NaN(ntimepts,ncells);
for i = 1:mapsparams.kfold
    for icell = 1:ncells
        Ypred_cst(cv.testsets{i},icell) = mean(spikeTrain(cv.trainsets{i},icell), 'omitnan');
    end
end

%Computing the percentage of explained variance and the log likelihood from
%the spike trains predicted by the place field model
EV = NaN(ncells,1);
LLH = NaN(ncells,1);
for icell = 1:ncells
    %Percentage of explained variance
    EV(icell) =  computeEV(spikeTrain(:, icell), Ypred(:, icell));
    %Log likelihood
    LLH(icell) = computeLLH_normal(spikeTrain(:, icell), Ypred(:, icell));
end

%Same for the spike train predicted by the mean constant model
EV_cst = NaN(ncells,1);
LLH_cst = NaN(ncells,1);
for icell = 1:ncells
    %Percentage of explained variance
    EV_cst(icell) =  computeEV(spikeTrain(:, icell), Ypred_cst(:, icell));
    %Log likelihood
    LLH_cst(icell) = computeLLH_normal(spikeTrain(:, icell), Ypred_cst(:, icell));
end

%Comparing the place field model to the constant mean model by performing a
%likelihood ratio test
LLH_pval = NaN(ncells,1);
goodidx = LLH > LLH_cst;
LLH_pval(~goodidx) = 1;
dof = sum(occmap(:) > 0) - 1;
if sum(goodidx) > 0
    [~, LLH_pval(goodidx)] = lratiotest(LLH(goodidx), LLH_cst(goodidx), dof);
end


%Computing a Jacknife estimate of the standard error
mapXY_SE = sqrt((mapsparams.kfold - 1)./mapsparams.kfold * sum((mapXY_cv - mapXY).^2, 4));

%%
%Populate the output structure with results to be saved
mapsparams.tidx = tidx;
Maps.mapsparams = mapsparams;

if nXbins > 1
    Maps.Xbincenters = mapsparams.Xbinedges(1:end-1) + diff(mapsparams.Xbinedges) / 2;
else
    Maps.Xbincenters = 1;
end
if nYbins > 1
    Maps.Ybincenters = mapsparams.Ybinedges(1:end-1) + diff(mapsparams.Ybinedges) / 2;
else
    Maps.Ybincenters = 1;
end

%Interpolating in case there are infinite values among the bin edges
%(probably at the beginiing or end in order to clamp de signal)
if nXbins > 1 && sum(~isinf(Maps.Xbincenters)) > 0
    xb = 1:nXbins;
    Maps.Xbincenters = interp1(xb(~isinf(Maps.Xbincenters)), Maps.Xbincenters(~isinf(Maps.Xbincenters)), xb, 'linear', 'extrap');
end
if nYbins > 1 && sum(~isinf(Maps.Ybincenters)) > 0
    yb = 1:nYbins;
    Maps.Ybincenters = interp1(yb(~isinf(Maps.Ybincenters)), Maps.Ybincenters(~isinf(Maps.Ybincenters)), yb, 'linear', 'extrap');
end

ncells_orig = size(Srep, 2);
Maps.map = NaN(ncells_orig, nYbins, nXbins);
Maps.map_cv = NaN(ncells_orig, nYbins, nXbins, mapsparams.kfold);
Maps.map_SE = NaN(ncells_orig, nYbins, nXbins);
Maps.occmap = NaN(1, nYbins, nXbins);
Maps.SI = NaN(ncells_orig, 1);
Maps.SparsityIndex = NaN(ncells_orig, 1);
Maps.SelectivityIndex = NaN(ncells_orig, 1);
Maps.DirectionalityIndex = NaN(ncells_orig, 1);
Maps.SI_pval = NaN(ncells_orig, 1);
Maps.SparsityIndex_pval = NaN(ncells_orig, 1);
Maps.SelectivityIndex_pval = NaN(ncells_orig, 1);
Maps.DirectionalityIndex_pval = NaN(ncells_orig, 1);
Maps.EV = NaN(ncells_orig,1);
Maps.EV_cst = NaN(ncells_orig,1);
Maps.LLH = NaN(ncells_orig,1);
Maps.LLH_cst = NaN(ncells_orig,1);
Maps.LLH_pval = NaN(ncells_orig,1);

Maps.map(cellidx,:,:,:) = mapXY;
Maps.map_cv(cellidx,:,:,:) = mapXY_cv;
Maps.map_SE(cellidx,:,:) = mapXY_SE;
Maps.occmap = occmap;
Maps.SI(cellidx,:) = SI;
Maps.SparsityIndex(cellidx,:) = SparsityIndex;
Maps.SelectivityIndex(cellidx,:) = SelectivityIndex;
Maps.DirectionalityIndex(cellidx) = DirectionalityIndex;
Maps.SI_pval(cellidx,:) = SI_pval;
Maps.SparsityIndex_pval(cellidx,:) = SparsityIndex_pval;
Maps.SelectivityIndex_pval(cellidx,:) = SelectivityIndex_pval;
Maps.DirectionalityIndex_pval(cellidx) = DirectionalityIndex_pval;
Maps.EV(cellidx) = EV;
Maps.EV_cst(cellidx) = EV_cst;
Maps.LLH(cellidx) = LLH;
Maps.LLH_cst(cellidx) = LLH_cst;
Maps.LLH_pval(cellidx) = LLH_pval;
end