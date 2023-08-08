function Maps = MapsAnalyses2D(Nav, Srep, mapsparams)
MapsAnalyses2D - Estimates two-dimensional place fields and their significance.

  Maps = MapsAnalyses2D(Nav, Srep, mapsparams)

  This function estimates two-dimensional maps and their significance
  using either shuffling or model comparison on cross-validated predictions.

  Inputs:
  - Nav: Struct containing navigation-related data variables on which
  responses should be mapped onto
  - Srep: array of responses (ntimes x ncells)
  - mapsparams: Struct containing parameters for place field estimation.

  Outputs:
  - Maps: Struct containing place field analysis results, including fields such as:
    * mapXY: Two-dimensional place fields
    * mapXY_cv: Two-dimensional place fields estimated using k-fold cross-validation
    * mapXY_SE: Jackknife estimate of standard error for place fields
    * mapsparams: Parameters used for analysis
    * Xbincenters: Bin centers along X-axis
    * Ybincenters: Bin centers along Y-axis
    * occmap: Occupancy map, a 1 x nXbins arrays (scaled by
      mapsparams.scalingFactor)
    * SI: Spatial information
    * SparsityIndex: Sparsity index for each cell.
    * SelectivityIndex: Selectivity index for each cell.
    * DirectionalityIndex: Directionality index for each cell.
    * SI_pval: P-values for spatial information, based on
      shuffle controls
    * SparsityIndex_pval: P-values for sparsity index, based on
      shuffle controls
    * SelectivityIndex_pval: P-values for selectivity index, based on
      shuffle controls
    * DirectionalityIndex_pval: P-values for directionality index, based on
      shuffle controls
    * EV: cross-validated percentage of explained variance from place field model
    * EV_cst: cross-validated percentage of explained variance from constant mean model
    * LLH: cross-validated Log likelihood from place field model
    * LLH_cst: cross-validated Log likelihood from constant mean model
    * LLH_pval: P-values  for model comparison from likelihood ratio test

  Usage:
   Nav = LoaddataNav(loadparams);
   Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
   Srep = Spk.spikeTrain;
   mapsparams = DefineMapsParams(Nav,Spk);
   Maps2D = MapsAnalyses2D(Nav, Spk.spikeTrain, mapsparams);

See also: Compute2DMap, GaussianSmooth, computeEV, computeLLH_normal,
          crossvalPartition, SpatialInfo, FieldSparsity, FieldSelectivity,
          FieldDirectionality, lratiotest

written by J.Fournier 08/2023 for the iBio Summer school
%

%%
if isempty(mapsparams.Xvariablename)
    mapsparams.Xvariablename = 'Xpos';
end
if isempty(mapsparams.Yvariablename)
    mapsparams.Yvariablename = 'XDir';
end

%Selecting time and cell indices over which to compute place fields
tidx = ismember(Nav.Condition, mapsparams.condition) &...
       ismember(Nav.XDir, mapsparams.dir) &...
       ismember(Nav.laptype, mapsparams.laptype) &...
       Nav.Spd >= mapsparams.spdthreshold &...
       ~isnan(Nav.(mapsparams.Xvariablename)) & ...
       ~isnan(Nav.(mapsparams.Yvariablename));

if islogical(mapsparams.cellidx)
    cellidx = find(mapsparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > mapsparams.nspk_th);
else
    cellidx = mapsparams.cellidx(sum(Srep(tidx,mapsparams.cellidx), 1, 'omitnan') > mapsparams.nspk_th);
end

%Subsetting spikeTrain, X and XDir.
spikeTrain = Srep(tidx,cellidx);
X = Nav.(mapsparams.Xvariablename)(tidx);
Y = Nav.(mapsparams.Yvariablename)(tidx);

%number of cells selected for place field analysis
ncells = size(spikeTrain, 2);

%number of position bins
nXbins = numel(mapsparams.Xbinedges) - 1;

%number of direction bins
nYbins = numel(mapsparams.Ybinedges) - 1;

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
occmap = Compute2DMap(Y_discrete, X_discrete, flat, nYbins, nXbins);

%Removing occupancy for positions below the occupancy threshold
occmap(occmap <= mapsparams.occ_th) = NaN;

%Smoothing the occupancy map with a 2D gaussian window
occmap = GaussianSmooth(occmap, [mapsparams.XsmthNbins mapsparams.YsmthNbins]);

%Computing and smoothing spike count map for each cell
scmap = NaN(ncells, nXbins, nYbins);
for icell = 1:ncells
    scmap(icell,:,:) = Compute2DMap(Y_discrete, X_discrete, spikeTrain(:,icell), nYbins, nXbins);
    scmap(icell,isnan(occmap)) = NaN;
    scmap(icell,:,:) = GaussianSmooth(squeeze(scmap(icell,:,:)), [mapsparams.XsmthNbins mapsparams.YsmthNbins]);
end

%Calculating the place field x direction maps by dividing scmap and occmap
occmap = permute(occmap, [3 1 2]);%permuting dimension for element-wise division
mapXY = scmap ./ occmap;
occmap = squeeze(occmap);

%%
%Quantifying position selectivity by computing the spatial information (SI),
%the sparsity index, the selectivity index and the directionality index.
SI = NaN(ncells, nYbins);
SparsityIndex = NaN(ncells, nYbins);
SelectivityIndex = NaN(ncells, nYbins);
DirectionalityIndex = NaN(ncells, 1);
for icell = 1:ncells
    for iy = 1:nYbins
        SI(icell,iy) = SpatialInfo(mapXY(icell,:,iy), occmap(:,iy));
        SparsityIndex(icell,iy) = FieldSparsity(mapXY(icell,:,iy), occmap(:,iy));
        SelectivityIndex(icell,iy) = FieldSelectivity(mapXY(icell,:,iy));
    end
    DirectionalityIndex(icell) = FieldDirectionality(mapXY(icell,:,1), mapXY(icell,:,2));
end

%%
%Computing shuffle controls by randomly shifting time bins of positions and
%calculate the selectivity metrics for each shuffle control
SI_Shf = NaN(ncells, nYbins, mapsparams.nShuffle);
SparsityIndex_Shf = NaN(ncells, nYbins, mapsparams.nShuffle);
SelectivityIndex_Shf = NaN(ncells, nYbins, mapsparams.nShuffle);
DirectionalityIndex_Shf = NaN(ncells, mapsparams.nShuffle);
%Initializing the random number generator for reproducibility purposes
s = RandStream('mt19937ar','Seed',0);
%Calculating the place field for each shuffle permutation
for iperm  = 1:mapsparams.nShuffle
    tshift = randi(s,ntimepts - 2 * mapsparams.sampleRate) + 1 * mapsparams.sampleRate;%(Check here)
    X_discrete_shf = circshift(X_discrete, round(tshift));%(Check here)
    Y_discrete_shf = circshift(Y_discrete, round(tshift));%(Check here)
    for icell = 1:ncells
        scmap_shf = Compute2DMap(Y_discrete_shf, X_discrete_shf, spikeTrain(:,icell), nYbins, nXbins);
        scmap_shf(isnan(occmap)) = NaN;
        scmap_shf = GaussianSmooth(scmap_shf, [mapsparams.XsmthNbins mapsparams.YsmthNbins]);
        mapX_shf = scmap_shf ./ occmap;
        
        %saving only the spatial selectivity metrics for each permutation
        for iy = 1:nYbins
            SI_Shf(icell,iy,iperm) = SpatialInfo(mapX_shf(:,iy), occmap(:,iy));
            SparsityIndex_Shf(icell,iy,iperm) = FieldSparsity(mapX_shf(:,iy), occmap(:,iy));
            SelectivityIndex_Shf(icell,iy,iperm) = FieldSelectivity(mapX_shf(:,iy));
        end
        DirectionalityIndex_Shf(icell,iperm) = FieldDirectionality(mapX_shf(:,1), mapX_shf(:,2));
    end
end

%Computing p-values from the distribution of selectivity measures obtained
%from the shuffle controls
SI_pval = sum(SI_Shf > SI, 3) / mapsparams.nShuffle;
SparsityIndex_pval = sum(SparsityIndex_Shf > SparsityIndex, 3) / mapsparams.nShuffle;
SelectivityIndex_pval = sum(SelectivityIndex_Shf > SelectivityIndex, 3) / mapsparams.nShuffle;
DirectionalityIndex_pval = sum(DirectionalityIndex_Shf > DirectionalityIndex, 2) / mapsparams.nShuffle;

%%
%Defining a partition of the data for k-fold cross-validation
cv = crossvalPartition(ntimepts, mapsparams.kfold);

%Computing the spike train predicted from the place field using k-fold 
%cross-validation
mapXY_cv = NaN(ncells, nXbins, nYbins, mapsparams.kfold);
Ypred = NaN(ntimepts,ncells);
for i = 1:mapsparams.kfold
    %Subsetting X and spiketrain according to the train set of the
    %current fold
    Xtraining = X_discrete(cv.trainsets{i});
    Dirtraining = Y_discrete(cv.trainsets{i});
    Spktraining = spikeTrain(cv.trainsets{i},:);
    
    %Computing occupancy map for the current fold
    flat = mapsparams.scalingFactor * ones(size(Xtraining));
    occmap_cv = Compute2DMap(Dirtraining, Xtraining, flat, nYbins, nXbins);
    occmap_cv(occmap_cv <= mapsparams.occ_th) = NaN;
    occmap_cv = GaussianSmooth(occmap_cv, [mapsparams.XsmthNbins mapsparams.YsmthNbins]);
    
    %Computing the spike count map and place field of each cell for the
    %current fold
    for icell = 1:ncells
        scmap_cv =Compute2DMap(Dirtraining, Xtraining, Spktraining(:,icell), nYbins, nXbins);
        scmap_cv(isnan(occmap_cv)) = NaN;
        scmap_cv = GaussianSmooth(scmap_cv, [mapsparams.XsmthNbins mapsparams.YsmthNbins]);
        mapXY_cv(icell,:,:,i) = scmap_cv ./ occmap_cv;
    end
    
    %Subsetting X and Y according to the test set of the current fold
    Xtest = X_discrete(cv.testsets{i});
    Ytest = Y_discrete(cv.testsets{i});
    
    %Computing the spike train predicted on the test set from the place 
    %computed from the train set
    for icell = 1:ncells
        XYlinidx = sub2ind([ncells, nXbins, nYbins, mapsparams.kfold], icell*ones(size(Xtest)), Xtest, Ytest, i*ones(size(Xtest)));
        try
        Ypred(cv.testsets{i},icell) = mapXY_cv(XYlinidx) * mapsparams.scalingFactor;%(Check here)
        catch
            keyboard
        end
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
[~, LLH_pval(goodidx)] = lratiotest(LLH(goodidx), LLH_cst(goodidx), dof);


%Computing a Jacknife estimate of the standard error
mapXY_SE = sqrt((mapsparams.kfold - 1)./mapsparams.kfold * sum((mapXY_cv - mapXY).^2, 4));

%%
%Populate the output structure with results to be saved
mapsparams.tidx = tidx;
Maps.mapsparams = mapsparams;

Maps.Xbincenters = mapsparams.Xbinedges(1:end-1) + mapsparams.Xbinsize / 2;
Maps.Ybincenters = mapsparams.Xbinedges(1:end-1) + mapsparams.Ybinsize / 2;

ncells_orig = size(Srep, 2);
nXbins = numel(Maps.Xbincenters);
Maps.mapXY = NaN(ncells_orig, nXbins, nYbins);
Maps.mapXY_cv = NaN(ncells_orig, nXbins, nYbins, mapsparams.kfold);
Maps.mapXY_SE = NaN(ncells_orig, nXbins, nYbins);
Maps.occmap = NaN(1, nXbins, nYbins);
Maps.SI = NaN(ncells_orig, nYbins);
Maps.SparsityIndex = NaN(ncells_orig, nYbins);
Maps.SelectivityIndex = NaN(ncells_orig, nYbins);
Maps.DirectionalityIndex = NaN(ncells_orig, 1);
Maps.SI_pval = NaN(ncells_orig, nYbins);
Maps.SparsityIndex_pval = NaN(ncells_orig, nYbins);
Maps.SelectivityIndex_pval = NaN(ncells_orig, nYbins);
Maps.DirectionalityIndex_pval = NaN(ncells_orig, 1);
Maps.EV = NaN(ncells_orig,1);
Maps.EV_cst = NaN(ncells_orig,1);
Maps.LLH = NaN(ncells_orig,1);
Maps.LLH_cst = NaN(ncells_orig,1);
Maps.LLH_pval = NaN(ncells_orig,1);

Maps.mapXY(cellidx,:,:,:) = mapXY;
Maps.mapXY_cv(cellidx,:,:,:) = mapXY_cv;
Maps.mapXY_SE(cellidx,:,:) = mapXY_SE;
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