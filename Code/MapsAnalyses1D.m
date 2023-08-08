function Maps = MapsAnalyses1D(Nav, Srep, mapsparams)
% Maps = MapsAnalyses1D(Nav, Srep, mapsparams) computes place fields and their significance
% either through shuffling or model comparison on cross-validated predictions in a
% one-dimensional spatial context.
%
% Input:
%   - Nav: Struct containing navigation-related data variables on which
%   responses should be mapped onto
%   Srep - ntime x ncells array of responses to map.
%   mapsparams - Structure containing parameters for place field estimation.
%
% Output:
%   Maps - Structure containing various results of the place field analysis, including:
%     - mapX: Place fields, a nCells x nXbin array.
%     - mapX_cv: Place fields estimated using k-fold cross-validation.
%     - mapX_SE: Jacknife estimate of the standard error for place fields.
%     - mapsparams: Parameters used for place field estimation.
%     - Xbincenters: Bin centers along the X-axis.
%     - occmap: Occupancy map, a 1 x nXbins arrays (scaled by
%       mapsparams.scalingFactor)
%     - SI: Spatial information for each cell.
%     - SparsityIndex: Sparsity index for each cell.
%     - SelectivityIndex: Selectivity index for each cell.
%     - SI_pval: P-values for spatial information significance, based on
%       shuffle controls
%     - SparsityIndex_pval: P-values for sparsity index significance, based on
%       shuffle controls
%     - SelectivityIndex_pval: P-values for selectivity index significance, based on
%       shuffle controls
%     - EV: cross-validated percentage of explained variance using place field model.
%     - EV_cst: cross-validated percentage of explained variance using constant mean model.
%     - LLH: cross-validated Log likelihood using place field model.
%     - LLH_cst: cross-validated Log likelihood using constant mean model.
%     - LLH_pval: P-values for model comparison using likelihood ratio test
%
% Usage example:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    mapsparams = DefineMapsParams(Nav,Spk);
%    Maps1D = MapsAnalyses1D(Nav, Spk.spikeTrain, mapsparams);
%
% See also: Compute1DMap, GaussianSmooth1D, crossvalPartition, computeEV,
%   computeLLH_normal, lratiotest
%
% written by J.Fournier 08/2023 for the iBio Summer school
%
% Note: This function provides place field analysis for a one-dimensional environment.
% For 2-dimensional analyses, consider using appropriate MapsAnalyses2D.

%%
if isempty(mapsparams.Xvariablename)
    mapsparams.Xvariablename = 'Xpos';
end

%Selecting time and cell indices over which to compute place fields
tidx = ismember(Nav.Condition, mapsparams.condition) &...
       ismember(Nav.XDir, mapsparams.dir) &...
       ismember(Nav.laptype, mapsparams.laptype) &...
       Nav.Spd >= mapsparams.spdthreshold &...
       ~isnan(Nav.(mapsparams.Xvariablename));

if islogical(mapsparams.cellidx)
    cellidx = find(mapsparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > mapsparams.nspk_th);
else
    cellidx = mapsparams.cellidx(sum(Srep(tidx,mapsparams.cellidx), 1, 'omitnan') > mapsparams.nspk_th);
end

%Subsetting spikeTrain and Xpos.
spikeTrain = Srep(tidx,cellidx);

Xpos = Nav.(mapsparams.Xvariablename)(tidx);

%number of cells selected for place field analysis
ncells = size(spikeTrain, 2);

%number of position bins
nbins = numel(mapsparams.Xbinedges) - 1;

%number of data points
ntimepts = size(spikeTrain, 1);

%%
%Discretizing position vectors according to mapsparams.binedges
Xpos_discrete = discretize(Xpos, mapsparams.Xbinedges);

%%
%Computing occupancy map (same for all cells)
flat =  mapsparams.scalingFactor * ones(size(Xpos_discrete));
occmap = Compute1DMap(Xpos_discrete, flat, nbins);

%Removing occupancy for positions below the occupancy threshold
occmap(occmap <= mapsparams.occ_th) = NaN;

%Smoothing the occupancy map with a gaussian window (mapsparams.XsmthNbins
%of sd)
occmap = GaussianSmooth1D(occmap, mapsparams.XsmthNbins);


%Computing and smoothing spike count map for each cell
scmap = NaN(ncells, nbins);
for icell = 1:ncells
    scmap(icell,:) = Compute1DMap(Xpos_discrete, spikeTrain(:,icell), nbins);
    scmap(icell,isnan(occmap)) = NaN;
    scmap(icell,:) = GaussianSmooth1D(scmap(icell,:), mapsparams.XsmthNbins);
end


%Calculating the place field maps by dividing scmap and occmap
mapX = scmap ./ occmap;


%%
%Quantifying position selectivity by computing the spatial information (SI),
%the sparsity index and the selectivity index.
SI = NaN(ncells, 1);
SparsityIndex = NaN(ncells, 1);
SelectivityIndex = NaN(ncells, 1);
for icell = 1:ncells
    SI(icell) = SpatialInfo(mapX(icell,:), occmap);
    SparsityIndex(icell) = FieldSparsity(mapX(icell,:), occmap);
    SelectivityIndex(icell) = FieldSelectivity(mapX(icell,:));
end

%%
%Computing shuffle controls by randomly shifting time bins of positions and
%calculate the selectivity metrics for each shuffle control
SI_Shf = NaN(ncells, mapsparams.nShuffle);
SparsityIndex_Shf = NaN(ncells, mapsparams.nShuffle);
SelectivityIndex_Shf = NaN(ncells, mapsparams.nShuffle);
%Initializing the random number generator for reproducibility purposes
s = RandStream('mt19937ar','Seed',0);
%Calculating the place field for each shuffle permutation
for iperm  = 1:mapsparams.nShuffle
    tshift = randi(s,ntimepts - 2 * mapsparams.sampleRate) + 1 * mapsparams.sampleRate;%(Check here)
    Xpos_discrete_shf = circshift(Xpos_discrete, round(tshift));%(Check here)
    for icell = 1:ncells
        scmap_shf = Compute1DMap(Xpos_discrete_shf, spikeTrain(:,icell), nbins);%(Check here)
        scmap_shf(isnan(occmap)) = NaN;
        scmap_shf = GaussianSmooth1D(scmap_shf, mapsparams.XsmthNbins);
        mapX_shf = scmap_shf ./ occmap;
        
        %saving only the spatial selectivity metrics for each permutation
        SI_Shf(icell,iperm) = SpatialInfo(mapX_shf, occmap);
        SparsityIndex_Shf(icell,iperm) = FieldSparsity(mapX_shf, occmap);
        SelectivityIndex_Shf(icell,iperm) = FieldSelectivity(mapX_shf);
    end
end

%Computing p-values from the distribution of selectivity measures obtained
%from the shuffle controls
SI_pval = sum(SI_Shf > SI, 2) / mapsparams.nShuffle;
SparsityIndex_pval = sum(SparsityIndex_Shf > SparsityIndex, 2) / mapsparams.nShuffle;
SelectivityIndex_pval = sum(SelectivityIndex_Shf > SelectivityIndex, 2) / mapsparams.nShuffle;

%%
%Defining a partition of the data for k-fold cross-validation
cv = crossvalPartition(ntimepts, mapsparams.kfold);

%Computing the spike train predicted from the place field using k-fold 
%cross-validation
mapX_cv = NaN(ncells, nbins, mapsparams.kfold);
Ypred = NaN(ntimepts,ncells);
for i = 1:mapsparams.kfold
    %Subsetting Xpos and spiketrain according to the train set of the
    %current fold
    Xtraining = Xpos_discrete(cv.trainsets{i});
    Spktraining = spikeTrain(cv.trainsets{i},:);
    
    %Computing occupancy map for the current fold
    flat = mapsparams.scalingFactor * ones(size(Xtraining));
    occmap_cv = Compute1DMap(Xtraining, flat, nbins);
    occmap_cv(occmap_cv <= mapsparams.occ_th) = NaN;
    occmap_cv = GaussianSmooth1D(occmap_cv, mapsparams.XsmthNbins);
    
    %Computing the spike count map and place field of each cell for the
    %current fold
    for icell = 1:ncells
        scmap_cv = Compute1DMap(Xtraining, Spktraining(:,icell), nbins);
        scmap_cv(isnan(occmap_cv)) = NaN;
        scmap_cv = GaussianSmooth1D(scmap_cv, mapsparams.XsmthNbins);
        mapX_cv(icell,:,i) = scmap_cv ./ occmap_cv;
    end
    
    %Subsetting Xpos according to the test set of the current fold
    Xtest = Xpos_discrete(cv.testsets{i});
    
    %Computing the spike train predicted on the test set from the place 
    %computed from the train set
    for icell = 1:ncells
        Ypred(cv.testsets{i},icell) = mapX_cv(icell,Xtest,i) * mapsparams.scalingFactor;%(Check here)
    end
end



%Now computing the spike train predicted from the mean firing rate of the 
%cell using the same k-fold partition as above
Ypred_cst = NaN(ntimepts,ncells);
for i = 1:mapsparams.kfold
    for icell = 1:ncells
        Ypred_cst(cv.testsets{i},icell) = nanmean(spikeTrain(cv.trainsets{i},icell));
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
dof = sum(occmap > 0) - 1;
[~, LLH_pval(goodidx)] = lratiotest(LLH(goodidx), LLH_cst(goodidx), dof);


%Computing a Jacknife estimate of the standard error (Check here)
mapX_SE = sqrt((mapsparams.kfold - 1)./mapsparams.kfold * sum((mapX_cv - mapX).^2, 3));

%%
%Populate the output structure with results to be saved
mapsparams.tidx = tidx;
Maps.mapsparams = mapsparams;

ncells_orig = size(Srep, 2);
Maps.Xbincenters = mapsparams.Xbinedges(1:end-1) + mapsparams.Xbinsize / 2;
nbins = numel(Maps.Xbincenters);
Maps.mapX = NaN(ncells_orig, nbins);
Maps.mapX_cv = NaN(ncells_orig, nbins, mapsparams.kfold);
Maps.occmap = NaN(1, nbins);
Maps.SI = NaN(ncells_orig, 1);
Maps.SparsityIndex = NaN(ncells_orig, 1);
Maps.SelectivityIndex = NaN(ncells_orig, 1);
Maps.SI_pval = NaN(ncells_orig, 1);
Maps.SparsityIndex_pval = NaN(ncells_orig, 1);
Maps.SelectivityIndex_pval = NaN(ncells_orig, 1);
Maps.EV = NaN(ncells_orig,1);
Maps.EV_cst = NaN(ncells_orig,1);
Maps.LLH = NaN(ncells_orig,1);
Maps.LLH_cst = NaN(ncells_orig,1);
Maps.LLH_pval = NaN(ncells_orig,1);

Maps.mapX(cellidx,:) = mapX;
Maps.mapX_cv(cellidx,:,:) = mapX_cv;
Maps.mapX_SE(cellidx,:) = mapX_SE;
Maps.occmap = occmap;
Maps.SI(cellidx) = SI;
Maps.SparsityIndex(cellidx) = SparsityIndex;
Maps.SelectivityIndex(cellidx) = SelectivityIndex;
Maps.SI_pval(cellidx) = SI_pval;
Maps.SparsityIndex_pval(cellidx) = SparsityIndex_pval;
Maps.SelectivityIndex_pval(cellidx) = SelectivityIndex_pval;
Maps.EV(cellidx) = EV;
Maps.EV_cst(cellidx) = EV_cst;
Maps.LLH(cellidx) = LLH;
Maps.LLH_cst(cellidx) = LLH_cst;
Maps.LLH_pval(cellidx) = LLH_pval;

%Then do the same for other conditions by changing subset indices in
%mapsparams and compare LtoR and RtoL laps, before and after training, etc

end