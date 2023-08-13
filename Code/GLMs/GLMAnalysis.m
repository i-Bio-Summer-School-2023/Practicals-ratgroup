function GLMs = GLMAnalysis(Nav, Srep, glmsparams)
% Estimates tuning curves to position and running speed using a Poisson GLM
% model. Model comparison is performed by comparing the likelihood of fitted
% data. This function computes GLMs for position-only, speed-only, and
% position x speed models, performs k-fold cross-validation, and evaluates
% model significance.
%
% Inputs:
%   Nav: Navigation data structure containing information about conditions,
%        positions and running speeds.
%   Srep: Responses, array where rows correspond to time bins and
%         columns correspond to cells. Each element represents the activity
%         of a cell in a time bin.
%   glmsparams: Structure containing GLM analysis parameters
%
% Outputs:
%   GLMs: Output structure containing the results and analysis of the GLM.
%       - glmsparams: Input parameters for the GLM analysis.
%       - Xbincenters: Bin centers for position tuning curves.
%       - Spdbincenters: Bin centers for speed tuning curves.
%       - mapX: Position tuning curves for selected cells.
%       - mapXcv: Cross-validated position tuning curves.
%       - mapX2: Speed tuning curves for selected cells.
%       - mapScv: Cross-validated speed tuning curves.
%       - pval_X1: P-values indicating significance of adding position to the model.
%       - pval_X2: P-values indicating significance of adding speed to the model.
%       - bestmodel: Best model indicator for each cell (0 for constant mean,
%         1 for position, 2 for speed, 3 for position x speed).
%       - LLH: Log likelihood values for different models (position, speed,
%         position x speed) for each cell.
%       - LLH_cst: Log likelihood values for constant mean model for each cell.
%       - mapX_SE: Standard error estimates for position tuning curves based on
%         k-fold cross-validation.
%       - tidx: Logical array indicating time indices used for analysis.
%       - ncells_orig: Number of original cells before cell selection.
%
% Example usage:
%   glmsparams = DefineGLMsParams(Nav, Spk);
%   GLMs = GLMAnalysis(Nav, Srep, glmsparams);
%
% See also: glmnet, crossvalPartition, GaussianSmooth1D, computeLLH_poisson
%
% written by J.Fournier 08/2023 for the iBio Summer school


%%
%options for glmnet

%Distribution of the GLM model. For Poisson, default link function is log
%so that means taht the fitted model is of the form y = a0*exp(X*b) where X
%is the design matrix b is the vector of coefficients to be estimated and
%a0 is a constant
family = 'poisson';

%Initializing the options for fitting the GLM with glmnet
options = glmnetSet;

%Number of lambda values used for the regularization path
options.nlambda = 100;

%Whether or not the predictors are randomized. The estimated coefficients
%are always returned on the original scale
options.standardize = true;

%Whether or not to include a constant term in the model (i.e. a0 in the
%equation above)
options.intr = glmsparams.intercept;

%Type of regularization: 1 for lasso, 0 for ridge and elastic net for in
%between
options.alpha = glmsparams.alpha;

%Maximum number of passes over the data for all lambda values
options.maxit = glmsparams.maxit;

%Convergence threshold for the gradient descent (fraction of the null
%deviance)
options.thresh = glmsparams.thresh;
%%

%Copy of the X variable of interest from the input structure (Nav)
X1 = Nav.(glmsparams.Xvariablename);

%If no X2 variable are indicated in glmsparams, we'll just compute
if ~isempty(glmsparams.Yvariablename)
    X2 = Nav.(glmsparams.Yvariablename);
else
    X2 = ones(size(X));
    glmsparams.Ybinedges = 1;
    glmsparams.YsmthNbins = 0;
end

%Selecting time indices over which to compute maps, according to parameters
%defined in glmsparams.subset
tidx = true(size(X));
pnames = fieldnames(glmsparams.subset);
fnames = fieldnames(Nav);
for i = 1:numel(pnames)
    if ismember(pnames{i},fnames)
        fn = str2func(glmsparams.subset.([pnames{i} '_op']));
        tidx = tidx & fn(Nav.(pnames{i}), glmsparams.subset.(pnames{i}));
    elseif ~strcmp(pnames{i}(end-2:end),'_op')
        warning('some fields of glmsparams.subset are not matching fields of Nav')
    end
end

%Selecting time indices where X and Y are in the range of bins
tidx = tidx &...
       X >= glmsparams.Xbinedges(1) & X <= glmsparams.Xbinedges(end) &...
       Y >= glmsparams.Ybinedges(1) & Y <= glmsparams.Ybinedges(end) &...
       ~isnan(X) & ~isnan(Y);
          
%Selecting time and cell indices over which to compute tuning curves
if islogical(glmsparams.cellidx)
    cellidx = find(glmsparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > glmsparams.nspk_th);
else
    cellidx = glmsparams.cellidx(sum(Srep(tidx,glmsparams.cellidx), 1, 'omitnan') > glmsparams.nspk_th);
end


%Subsetting spike trains, X1 and X2.
spikeTrain = Srep(tidx,cellidx);
X1 = X1(tidx);
X2 = X2(tidx);

%Number of cells selected for estimating tuning curves
ncells = size(spikeTrain, 2);

%Number of bins for variable X1
nX1bins = numel(glmsparams.Xbinedges) - 1;

%Number of bins for variable X2
nX2bins = numel(glmsparams.Spdbinedges) - 1;

%Number of data points
ntimepts = size(spikeTrain,1);

%Defining a partition of the data for k-fold cross-validation
cv = crossvalPartition(ntimepts, glmsparams.kfold);

%%
%Discretizing X1 vector according to glmsparams.binedges{1}
X1_discrete = discretize(X1, glmsparams.Xbinedges);

%Discretizing X2 vector according to glmsparams.binedges{2}
X2_discrete = discretize(X2, glmsparams.Spdbinedges);

%%
%Build the matrix of predictors (design matrix) for the first variable.
X1mat = full(sparse(1:ntimepts, X1_discrete, ones(ntimepts,1), ntimepts, nX1bins));

%Removing predictors for which occupancy is below the threshold
validX1 = find(sum(X1mat, 1) > glmsparams.occ_th);
X1mat = X1mat(:,validX1);

%Number of valid predictors for positions
nvalidX1 = numel(validX1);

%%
%Quick note aside: if we wanted to compute place fields by linear
%regression, solving the system of linear equation X1mat * b = spikeTrain,
%we'll do this:
mapXlin = NaN(ncells, nX1bins);
for icell = 1:ncells
    mapXlin(icell,validX1) = spikeTrain(:,icell) \ X1mat;
end

%%
%Computing tuning curves with glmnet for the X1-only model.
%Quick technical detail: Since glmnet in matlab doesn't automatically 
%include lambda = 0 in the regularization path, we have to do it manually 
%by calling glmnet twice so that we can use the regularization path defined by
%glmnet (glmnet is a bit sensitive to the set of lambda values used in the 
%regularization path).
mapX1 = NaN(ncells, nX1bins);
mapX1_cv = NaN(ncells, nX1bins, glmsparams.kfold);
LLH_X1 = NaN(ncells, 1);
parfor icell = 1:ncells
    s = spikeTrain(:,icell);
    opt = options;
    %First call to get the set of lambda values used for the regularization
    %path
    opt.lambda = [];
    gfit = glmnet(X1mat, s, family, opt);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    opt.lambda = [gfit.lambda ; 0];
    gfit = glmnet(X1mat, s, family, opt);
    
    %The tuning curve we're interested in is the one with no
    %regularization, i.e. gfit.beta(:,end). We'll smooth that tuning curve
    %since it usually gives better prediction performances.
    gfitsmth = NaN(nX1bins, 1);
    gfitsmth(validX1) = gfit.beta(:,end);
    gfitsmth = GaussianSmooth(gfitsmth, glmsparams.smthNbins{1});
    gfit.beta(:,end) = gfitsmth(validX1);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp = NaN(1, nX1bins);
    maptemp(validX1) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    
    %Finally returning the result into mapX1. 
    mapX1(icell,:) = maptemp;
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptemp_cv = NaN(nX1bins, glmsparams.kfold);
    for i = 1:glmsparams.kfold
        %Fitting on the training set
        X1mattraining = X1mat(cv.trainsets{i}, :);
        Spktraining = s(cv.trainsets{i});
        opt.lambda = [];
        gfit = glmnet(X1mattraining, Spktraining, family, opt);
        opt.lambda = [gfit.lambda ; 0];
        gfit = glmnet(X1mattraining, Spktraining, family, opt);
        
        %Smothing the fitted coefficient (because it usually gives better
        %predictive performance).
        gfitsmth = NaN(nX1bins, 1);
        gfitsmth(validX1) = gfit.beta(:,end);
        gfitsmth = GaussianSmooth(gfitsmth, glmsparams.smthNbins{1});
        gfit.beta(:,end) = gfitsmth(validX1);
        
        %Constructing the response for the test set
        X1mattest = X1mat(cv.testsets{i}, :);
        ypred(cv.testsets{i}) = exp(gfit.a0(end) + X1mattest * gfit.beta(:,end));
        
        %Converting the fitted coefficients into tuning curve equivalents.
        maptemp_cv(validX1,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    end

    %Saving the cross-validated tuning curves into mapX1_cv.
    mapX1_cv(icell,:,:) = maptemp_cv;

    %Computing the log likelihood of the predicted data under the model
    LLH_X1(icell) = computeLLH_poisson(s, ypred);
end

%%
%Build the matrix of predictors (design matrix) for speed
X2mat = full(sparse(1:ntimepts, X2_discrete, ones(ntimepts,1), ntimepts, nX2bins));

%Removing speed predictors for which occupancy is below the threshold
validX2 = find(sum(X2mat, 1) > glmsparams.occ_th);
X2mat = X2mat(:,validX2);

%Number of valid predictors for speed
nvalidX2 = numel(validX2);

%%
%Computing speed tuning curves with glmnet for the speed-only model. Same
%procedure as for the position-only model
mapX2 = NaN(ncells, nX2bins);
mapS_cv = NaN(ncells, nX2bins, glmsparams.kfold);
LLH_X2 = NaN(ncells, 1);
parfor icell = 1:ncells
    s = spikeTrain(:,icell);
    opt = options;
    %First call to get the set of lambda values used for the regularization
    %path
    opt.lambda = [];
    gfit = glmnet(X2mat, s, family, opt);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    opt.lambda = [gfit.lambda ; 0];
    gfit = glmnet(X2mat, s, family, opt);
    
    %The tuning curve we're interested in is the one with no
    %regularization, i.e. gfit.beta(:,end). We'll smooth that tuning curve
    %since it usually gives better prediction performances.
    gfitsmth = NaN(nX2bins, 1);
    gfitsmth(validX2) = gfit.beta(:,end);
    gfitsmth = GaussianSmooth(gfitsmth, glmsparams.smthNbins{2});
    gfit.beta(:,end) = gfitsmth(validX2);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp = NaN(1, nX2bins);
    maptemp(validX2) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    
    %Finally returning the result into mapX2. 
    mapX2(icell,:) = maptemp;
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptemp_cv = NaN(nX2bins, glmsparams.kfold);
    for i = 1:glmsparams.kfold
        %Fitting on the training set
        X2mattraining = X2mat(cv.trainsets{i}, :);
        Spktraining = s(cv.trainsets{i});
        opt.lambda = [];
        gfit = glmnet(X2mattraining, Spktraining, family, opt);
        opt.lambda = [gfit.lambda ; 0];
        gfit = glmnet(X2mattraining, Spktraining, family, opt);
        
        %Smothing the fitted coefficient (because it usually gives better
        %predictive performance).
        gfitsmth = NaN(nX2bins, 1);
        gfitsmth(validX2) = gfit.beta(:,end);
        gfitsmth = GaussianSmooth(gfitsmth, glmsparams.smthNbins{2});
        gfit.beta(:,end) = gfitsmth(validX2);
        
        %Constructing the response for the test set
        X2mattest = X2mat(cv.testsets{i}, :);
        ypred(cv.testsets{i}) = exp(gfit.a0(end) + X2mattest * gfit.beta(:,end));
        
        %Converting the fitted coefficients into tuning curve equivalents.
        maptemp_cv(validX2,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    end

    %Saving the cross-validated tuning curves into mapS_cv.
    mapX2_cv(icell,:,:) = maptemp_cv;

    %Computing the log likelihood of the predicted data under the model
    LLH_X2(icell) = computeLLH_poisson(s, ypred);
end

%%
%Finally doing the same for the position x speed model

%Build the design matrix for the Position x Speed model
X12mat = cat(2, X1mat, X2mat);

%Indices of predictors related to position
X1idx = 1:nvalidX1;

%Indices of predictors related to speed
X2idx = (nvalidX1 + 1):(nvalidX1 + nvalidX2);

%Computing position and speed tuning curves with glmnet for the position x
%speed model. Same procedure as for the single-variable models
mapX12_X = NaN(ncells, nX1bins);
mapX12_Xcv = NaN(ncells, nX1bins, glmsparams.kfold);
mapX12_S = NaN(ncells, nX2bins);
mapX12_Scv = NaN(ncells, nX2bins, glmsparams.kfold);
LLH_X12 = NaN(ncells, 1);
parfor icell = 1:ncells
    s = spikeTrain(:,icell);
    opt = options;
    %First call to get the set of lambda values used for the regularization
    %path
    opt.lambda = [];
    gfit = glmnet(X12mat, s, family, opt);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    opt.lambda = [gfit.lambda ; 0];
    gfit = glmnet(X12mat, s, family, opt);
    
    %The tuning curve we're interested in is the one with no
    %regularization, i.e. gfit.beta(:,end). We'll smooth that tuning curve
    %since it usually gives better prediction performances.
    gfitsmth_X1 = NaN(nX1bins, 1);
    gfitsmth_X1(validX1) = gfit.beta(Xidx,end);
    gfitsmth_X1 = GaussianSmooth(gfitsmth_X1, glmsparams.smthNbins{1});
    gfit.beta(X1idx,end) = gfitsmth_X(validX1);

    gfitsmth_X2 = NaN(nX2bins, 1);
    gfitsmth_X2(validX2) = gfit.beta(X2idx,end);
    gfitsmth_X2 = GaussianSmooth1D(gfitsmth_X2, glmsparams.SpdsmthNbins);
    gfit.beta(X2idx,end) = gfitsmth_S(validX2);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp_X1 = NaN(1, nX1bins);
    maptemp_X1(validX1) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Xidx,end));
    maptemp_X2 = NaN(1, nX2bins);
    maptemp_X2(validX2) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Sidx,end));
    
    %Finally returning the result into mapX12_X1 and mapX12_X2. 
    mapX12_X1(icell,:) = maptemp_X1;
    mapX12_X2(icell,:) = maptemp_X2;
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptempX1_cv = NaN(nX1bins, glmsparams.kfold);
    maptempX2_cv = NaN(nX2bins, glmsparams.kfold);
    for i = 1:glmsparams.kfold
        %Fitting on the training set
        X12mattraining = X12mat(cv.trainsets{i}, :);
        Spktraining = s(cv.trainsets{i});
        opt.lambda = [];
        gfit = glmnet(X12mattraining, Spktraining, family, opt);
        opt.lambda = [gfit.lambda ; 0];
        gfit = glmnet(X12mattraining, Spktraining, family, opt);
        
        %Smothing the fitted coefficient (because it usually gives better
        %predictive performance).
        gfitsmth_X1 = NaN(nX1bins, 1);
        gfitsmth_X1(validX1) = gfit.beta(X1idx,end);
        gfitsmth_X1 = GaussianSmooth1D(gfitsmth_X1, glmsparams.smthNbins{1});
        gfit.beta(X1idx,end) = gfitsmth_X1(validX1);

        gfitsmth_X2 = NaN(nX2bins, 1);
        gfitsmth_X2(validX2) = gfit.beta(X2idx,end);
        gfitsmth_X2 = GaussianSmooth(gfitsmth_X2, glmsparams.smthNbins{2});
        gfit.beta(X2idx,end) = gfitsmth_X2(validX2);
        
        %Constructing the response for the test set
        X12mattest = X12mat(cv.testsets{i}, :);
        ypred(cv.testsets{i}) = exp(gfit.a0(end) + X12mattest * gfit.beta(:,end));
        
        %Converting the fitted coefficients into tuning curve equivalents.
        maptempX1_cv(validX1,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(X1idx,end));
        maptempX2_cv(validX2,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(X2idx,end));
    end

    %Saving the cross-validated tuning curves into mapX12_X1cv and mapX12_X2cv.
    mapX12_X1cv(icell,:,:) = maptempX_cv;
    mapX12_X2cv(icell,:,:) = maptempS_cv;

    %Computing the log likelihood of the predicted data under the model
    LLH_X12(icell) = computeLLH_poisson(s, ypred);
end

%%
%Computing the log likelihood of the data under a constant mean model
LLH_cst = NaN(ncells, 1);
for icell = 1:ncells
    ypred = NaN(size(spikeTrain(:,icell)));
    for i = 1:glmsparams.kfold
        ypred(cv.testsets{i}) = nanmean(spikeTrain(cv.trainsets{i},icell));
    end
    LLH_cst(icell) = computeLLH_poisson(spikeTrain(:,icell), ypred);
end

%%
%Model selection: which of he tree models estimated above is best for each
%cell

%First identifying which single-variable model is best for each cell
[~, bestsinglemodel] = max([LLH_X1 LLH_X2], [], 2);

bestmodel = zeros(ncells, 1);
pval_X1 = NaN(ncells, 1);
pval_X2 = NaN(ncells, 1);
%First considering cells for which position alone is better than speed
%alone to explain the response.
bestX1cellidx = bestsinglemodel == 1;

%Comparing position-only model to the constant mean model and to the
%full model to test the significance of adding each variable to the model
if sum(bestX1cellidx) > 0
    %Comparing position-only model to the constant mean model to get the
    %significance of adding position to fit the cell response.
    dof = nvalidX1;
    goodidx = LLH_X1 > LLH_cst;
    pval_X1(bestX1cellidx & ~goodidx) = 1;
    [~, pval_X1(bestX1cellidx & goodidx)] = lratiotest(LLH_X1(bestX1cellidx & goodidx), LLH_cst(bestX1cellidx & goodidx), dof);
    
    %Comparing position-only model to the full model to get the
    %significance of adding speed to the position-only model.
    dof = (1 + nvalidX1 + nvalidX2) - (1 + nvalidX1);
    goodidx = LLH_X12 > LLH_X1;
    pval_X2(bestX1cellidx & ~goodidx) = 1;
    [~, pval_X2(bestX1cellidx & goodidx)] = lratiotest(LLH_X12(bestX1cellidx & goodidx), LLH_X1(bestX1cellidx & goodidx), dof);
    
    %If adding speed leads to a significant improvement of the model, we
    %consider the full model as the best model.
    bestmodel(bestX1cellidx & pval_X2 <= glmsparams.pval_th) = 3;
    %Otherwise, we consider the position-only model as the best one.
    bestmodel(bestX1cellidx & pval_X2 > glmsparams.pval_th) = 1;
    %All of that provided that the position-only model is significantly
    %better than the constant mean model in the first place.
    bestmodel(bestX1cellidx & pval_X1 > glmsparams.pval_th) = 0;
end

%Now doing the same thing for cells for which speed alone is better than
%positin alone to explain the response.
bestX2cellidx = bestsinglemodel == 2;

%Comparing speed-only model to the constant mean model and to the
%full model to test the significance of adding each variable to the model
if sum(bestX2cellidx) > 0
    %Comparing speed-only model to the constant mean model to get the
    %significance of adding speed to fit the cell response.
    dof = nvalidX2;
    goodidx = LLH_X2 > LLH_cst;
    pval_X2(bestX2cellidx & ~goodidx) = 1;
    [~, pval_X2(bestX2cellidx & goodidx)] = lratiotest(LLH_X2(bestX2cellidx & goodidx), LLH_cst(bestX2cellidx & goodidx), dof);
    
    %Comparing speed-only model to the full model to get the
    %significance of adding position to the speed-only model.
    dof = (1 + nvalidX1 + nvalidX2) - (1 + nvalidX2);
    goodidx = LLH_X12 > LLH_X2;
    pval_X1(bestX2cellidx & ~goodidx) = 1;
    [~, pval_X1(bestX2cellidx & goodidx)] = lratiotest(LLH_X12(bestX2cellidx & goodidx), LLH_X2(bestX2cellidx & goodidx), dof);
    
    %If adding position leads to a significant improvement of the model, we
    %consider the full model as the best model.
    bestmodel(bestX2cellidx & pval_X1 <= glmsparams.pval_th) = 3;
    %Otherwise, we consider the speed-only model as the best one.
    bestmodel(bestX2cellidx & pval_X1 > glmsparams.pval_th) = 2;
    %All of that provided that the speed-only model is significantly
    %better than the constant mean model in the first place.
    bestmodel(bestX2cellidx & pval_X2 <= glmsparams.pval_th) = 0;
end


%%
%Populating the output structure with the results to be saved
glmsparams.tidx = tidx;
GLMs.glmsparams = glmsparams;
GLMs.Xbincenters = glmsparams.Xbinedges(1:end-1) + glmsparams.Xbinsize / 2;
GLMs.Spdbincenters = glmsparams.Spdbinedges(1:end-1) + glmsparams.Spdbinsize / 2;

ncells_orig = size(Srep, 2);
GLMs.mapX1 = NaN(ncells_orig, nX1bins);
GLMs.mapX1cv = NaN(ncells_orig, nX1bins, glmsparams.kfold);
GLMs.mapX2 = NaN(ncells_orig, nX2bins);
GLMs.map2cv = NaN(ncells_orig, nX2bins, glmsparams.kfold);
GLMs.pval_X1 = NaN(ncells_orig, 1);
GLMs.pval_X2 = NaN(ncells_orig, 1);
GLMs.bestmodel = NaN(ncells_orig, 1);
GLMs.LLH = NaN(ncells_orig, 3);
GLMs.LLH_cst = NaN(ncells_orig, 1);

GLMs.pval_X1(cellidx) = pval_X1;
GLMs.pval_X2(cellidx) = pval_X2;
GLMs.bestmodel(cellidx) = bestmodel;
GLMs.LLH(cellidx,:) = [LLH_X1 LLH_X2 LLH_X12];
GLMs.LLH_cst(cellidx) = LLH_cst;
for icell = 1:ncells
    switch bestmodel(icell)
        case 1
            GLMs.mapX1(cellidx(icell),:) = mapX1(icell,:);
            GLMs.mapX1cv(cellidx(icell),:,:) = mapX1_cv(icell,:,:);
        case 2
            GLMs.mapX2(cellidx(icell),:) = mapX2(icell,:);
            GLMs.map2cv(cellidx(icell),:,:) = mapX2_cv(icell,:,:);
        case 3
            GLMs.mapX1(cellidx(icell),:) = mapX12_X1(icell,:);
            GLMs.mapX2(cellidx(icell),:) = mapX12_X2(icell,:);
            GLMs.mapX1cv(cellidx(icell),:,:) = mapX12_X1cv(icell,:,:);
            GLMs.mapX2cv(cellidx(icell),:,:) = mapX12_X2cv(icell,:,:);
    end
end
%Computing a Jacknife estimate of the standard error
GLMs.mapX1_SE = sqrt((glmsparams.kfold - 1)./glmsparams.kfold * sum((GLMs.mapX1cv - GLMs.mapX1).^2, 3));
GLMs.mapX2_SE = sqrt((glmsparams.kfold - 1)./glmsparams.kfold * sum((GLMs.mapX2cv - GLMs.mapX2).^2, 3));

end