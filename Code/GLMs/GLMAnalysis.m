function GLMs = GLMAnalysis(Nav, Srep, glmsparams)
% GLMs = GLMAnalysis(Nav, Srep, glmsparams)
% Estimates tuning curves for up to two variables using a Poisson GLM
% model. Model comparison is performed by comparing the likelihood of
% held-out data under the model. This function computes GLMs for 
% single-variable models and for the full models, performs k-fold 
% cross-validation, and evaluates model significance.
%
% Inputs:
% - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables
% - Srep: Responses, array where rows correspond to time bins and
%         columns correspond to cells. Each element represents the activity
%         of a cell in a time bin.
% - glmsparams: Structure containing GLM analysis parameters
%
% Outputs:
%   GLMs: Output structure containing the results and analysis of the GLM.
%       - glmsparams: Input parameters for the GLM analysis.
%       - bestmodel: Best model indicator for each cell (0 for constant mean,
%         1 for position, 2 for speed, 3 for position x speed).
%       - LLH: Log likelihood values for different models (position, speed,
%         position x speed) for each cell.
%       - LLH_cst: Log likelihood values for constant mean model for each cell.
%       - tuning: a array of structures with as many elements as the number
%       of input variable (max 2 in this code). tuning has the following
%       fields:
%               * bincenters: Bin centers for the corresponding variable.
%               * map: tuning curves for selected cells along the corresponding variable.
%               * mapcv: Cross-validated tuning curves.
%               * map_SE: Standard error estimates for tuning curves (based on
%                k-fold cross-validation).
%               * pval: P-values indicating significance of adding varaible to the model.
%       
%
% Example usage:
%   Nav = LoaddataNav(loadparams);
%   Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%   Srep = Spk.spikeTrain;
%   glmsparams = SetGLMsParams(Nav, Srep);
%   %change glmsparams here is necessary
%   GLMs = GLMAnalysis(Nav, Srep, glmsparams);
%
% See also: glmnet, crossvalPartition, GaussianSmooth, computeLLH_poisson
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological recordings"


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

%If no X2 variable are indicated in glmsparams, we'll just compute single
%variable model
if ~isempty(glmsparams.variablename{1})
    X1 = Nav.(glmsparams.variablename{1});
else
    X1 = ones(size(Nav.sampleTimes));
    glmsparams.binedges{1} = 1;
    glmsparams.smthNbins{1} = 0;
end

%If no X2 variable are indicated in glmsparams, we'll just compute single
%variable model
if ~isempty(glmsparams.variablename{2})
    X2 = Nav.(glmsparams.variablename{2});
else
    X2 = ones(size(Nav.sampleTimes));
    glmsparams.binedges{2} = 1;
    glmsparams.smthNbins{2} = 0;
end

%Selecting time indices over which to estimate GLMs, according to parameters
%defined in glmsparams.subset
tidx = true(size(X1));
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
       X1 >= glmsparams.binedges{1}(1) & X1 <= glmsparams.binedges{1}(end) &...
       X2 >= glmsparams.binedges{2}(1) & X2 <= glmsparams.binedges{2}(end) &...
       ~isnan(X1) & ~isnan(X2);
          
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
nX1bins = numel(glmsparams.binedges{1}) - 1;

%Number of bins for variable X2
nX2bins = numel(glmsparams.binedges{2}) - 1;

%Number of data points
ntimepts = size(spikeTrain,1);

%Defining a partition of the data for k-fold cross-validation
cv = crossvalPartition(ntimepts, glmsparams.kfold);

%%
%Discretizing X1 vector according to glmsparams.binedges{1}
X1_discrete = discretize(X1, glmsparams.binedges{1});

%Discretizing X2 vector according to glmsparams.binedges{2}
X2_discrete = discretize(X2, glmsparams.binedges{2});

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

    %Constructing the response for the train set
    ytrain = exp(gfit.a0(end) + X1mat * gfit.beta(:,end));
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptemp_cv = NaN(nX1bins, glmsparams.kfold);
    if glmsparams.kfold > 1
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
    else
        %if glmsparams.kfold <= 1, we compute the log likelihood of the
        %fittied data under the model
        LLH_X1(icell) = computeLLH_poisson(s, ytrain);
    end
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
mapX2_cv = NaN(ncells, nX2bins, glmsparams.kfold);
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
    
    %Constructing the fitted response
    ytrain = exp(gfit.a0(end) + X2mat * gfit.beta(:,end));
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptemp_cv = NaN(nX2bins, glmsparams.kfold);
    if glmsparams.kfold > 1
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
    else
        %if glmsparams.kfold <= 1, we compute the log likelihood of the
        %fitted data under the model
        LLH_X2(icell) = computeLLH_poisson(s, ytrain);
    end
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
mapX12_X1 = NaN(ncells, nX1bins);
mapX12_X1cv = NaN(ncells, nX1bins, glmsparams.kfold);
mapX12_X2 = NaN(ncells, nX2bins);
mapX12_X2cv = NaN(ncells, nX2bins, glmsparams.kfold);
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
    gfitsmth_X1(validX1) = gfit.beta(X1idx,end);
    gfitsmth_X1 = GaussianSmooth(gfitsmth_X1, glmsparams.smthNbins{1});
    gfit.beta(X1idx,end) = gfitsmth_X1(validX1);

    gfitsmth_X2 = NaN(nX2bins, 1);
    gfitsmth_X2(validX2) = gfit.beta(X2idx,end);
    gfitsmth_X2 = GaussianSmooth(gfitsmth_X2, glmsparams.smthNbins{2});
    gfit.beta(X2idx,end) = gfitsmth_X2(validX2);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp_X1 = NaN(1, nX1bins);
    maptemp_X1(validX1) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(X1idx,end));
    maptemp_X2 = NaN(1, nX2bins);
    maptemp_X2(validX2) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(X2idx,end));
    
    %Finally returning the result into mapX12_X1 and mapX12_X2. 
    mapX12_X1(icell,:) = maptemp_X1;
    mapX12_X2(icell,:) = maptemp_X2;

    %Reconstructing the fitted data
    ytrain = exp(gfit.a0(end) + X12mat * gfit.beta(:,end));
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptempX1_cv = NaN(nX1bins, glmsparams.kfold);
    maptempX2_cv = NaN(nX2bins, glmsparams.kfold);
    if glmsparams.kfold > 1
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
            gfitsmth_X1 = GaussianSmooth(gfitsmth_X1, glmsparams.smthNbins{1});
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
        mapX12_X1cv(icell,:,:) = maptempX1_cv;
        mapX12_X2cv(icell,:,:) = maptempX2_cv;

        %Computing the log likelihood of the predicted data under the model
        LLH_X12(icell) = computeLLH_poisson(s, ypred);
    else
        %if glmsparams.kfold <= 1, we compute the log likelihood of the
        %fitted data under the model
        LLH_X12(icell) = computeLLH_poisson(s, ytrain);
    end
end

%%
%Computing the log likelihood of the data under a constant mean model
LLH_cst = NaN(ncells, 1);
for icell = 1:ncells
    ypred = NaN(size(spikeTrain(:,icell)));
    if glmsparams.kfold > 1
        for i = 1:glmsparams.kfold
            ypred(cv.testsets{i}) = mean(spikeTrain(cv.trainsets{i},icell),'omitnan');
        end
        LLH_cst(icell) = computeLLH_poisson(spikeTrain(:,icell), ypred);
    else
        LLH_cst(icell) = computeLLH_poisson(spikeTrain(:,icell), mean(spikeTrain(:,icell)) * ones(size(spikeTrain(:,icell))));
    end
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
    if sum(bestX1cellidx & goodidx) > 0
        [~, pval_X1(bestX1cellidx & goodidx)] = lratiotest(LLH_X1(bestX1cellidx & goodidx), LLH_cst(bestX1cellidx & goodidx), dof);
    end

    %Comparing position-only model to the full model to get the
    %significance of adding speed to the position-only model.
    dof = (1 + nvalidX1 + nvalidX2) - (1 + nvalidX1);
    goodidx = LLH_X12 > LLH_X1;
    pval_X2(bestX1cellidx & ~goodidx) = 1;
    if sum(bestX1cellidx & goodidx) > 0
        [~, pval_X2(bestX1cellidx & goodidx)] = lratiotest(LLH_X12(bestX1cellidx & goodidx), LLH_X1(bestX1cellidx & goodidx), dof);
    end

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
    if sum(bestX2cellidx & goodidx) > 0
        [~, pval_X2(bestX2cellidx & goodidx)] = lratiotest(LLH_X2(bestX2cellidx & goodidx), LLH_cst(bestX2cellidx & goodidx), dof);
    end

    %Comparing speed-only model to the full model to get the
    %significance of adding position to the speed-only model.
    dof = (1 + nvalidX1 + nvalidX2) - (1 + nvalidX2);
    goodidx = LLH_X12 > LLH_X2;
    pval_X1(bestX2cellidx & ~goodidx) = 1;
    if sum(bestX2cellidx & goodidx) > 0
        [~, pval_X1(bestX2cellidx & goodidx)] = lratiotest(LLH_X12(bestX2cellidx & goodidx), LLH_X2(bestX2cellidx & goodidx), dof);
    end

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

if nX1bins > 1
    GLMs.tuning(1).bincenters = glmsparams.binedges{1}(1:end-1) + diff(glmsparams.binedges{1}) / 2;
else
    GLMs.tuning(1).bincenters = 1;
end
if nX2bins > 1
    GLMs.tuning(2).bincenters = glmsparams.binedges{2}(1:end-1) + diff(glmsparams.binedges{2}) / 2;
else
    GLMs.tuning(2).bincenters = 1;
end

%Interpolating in case there are infinite values among the bin edges
%(probably at the beginiing or end in order to clamp de signal)
if nX1bins > 1 && sum(isinf(GLMs.tuning(2).bincenters)) > 0
    x1 = 1:nX1bins;
    GLMs.tuning(1).bincenters = interp1(x1(~isinf(GLMs.tuning(1).bincenters)), GLMs.tuning(1).bincenters(~isinf(GLMs.tuning(1).bincenters)), x1, 'linear', 'extrap');
end
if nX2bins > 1 && sum(isinf(GLMs.tuning(2).bincenters)) > 0
    x2 = 1:nX2bins;
    GLMs.tuning(2).bincenters = interp1(x2(~isinf(GLMs.tuning(2).bincenters)), GLMs.tuning(2).bincenters(~isinf(GLMs.tuning(2).bincenters)), x2, 'linear', 'extrap');
end

ncells_orig = size(Srep, 2);
for i = 1:2
    nbins = numel(GLMs.tuning(i).bincenters);
    GLMs.tuning(i).map = NaN(ncells_orig, nbins);
    GLMs.tuning(i).mapcv = NaN(ncells_orig, nbins, glmsparams.kfold);
    GLMs.tuning(i).pval = NaN(ncells_orig, 1);
end
GLMs.bestmodel = NaN(ncells_orig, 1);
GLMs.LLH = NaN(ncells_orig, 3);
GLMs.LLH_cst = NaN(ncells_orig, 1);

GLMs.tuning(1).pval(cellidx) = pval_X1;
GLMs.tuning(2).pval(cellidx) = pval_X2;
GLMs.bestmodel(cellidx) = bestmodel;
GLMs.LLH(cellidx,:) = [LLH_X1 LLH_X2 LLH_X12];
GLMs.LLH_cst(cellidx) = LLH_cst;
for icell = 1:ncells
    switch bestmodel(icell)
        case 1
            GLMs.tuning(1).map(cellidx(icell),:) = mapX1(icell,:);
            GLMs.tuning(1).mapcv(cellidx(icell),:,:) = mapX1_cv(icell,:,:);
        case 2
            GLMs.tuning(2).map(cellidx(icell),:) = mapX2(icell,:);
            GLMs.tuning(2).mapcv(cellidx(icell),:,:) = mapX2_cv(icell,:,:);
        case 3
            GLMs.tuning(1).map(cellidx(icell),:) = mapX12_X1(icell,:);
            GLMs.tuning(2).map(cellidx(icell),:) = mapX12_X2(icell,:);
            GLMs.tuning(1).mapcv(cellidx(icell),:,:) = mapX12_X1cv(icell,:,:);
            GLMs.tuning(2).mapcv(cellidx(icell),:,:) = mapX12_X2cv(icell,:,:);
    end
end
%Computing a Jacknife estimate of the standard error
GLMs.tuning(1).map_SE = sqrt((glmsparams.kfold - 1)./glmsparams.kfold * sum((GLMs.tuning(1).mapcv - GLMs.tuning(1).map).^2, 3));
GLMs.tuning(2).map_SE = sqrt((glmsparams.kfold - 1)./glmsparams.kfold * sum((GLMs.tuning(1).mapcv - GLMs.tuning(1).map).^2, 3));

end