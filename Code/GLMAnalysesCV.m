function GLMs = GLMAnalysesCV(Nav, Srep, glmsparams)
%Estimates tuning curves to position and running speed using a Poisson GLM
%model. Ideally, the model comparison should be done on cross-validated
%data; however, the procedure for fitting GLMs takes too long to perform it
%over kfold. Therefore, we'll do it here by comparing the likelihood of 
%fitted data rather than cross-validated data. (Maybe we should do it
%offline and provide the GLM output structure computed beforehand before
%looking at the results).

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

%Subset of time indices over which GLMs will be computed
tidx = ismember(Nav.Condition, glmsparams.condition) &...
       ismember(Nav.XDir, glmsparams.dir) &...
       Nav.Spd >= glmsparams.spdthreshold &...
       ~isnan(Nav.Xpos);
          
%Selecting time and cell indices over which to compute tuning curves
if islogical(glmsparams.cellidx)
    cellidx = find(glmsparams.cellidx(:)' & sum(Srep(tidx,:), 1, 'omitnan') > glmsparams.nspk_th);
else
    cellidx = glmsparams.cellidx(sum(Srep(tidx,glmsparams.cellidx), 1, 'omitnan') > glmsparams.nspk_th);
end


%Subsetting spike trains, positions and speeds.
spikeTrain = Srep(tidx,cellidx);
Xpos = Nav.Xpos(tidx);
Spd = Nav.Spd(tidx);

%Number of cells selected for estimating tuning curves
ncells = size(spikeTrain, 2);

%Number of position bins
nXbins = numel(glmsparams.Xbinedges) - 1;

%Number of speed bins
nSpdbins = numel(glmsparams.Spdbinedges) - 1;

%Number of data points
ntimepts = size(spikeTrain,1);

%Defining a partition of the data for k-fold cross-validation
cv = crossvalPartition(ntimepts, glmsparams.kfold);

%%
%Discretizing position vector according to glmsparams.Xbinedges
Xpos_discrete = discretize(Xpos, glmsparams.Xbinedges);

%Discretizing speed vector according to glmsparams.Spdbinedges
%Speeds are first clamped within the speed range
Spd(Spd > glmsparams.Spdbinedges(end)) = glmsparams.Spdbinedges(end);
Spd(Spd < glmsparams.Spdbinedges(1)) = glmsparams.Spdbinedges(1);

Spd_discrete = discretize(Spd, glmsparams.Spdbinedges);

%%
%Build the matrix of predictors (design matrix) for positions
Xmat = full(sparse(1:ntimepts, Xpos_discrete, ones(ntimepts,1), ntimepts, nXbins));

%Removing positions for which occupancy is below the threshold
validXpos = find(sum(Xmat, 1) > glmsparams.occ_th);
Xmat = Xmat(:,validXpos);

%Number of valid predictors for positions
nvalidXpos = numel(validXpos);

%%
%Quick note aside: if we wanted to compute place fields by linear
%regression, solving the system of linear equation Xmat * b = spikeTrain,
%we'll do this:
mapXlin = NaN(ncells, nXbins);
for icell = 1:ncells
    mapXlin(icell,validXpos) = spikeTrain(:,icell) \ Xmat;
end

%%
%Computing position tuning curves with glmnet for the position-only model.
%Quick technical detail: Since glmnet in matlab doesn't automatically 
%include lambda = 0 in the regularization path, we have to do it manually 
%by calling glmnet twice so that we use the regularization path defined by
%glmnet (glmnet is a bit sensitive to the set of lambda values used in the 
%regularization path).
mapX = NaN(ncells, nXbins);
mapX_cv = NaN(ncells, nXbins, glmsparams.kfold);
LLH_X = NaN(ncells, 1);
parfor icell = 1:ncells
    s = spikeTrain(:,icell);
    opt = options;
    %First call to get the set of lambda values used for the regularization
    %path
    opt.lambda = [];
    gfit = glmnet(Xmat, s, family, opt);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    opt.lambda = [gfit.lambda ; 0];
    gfit = glmnet(Xmat, s, family, opt);
    
    %The tuning curve we're interested in is the one with no
    %regularization, i.e. gfit.beta(:,end). We'll smooth that tuning curve
    %since it usually gives better prediction performances.
    gfitsmth = NaN(nXbins, 1);
    gfitsmth(validXpos) = gfit.beta(:,end);
    gfitsmth = GaussianSmooth1D(gfitsmth, glmsparams.XsmthNbins);
    gfit.beta(:,end) = gfitsmth(validXpos);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp = NaN(1, nXbins);
    maptemp(validXpos) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    
    %Finally returning the result into mapX. 
    mapX(icell,:) = maptemp;
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptemp_cv = NaN(nXbins, glmsparams.kfold);
    for i = 1:glmsparams.kfold
        %Fitting on the training set
        Xmattraining = Xmat(cv.trainsets{i}, :);
        Spktraining = s(cv.trainsets{i});
        opt.lambda = [];
        gfit = glmnet(Xmattraining, Spktraining, family, opt);
        opt.lambda = [gfit.lambda ; 0];
        gfit = glmnet(Xmattraining, Spktraining, family, opt);
        
        %Smothing the fitted coefficient (because it usually gives better
        %predictive performance).
        gfitsmth = NaN(nXbins, 1);
        gfitsmth(validXpos) = gfit.beta(:,end);
        gfitsmth = GaussianSmooth1D(gfitsmth, glmsparams.XsmthNbins);
        gfit.beta(:,end) = gfitsmth(validXpos);
        
        %Constructing the response for the test set
        Xmattest = Xmat(cv.testsets{i}, :);
        ypred(cv.testsets{i}) = exp(gfit.a0(end) + Xmattest * gfit.beta(:,end));
        
        %Converting the fitted coefficients into tuning curve equivalents.
        maptemp_cv(validXpos,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    end

    %Saving the cross-validated tuning curves into mapX_cv.
    mapX_cv(icell,:,:) = maptemp_cv;

    %Computing the log likelihood of the predicted data under the model
    LLH_X(icell) = computeLLH_poisson(s, ypred);
end

%%
%Build the matrix of predictors (design matrix) for speed
Smat = full(sparse(1:ntimepts, Spd_discrete, ones(ntimepts,1), ntimepts, nSpdbins));

%Removing speed predictors for which occupancy is below the threshold
validSpd = find(sum(Smat, 1) > glmsparams.occ_th);
Smat = Smat(:,validSpd);

%Number of valid predictors for speed
nvalidSpd = numel(validSpd);

%%
%Computing speed tuning curves with glmnet for the speed-only model. Same
%procedure as for the position-only model
mapS = NaN(ncells, nSpdbins);
mapS_cv = NaN(ncells, nSpdbins, glmsparams.kfold);
LLH_S = NaN(ncells, 1);
parfor icell = 1:ncells
    s = spikeTrain(:,icell);
    opt = options;
    %First call to get the set of lambda values used for the regularization
    %path
    opt.lambda = [];
    gfit = glmnet(Smat, s, family, opt);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    opt.lambda = [gfit.lambda ; 0];
    gfit = glmnet(Smat, s, family, opt);
    
    %The tuning curve we're interested in is the one with no
    %regularization, i.e. gfit.beta(:,end). We'll smooth that tuning curve
    %since it usually gives better prediction performances.
    gfitsmth = NaN(nSpdbins, 1);
    gfitsmth(validSpd) = gfit.beta(:,end);
    gfitsmth = GaussianSmooth1D(gfitsmth, glmsparams.SpdsmthNbins);
    gfit.beta(:,end) = gfitsmth(validSpd);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp = NaN(1, nSpdbins);
    maptemp(validSpd) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    
    %Finally returning the result into mapS. 
    mapS(icell,:) = maptemp;
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptemp_cv = NaN(nSpdbins, glmsparams.kfold);
    for i = 1:glmsparams.kfold
        %Fitting on the training set
        Smattraining = Smat(cv.trainsets{i}, :);
        Spktraining = s(cv.trainsets{i});
        opt.lambda = [];
        gfit = glmnet(Smattraining, Spktraining, family, opt);
        opt.lambda = [gfit.lambda ; 0];
        gfit = glmnet(Smattraining, Spktraining, family, opt);
        
        %Smothing the fitted coefficient (because it usually gives better
        %predictive performance).
        gfitsmth = NaN(nSpdbins, 1);
        gfitsmth(validSpd) = gfit.beta(:,end);
        gfitsmth = GaussianSmooth1D(gfitsmth, glmsparams.SpdsmthNbins);
        gfit.beta(:,end) = gfitsmth(validSpd);
        
        %Constructing the response for the test set
        Smattest = Smat(cv.testsets{i}, :);
        ypred(cv.testsets{i}) = exp(gfit.a0(end) + Smattest * gfit.beta(:,end));
        
        %Converting the fitted coefficients into tuning curve equivalents.
        maptemp_cv(validSpd,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    end

    %Saving the cross-validated tuning curves into mapS_cv.
    mapS_cv(icell,:,:) = maptemp_cv;

    %Computing the log likelihood of the predicted data under the model
    LLH_S(icell) = computeLLH_poisson(s, ypred);
end

%%
%Finally doing the same for the position x speed model

%Build the design matrix for the Position x Speed model
XSmat = cat(2, Xmat, Smat);

%Indices of predictors related to position
Xidx = 1:nvalidXpos;

%Indices of predictors related to speed
Sidx = (nvalidXpos + 1):(nvalidXpos + nvalidSpd);

%Computing position and speed tuning curves with glmnet for the position x
%speed model. Same procedure as for the single-variable models
mapXS_X = NaN(ncells, nXbins);
mapXS_Xcv = NaN(ncells, nXbins, glmsparams.kfold);
mapXS_S = NaN(ncells, nSpdbins);
mapXS_Scv = NaN(ncells, nSpdbins, glmsparams.kfold);
LLH_XS = NaN(ncells, 1);
parfor icell = 1:ncells
    s = spikeTrain(:,icell);
    opt = options;
    %First call to get the set of lambda values used for the regularization
    %path
    opt.lambda = [];
    gfit = glmnet(XSmat, s, family, opt);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    opt.lambda = [gfit.lambda ; 0];
    gfit = glmnet(XSmat, s, family, opt);
    
    %The tuning curve we're interested in is the one with no
    %regularization, i.e. gfit.beta(:,end). We'll smooth that tuning curve
    %since it usually gives better prediction performances.
    gfitsmth_X = NaN(nXbins, 1);
    gfitsmth_X(validXpos) = gfit.beta(Xidx,end);
    gfitsmth_X = GaussianSmooth1D(gfitsmth_X, glmsparams.XsmthNbins);
    gfit.beta(Xidx,end) = gfitsmth_X(validXpos);

    gfitsmth_S = NaN(nSpdbins, 1);
    gfitsmth_S(validSpd) = gfit.beta(Sidx,end);
    gfitsmth_S = GaussianSmooth1D(gfitsmth_S, glmsparams.SpdsmthNbins);
    gfit.beta(Sidx,end) = gfitsmth_S(validSpd);

    %Converting the fitted coefficient into tuning curves equivalent
    %expressed in spike / s. This is done in a temporary map for indexing
    %purposes due to the parfor loop.
    maptemp_X = NaN(1, nXbins);
    maptemp_X(validXpos) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Xidx,end));
    maptemp_S = NaN(1, nSpdbins);
    maptemp_S(validSpd) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Sidx,end));
    
    %Finally returning the result into mapXS_X and mapXS_S. 
    mapXS_X(icell,:) = maptemp_X;
    mapXS_S(icell,:) = maptemp_S;
    
    %Doing the same thing for all k-fold partitions in order to compute the
    %predicted response
    ypred = NaN(size(s));
    maptempX_cv = NaN(nXbins, glmsparams.kfold);
    maptempS_cv = NaN(nSpdbins, glmsparams.kfold);
    for i = 1:glmsparams.kfold
        %Fitting on the training set
        XSmattraining = XSmat(cv.trainsets{i}, :);
        Spktraining = s(cv.trainsets{i});
        opt.lambda = [];
        gfit = glmnet(XSmattraining, Spktraining, family, opt);
        opt.lambda = [gfit.lambda ; 0];
        gfit = glmnet(XSmattraining, Spktraining, family, opt);
        
        %Smothing the fitted coefficient (because it usually gives better
        %predictive performance).
        gfitsmth_X = NaN(nXbins, 1);
        gfitsmth_X(validXpos) = gfit.beta(Xidx,end);
        gfitsmth_X = GaussianSmooth1D(gfitsmth_X, glmsparams.XsmthNbins);
        gfit.beta(Xidx,end) = gfitsmth_X(validXpos);

        gfitsmth_S = NaN(nSpdbins, 1);
        gfitsmth_S(validSpd) = gfit.beta(Sidx,end);
        gfitsmth_S = GaussianSmooth1D(gfitsmth_S, glmsparams.SpdsmthNbins);
        gfit.beta(Sidx,end) = gfitsmth_S(validSpd);
        
        %Constructing the response for the test set
        XSmattest = XSmat(cv.testsets{i}, :);
        ypred(cv.testsets{i}) = exp(gfit.a0(end) + XSmattest * gfit.beta(:,end));
        
        %Converting the fitted coefficients into tuning curve equivalents.
        maptempX_cv(validXpos,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Xidx,end));
        maptempS_cv(validSpd,i) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Sidx,end));
    end

    %Saving the cross-validated tuning curves into mapXS_Xcv and mapXS_Scv.
    mapXS_Xcv(icell,:,:) = maptempX_cv;
    mapXS_Scv(icell,:,:) = maptempS_cv;

    %Computing the log likelihood of the predicted data under the model
    LLH_XS(icell) = computeLLH_poisson(s, ypred);
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
[~, bestsinglemodel] = max([LLH_X LLH_S], [], 2);

bestmodel = zeros(ncells, 1);
pval_X = NaN(ncells, 1);
pval_S = NaN(ncells, 1);
%First considering cells for which position alone is better than speed
%alone to explain the response.
bestXcellidx = bestsinglemodel == 1;

%Comparing position-only model to the constant mean model and to the
%full model to test the significance of adding each variable to the model
if sum(bestXcellidx) > 0
    %Comparing position-only model to the constant mean model to get the
    %significance of adding position to fit the cell response.
    dof = nvalidXpos;
    goodidx = LLH_X > LLH_cst;
    pval_X(bestXcellidx & ~goodidx) = 1;
    [~, pval_X(bestXcellidx & goodidx)] = lratiotest(LLH_X(bestXcellidx & goodidx), LLH_cst(bestXcellidx & goodidx), dof);
    
    %Comparing position-only model to the full model to get the
    %significance of adding speed to the position-only model.
    dof = (1 + nvalidXpos + nvalidSpd) - (1 + nvalidXpos);
    goodidx = LLH_XS > LLH_X;
    pval_S(bestXcellidx & ~goodidx) = 1;
    [~, pval_S(bestXcellidx & goodidx)] = lratiotest(LLH_XS(bestXcellidx & goodidx), LLH_X(bestXcellidx & goodidx), dof);
    
    %If adding speed leads to a significant improvement of the model, we
    %consider the full model as the best model.
    bestmodel(bestXcellidx & pval_S <= glmsparams.pval_th) = 3;
    %Otherwise, we consider the position-only model as the best one.
    bestmodel(bestXcellidx & pval_S > glmsparams.pval_th) = 1;
    %All of that provided that the position-only model is significantly
    %better than the constant mean model in the first place.
    bestmodel(bestXcellidx & pval_X > glmsparams.pval_th) = 0;
end

%Now doing the same thing for cells for which speed alone is better than
%positin alone to explain the response.
bestScellidx = bestsinglemodel == 2;

%Comparing speed-only model to the constant mean model and to the
%full model to test the significance of adding each variable to the model
if sum(bestScellidx) > 0
    %Comparing speed-only model to the constant mean model to get the
    %significance of adding speed to fit the cell response.
    dof = nvalidSpd;
    goodidx = LLH_S > LLH_cst;
    pval_S(bestScellidx & ~goodidx) = 1;
    [~, pval_S(bestScellidx & goodidx)] = lratiotest(LLH_S(bestScellidx & goodidx), LLH_cst(bestScellidx & goodidx), dof);
    
    %Comparing speed-only model to the full model to get the
    %significance of adding position to the speed-only model.
    dof = (1 + nvalidXpos + nvalidSpd) - (1 + nvalidSpd);
    goodidx = LLH_XS > LLH_S;
    pval_X(bestScellidx & ~goodidx) = 1;
    [~, pval_X(bestScellidx & goodidx)] = lratiotest(LLH_XS(bestScellidx & goodidx), LLH_S(bestScellidx & goodidx), dof);
    
    %If adding position leads to a significant improvement of the model, we
    %consider the full model as the best model.
    bestmodel(bestScellidx & pval_X <= glmsparams.pval_th) = 3;
    %Otherwise, we consider the speed-only model as the best one.
    bestmodel(bestScellidx & pval_X > glmsparams.pval_th) = 2;
    %All of that provided that the speed-only model is significantly
    %better than the constant mean model in the first place.
    bestmodel(bestScellidx & pval_S <= glmsparams.pval_th) = 0;
end


%%
%Populating the output structure with the results to be saved
glmsparams.tidx = tidx;
GLMs.glmsparams = glmsparams;
GLMs.Xbincenters = glmsparams.Xbinedges(1:end-1) + glmsparams.Xbinsize / 2;
GLMs.Spdbincenters = glmsparams.Spdbinedges(1:end-1) + glmsparams.Spdbinsize / 2;

ncells_orig = size(Srep, 2);
GLMs.mapX = NaN(ncells_orig, nXbins);
GLMs.mapXcv = NaN(ncells_orig, nXbins, glmsparams.kfold);
GLMs.mapS = NaN(ncells_orig, nSpdbins);
GLMs.mapScv = NaN(ncells_orig, nSpdbins, glmsparams.kfold);
GLMs.pval_X = NaN(ncells_orig, 1);
GLMs.pval_S = NaN(ncells_orig, 1);
GLMs.bestmodel = NaN(ncells_orig, 1);
GLMs.LLH = NaN(ncells_orig, 3);
GLMs.LLH_cst = NaN(ncells_orig, 1);

GLMs.pval_X(cellidx) = pval_X;
GLMs.pval_S(cellidx) = pval_S;
GLMs.bestmodel(cellidx) = bestmodel;
GLMs.LLH(cellidx,:) = [LLH_X LLH_S LLH_XS];
GLMs.LLH_cst(cellidx) = LLH_cst;
for icell = 1:ncells
    switch bestmodel(icell)
        case 1
            GLMs.mapX(cellidx(icell),:) = mapX(icell,:);
            GLMs.mapXcv(cellidx(icell),:,:) = mapX_cv(icell,:,:);
        case 2
            GLMs.mapS(cellidx(icell),:) = mapS(icell,:);
            GLMs.mapScv(cellidx(icell),:,:) = mapS_cv(icell,:,:);
        case 3
            GLMs.mapX(cellidx(icell),:) = mapXS_X(icell,:);
            GLMs.mapS(cellidx(icell),:) = mapXS_S(icell,:);
            GLMs.mapXcv(cellidx(icell),:,:) = mapXS_Xcv(icell,:,:);
            GLMs.mapScv(cellidx(icell),:,:) = mapXS_Scv(icell,:,:);
    end
end
%Computing a Jacknife estimate of the standard error
GLMs.mapX_SE = sqrt((glmsparams.kfold - 1)./glmsparams.kfold * sum((GLMs.mapXcv - GLMs.mapX).^2, 3));

end