function GLMs = GLMAnalyses(Nav, Srep, glmsparams)
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

%Subsetting spike trains, positions and speeds
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
LLH_X = NaN(ncells, 1);
for icell = 1:ncells
    %First call to get the set of lambda values used for the regularization
    %path
    options.lambda = [];
    gfit = glmnet(Xmat, spikeTrain(:,icell), family, options);
    
    %Calling glmnet again after adding lambda = 0 as an additional
    %regularization value.
    options.lambda = [gfit.lambda ; 0];
    gfit = glmnet(Xmat, spikeTrain(:,icell), family, options);

    %Computing the output of the full model without regularization.
    yfit = exp(gfit.a0(end) + Xmat * gfit.beta(:,end));

    %Alternatively, we could call the glmnetPredict function to get the
    %reconstructed response.
    %yfit = glmnetPredict(gfit, Xmat, 0, 'response');

    %Computing the log likelihood of the fitted data under the model
    LLH_X(icell) = computeLLH_poisson(spikeTrain(:,icell), yfit);

    %Getting the tuning curve estimated for the last lambda value (i.e.
    %lambda = 0)
    mapX(icell,validXpos) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));

    %Smoothing the tuning curve with a Gaussian
    mapX(icell,:) = GaussianSmooth1D(mapX(icell,:), glmsparams.XsmthNbins);
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
LLH_S = NaN(ncells, 1);
for icell = 1:ncells
    %Estimating speed tuning curve with glmnet
    options.lambda = [];
    gfit = glmnet(Smat, spikeTrain(:,icell), family, options);
    options.lambda = [gfit.lambda ; 0];
    gfit = glmnet(Smat, spikeTrain(:,icell), family, options);

    %Computing the output of the full model without regularization.
    yfit = exp(gfit.a0(end) + Smat * gfit.beta(:,end));

    %Alternatively, we could call the glmnetPredict function to get the
    %reconstructed response.
    %yfit = glmnetPredict(gfit, Smat, 0, 'response');

    %Computing the log likelihood of the fitted data under the model
    LLH_S(icell) = computeLLH_poisson(spikeTrain(:,icell), yfit);
    
    %Getting speed tuning curve with no regularization
    mapS(icell,validSpd) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(:,end));
    
    %Smoothing the tuning curve with a Gaussian
    mapS(icell,:) = GaussianSmooth1D(mapS(icell,:), glmsparams.SpdsmthNbins);
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
mapXS_S = NaN(ncells, nSpdbins);
LLH_XS = NaN(ncells, 1);
for icell = 1:ncells
    %Estimating the full model with glmnet
    options.lambda = [];
    gfit = glmnet(XSmat, spikeTrain(:,icell), family, options);
    options.lambda = [gfit.lambda ; 0];
    gfit = glmnet(XSmat, spikeTrain(:,icell), family, options);

    %Computing the output of the full model without regularization.
    yfit = exp(gfit.a0(end) + XSmat * gfit.beta(:,end));

    %Alternatively, we could call the glmnetPredict function to get the
    %reconstructed response.
    %yfit = glmnetPredict(gfit, XSmat, 0, 'response');
    
    %Computing the log likelihood of the fitted data under the model
    LLH_XS(icell) = computeLLH_poisson(spikeTrain(:,icell), yfit);
    
    %Getting the position ansd speed tuning curves with no regularization.
    %To get curves equivalent to tuning curves, we marginalize out the
    %effect of other variables by adding their mean to the tuning curve
    %coefficients. Then because it is a Poisson distribution that we're 
    %using we need to invert the log link function (taking the exponential)
    %and we convert it to spike / s by multiplying with the sampling rate.
    mapXS_X(icell,validXpos) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Xidx,end) + mean(gfit.beta(Sidx,end)));
    mapXS_S(icell,validSpd) = 1 / glmsparams.scalingFactor * exp(gfit.a0(end) + gfit.beta(Sidx,end) + mean(gfit.beta(Xidx,end)));
    
    %Smoothing the tuning curves with a Gaussian
    mapXS_X(icell,:) = GaussianSmooth1D(mapXS_X(icell,:), glmsparams.XsmthNbins);
    mapXS_S(icell,:) = GaussianSmooth1D(mapXS_S(icell,:), glmsparams.SpdsmthNbins);
end

%%
%Computing the log likelihood of the data under a constant mean model
LLH_cst = NaN(ncells, 1);
for icell = 1:ncells
    yfit = nanmean(spikeTrain(:,icell));
    LLH_cst(icell) = computeLLH_poisson(spikeTrain(:,icell), yfit);
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
    [~, pval_X(bestXcellidx)] = lratiotest(LLH_X(bestXcellidx), LLH_cst(bestXcellidx), dof);
    
    %Comparing position-only model to the full model to get the
    %significance of adding speed to the position-only model.
    dof = (1 + nvalidXpos + nvalidSpd) - (1 + nvalidXpos);
    [~, pval_S(bestXcellidx)] = lratiotest(LLH_XS(bestXcellidx), LLH_X(bestXcellidx), dof);
    
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
    [~, pval_S(bestScellidx)] = lratiotest(LLH_S(bestScellidx), LLH_cst(bestScellidx), dof);
    
    %Comparing speed-only model to the full model to get the
    %significance of adding position to the speed-only model.
    dof = (1 + nvalidXpos + nvalidSpd) - (1 + nvalidSpd);
    [~, pval_X(bestScellidx)] = lratiotest(LLH_XS(bestScellidx), LLH_S(bestScellidx), dof);
    
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
ncells_orig = size(Srep, 2);
GLMs.mapX = NaN(ncells_orig, nXbins);
GLMs.mapS = NaN(ncells_orig, nSpdbins);
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
        case 2
            GLMs.mapS(cellidx(icell),:) = mapS(icell,:);
        case 3
            GLMs.mapX(cellidx(icell),:) = mapXS_X(icell,:);
            GLMs.mapS(cellidx(icell),:) = mapXS_S(icell,:);
    end
end

GLMs.Xbincenters = glmsparams.Xbinedges(1:end-1) + glmsparams.Xbinsize / 2;
GLMs.Spdbincenters = glmsparams.Spdbinedges(1:end-1) + glmsparams.Spdbinsize / 2;
GLMs.glmsparams = glmsparams;

end