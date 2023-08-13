function glmsparams = SetGLMsParams(Nav,Spk)
% Define a set of parameters needed to estimate Generalized Linear Models (GLMs)
% using GLMAnalyses or GLMAnalysesCV functions.
%
% Inputs:
%   Nav: Structure containing navigation-related data (e.g., spatial positions, speeds)
%   Spk: Structure containing spike train data
%
% Outputs:
%   glmsparams: Structure containing the defined GLM parameters
%
% The output structure glmsparams contains the following fields:
%    - subset: a structure with parameters to subset the data across time,
%      according to values taken by certain fields of Nav.
%      Fields in subset have names matching the fields in Nav that we want 
%      to apply the condition on. The values and operators defined in 
%      subset are used to subset the data according to the fields of Nav.
%      data. For instance, if we want to subset data corresponding to
%      Nav.Condition = 1 or 3 and Nav.Spd >= 5, glmsparams.subset should be
%      defined as:
%      glmsparams.subset.Condition = [1 3];
%      glmsparams.subset.Condition_op = 'ismember';
%      glmsparams.subset.Spd = 5;
%      glmsparams.subset.Spd_op = '>=';
%   - cellidx: Logical array indicating the subset of cells for which GLMs will be computed.
%   - sampleRate: Sampling rate of the data.
%   - scalingFactor: Scaling factor on the response data.
%   - XsmthNbins: Size of the Gaussian window for smoothing along X (in bins).
%   - Xbinedges: Edges of position bins used to discretize positions.
%   - YsmthNbins: Size of the Gaussian window for smoothing speed tuning curves (in bins).
%   - Ybinedges: Edges of speed bins used to discretize speed.
%   - occ_th: Occupancy threshold above which predictors are included in the GLM estimate.
%   - nspk_th: Minimal number of spikes to consider a cell for GLM estimation.
%   - kfold: Number of folds for cross-validation.
%   - intercept: Flag indicating whether the GLM model should include a constant term.
%   - alpha: Type of regularization in glmnet: 1 for lasso, 0 for ridge, elastic net for in-between.
%   - maxit: Maximal number of passes over the data before reaching convergence.
%   - thresh: Convergence threshold for coordinate descent, as a fraction of the null deviance.
%   - pval_th: P-value threshold for considering a predictor as significantly contributing to the cell response.
%
% USAGE:
%   glmsparams = SetGLMsParams(Nav, Spk);
%
% Note: Modify the values of the fields as needed for your analysis.
%
% written by J.Fournier 08/2023 for the iBio Summer school


%Conditions over the fields of Nav for which place fields will be estimated
%glmsparams.subset should be a structure where fields have names of the 
%fields of Nav to which the condition should apply to.
glmsparams.subset = [];

%For instance, for the example data set, we define the following fields
glmsparams.subset.Condition = [1 3 5];
glmsparams.subset.Condition_op = 'ismember';

glmsparams.subset.XDir = [1];
glmsparams.subset.XDir_op = 'ismember';

glmsparams.subset.laptype = [-1 0 1];
glmsparams.subset.laptype_op = 'ismember';

glmsparams.subset.Spd =  2.5;
glmsparams.subset.Spd_op = '>=';

%Subset of cells for which GLMs will be computed
glmsparams.cellidx = true(1, size(Spk.spikeTrain, 2));

%Sampling rate of the data
glmsparams.sampleRate = 1 / nanmean(diff(Nav.sampleTimes));

%Scaling factor on the response data (default is 1 / samplingRate so that
%spiking data are returned in spike / s)
glmsparams.scalingFactor = 1 / glmsparams.sampleRate;

%Name of the field of Nav that'll be used as the first predictor. Default
%is Xpos
glmsparams.variablename{1} = 'Xpos';

%Size of the gaussian window for smoothing tuning curves along the first
%predictor (in bins).
glmsparams.smthNbins{1} = 1;

%Edges of bins used to discretize the first predictor.
glmsparams.binedges{1} = 0 : 4: 100;

%Name of the field of Nav that'll be used as the second predictor. Default
%is Spd
glmsparams.variablename{2} = 'Spd';

%Size of the gaussian window for smoothing tuning curves along the second
%predictor (in bins).
glmsparams.smthNbins{2} = 1;

%Edges of bins used to discretize the second predictor
glmsparams.binedges{2} = [2.5:5:47.5 inf];

%Occupancy threshold above which predictors are included in the GLM
%estimate.
glmsparams.occ_th = 0;

%Minimal number of spikes to consider a cell for GLM estimation
glmsparams.nspk_th = 10;

%Number of folds to consider for cross-validation
glmsparams.kfold = 10;

%Whether or not the GLM model should include a constant term
glmsparams.intercept = true;

%Type of regularization in glmnet: 1 for lasso, 0 for ridge, elsticnet for
%everything in between
glmsparams.alpha = 1;

%Maximal number of passes over the data before reaching convergence
glmsparams.maxit = 10^3;

%Convergence threshold for coordinate descent, expressed as a fraction of
%the null deviance.
glmsparams.thresh = 10^-3;

%P-value threshold for considering a predictor as significantly
%contributing to the cell response.
glmsparams.pval_th = 0.05;
end