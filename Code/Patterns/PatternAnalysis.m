function Pat = PatternAnalysis(Nav, Srep, pattparams)
% PatternAnalysis - Identify cell assemblies in Srep, based on independent component analysis.
%
%   Pat = PatternAnalysis(Nav, Srep, pattparams) identifies cell assemblies
%   in spike response data using independent component analysis (ICA).
%
% INPUTS:
% - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables
% - Srep: Spike train data for each neuron (timepoints x neurons).
% - pattparams: Structure with pattern analysis parameters. See
%   SetPattParams for more details.
%
% OUTPUT:
% - Pat: Structure containing pattern analysis results with the following fields:
%   - pattparams: Structure of parameters used for pattern detection.
%   - weights: Weight matrix representing the contribution of each neuron to each assembly pattern.
%   - cellAssemblies: Logical matrix indicating which neurons belong to each detected cell assembly.
%   - Sparity: Sparsity measure of each assembly pattern (between 0 and 1).
%   - strength: Expression strength of each assembly pattern over time.
%   - activation: Events of assembly pattern activation.
%
% SEE ALSO:
%   SetPattParams, fastICA, GaussianSmooth, findpeaks
%
% Method based on the description in the supplementary information of
% van de Ven et al., 2016 in Neuron.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    pattparams = SetPattParams(Nav, Srep);
%    %change parameters in pattparams here if needed. For instance:
%    %pattparams.subset.Condition = [3]; %identification of patterns on specific condition
%    %pattparams.NoiseCov = true;%subtract signal correlations from covariance matrix
%    %pattparams.Marcenko = false;%use shuffling procedure instead of Marcenko-Pastur law to select principal components
%    Pat = PatternAnalysis(Nav, Srep, pattparams)
%
% Written by J Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"

%%
%Smoothing the spike train over the coincidence window to get spike counts
decbinwin = 2 * floor(0.5 * pattparams.timewin * pattparams.sampleRate) + 1;
if decbinwin > 1
    spkCount = zeros(size(Srep));
    for icell = 1:size(Srep,2)
        spkCount(:,icell) = smooth(Srep(:,icell), decbinwin) * decbinwin;
    end
else
    spkCount = Srep;
end

%Selecting time indices over which to the patterns will be identified, 
% according to parameters defined in pattparams.subset.
traintidx = true(size(Nav.sampleTimes));
pnames = fieldnames(pattparams.subset);
fnames = fieldnames(Nav);
for i = 1:numel(pnames)
    if ismember(pnames{i},fnames)
        fn = str2func(pattparams.subset.([pnames{i} '_op']));
        traintidx = traintidx & fn(Nav.(pnames{i}), pattparams.subset.(pnames{i}));
    elseif ~strcmp(pnames{i}(end-2:end),'_op')
        warning('some fields of pattparams.subset are not matching fields of Nav')
    end
end

%Selecting cell indices over which to detect assemblies
if islogical(pattparams.cellidx)
    cellidx = find(pattparams.cellidx(:)' & sum(Srep(traintidx,:), 1, 'omitnan') > pattparams.nspk_th);
else
    cellidx = pattparams.cellidx(sum(Srep(traintidx,pattparams.cellidx), 1, 'omitnan') > pattparams.nspk_th);
end

%Subsetting spike trains across cells.
spikeTrain = Srep(:,cellidx);
spkCount = spkCount(:,cellidx);
ncells = size(spkCount, 2);

%%
%z-scoring spike counts
mu = mean(spkCount(traintidx,:), 1);
sd = std(spkCount(traintidx,:), 1);
Z = (spkCount(traintidx,:) - mu) ./ sd;

%Computing the covariance matrix
ntimepts = sum(traintidx);
Corrmat = Z' * Z / ntimepts;

%Computing the principal components (ie the eigenvectors of the covariance
%matrix) from the covariance matrix
[P,D] = eig(Corrmat);
eigval = diag(D);

%To select the most significant eigenvalue, we define a threshold on the 
%eigenvalues either using the Marcenko-Pastur law or a shuffle
%distribution
if pattparams.Marcenko
    %Selecting only principal components with an eigenvalue higher than
    %expected according to the Marčenko-Pastur law (see van de Ven et al., 2016).
    eigval_th = (1 + ncells / ntimepts)^2;
end
if ~pattparams.Marcenko || pattparams.NoiseCov
    %The shuffling procedure consists in establishing a distribution of
    %eigenvalues obtained after shuffling time points within bins of the
    %varaibles provided in pattparams.variablenames.
    %We first start by discretizing these variables.
    nVars = numel(pattparams.variablenames);
    if nVars > 0
        vars_discrete = cell(1,nVars);
        sz = cellfun(@numel,pattparams.binedges) - 1;
        for i = 1:nVars
            vars_discrete{i} = discretize(Nav.(pattparams.variablenames{i})(traintidx), pattparams.binedges{i});
        end
        %linearizing the indices across all variables
        varlin_discrete = sub2ind(sz, vars_discrete{:});
        nbins = prod(sz);
    else
        varlin_discrete = NaN(ntimepts, 1);
        nbins = 0;
    end
    %Initializing the random number generator for reproducibility purposes
    eigvalshf = NaN(size(Z,2), pattparams.nShuffle);
    CorrmatSignal = 0;
    parfor ishf = 1:pattparams.nShuffle
        s = RandStream('mt19937ar','Seed',ishf);
        Zshf = NaN(size(Z));
        if nVars > 0
            for k = 1:nbins
                idx = find(varlin_discrete == k);
                for icell = 1:size(Z, 2)
                    idxshf = idx(randperm(s,numel(idx)));
                    Zshf(idx,icell) = Z(idxshf,icell);
                end
            end
        else
            %Circularly shifting each spike count by a random amount > 1 second
            %if no variables have been provided for the shuffle controls
            for icell = 1:size(Z, 2)
                tshift = randi(s,ntimepts - 2 * pattparams.sampleRate) + 1 * pattparams.sampleRate;
                Zshf(:,icell) = circshift(Z(:,icell), tshift);
            end
        end

        %Covariance matrix after shuffling
        CorrmatShf = Zshf' * Zshf / ntimepts;
        
        %Eigenvector decomposition from shuffled covariance
        [~,Dshf] = eig(CorrmatShf);
        eigvalshf(:,ishf) = diag(Dshf);

        %Averaging shuffled covariance matrix, as an estimate of the signal
        %covariance.
        CorrmatSignal = CorrmatSignal + CorrmatShf / pattparams.nShuffle;
    end
    if ~pattparams.Marcenko
        eigval_th = min(eigval(sum(eigval < eigvalshf, 2) / pattparams.nShuffle <= pattparams.pvalshf_th));
    end
end

%Indices of eigenvectors to keep.
sigeigidx = find(eigval > eigval_th, 1, 'first'):numel(eigval);

%pattparams.NoiseCov is true, we remove the signal covariance from the
%observed covariance to eventually obtain assemblies in which coincidental
%firing can't be explained by shared selectivity to the variables indicated
%in pattparams.variablenames.
if pattparams.NoiseCov
    Corrmat(~eye(size(Corrmat))) = Corrmat(~eye(size(Corrmat))) - CorrmatSignal(~eye(size(Corrmat)));
    [P,~] = eig(Corrmat);
end

%Selecting eigenvalues that are above the threshold.
Psign = P(:, sigeigidx);

%Projecting the data onto the selected principal components
Zproj = Z * Psign;

%%
%Doing the ICA on the data projected onto the PCA subspace. W is the mixing
%matrix whose rows contain the weights to reconstruct each independent
%component.
Ncomp = size(Zproj, 2);
[~, W, ~, ~] = fastICA(Zproj', Ncomp, 'kurtosis', 0);

%Expressing the independent components in the original space, ie with a
%weight for each neuron describing its contribution to the corresponding
%component
Vpatt = Psign * W';

%Scaling the weights in Vpatt to unit norm.
Vpatt = Vpatt ./ sqrt(sum(Vpatt.^2, 1));

%By convention, enforcing that the largest element of each vector has a
%positive sign.
[p, d] = size(Vpatt);
[~, maxind] = max(abs(Vpatt), [], 1);
colsign = sign(Vpatt(maxind + (0:p:(d-1)*p)));
Vpatt = bsxfun(@times, Vpatt, colsign);

%%
%Computing the sparsity of each assembly pattern. 1 if only one neuron
%contributes the pattern's weights; 0 if all neurons contribute equally.
PSparity = ((sqrt(ncells) - sum(abs(Vpatt), 1)) ./ (sqrt(ncells) - 1));

%Reordering the patterns from highest to lowest sparsity.
[PSparity, idx] = sort(PSparity, 'descend');
Vpatt = Vpatt(:,idx);

%Defining cell assemblies as neurons whose weight exceed the mean weight by
%two s.d. within a given pattern
m = mean(Vpatt, 1);
s = std(Vpatt, 1);
cellAss = Vpatt > m + 2 * s;

%%
%Computing the expression of each assembly pattern over time
%First smoothing the zscored spike trains by a Gaussian window having the 
%same s.d. as a square time window (??)
Zall = spikeTrain;
smthNbins = (pattparams.timewin * pattparams.sampleRate) / sqrt(12);
Zall = GaussianSmooth(Zall, [smthNbins 0]);

%z-scoring the smoothed spiketrains
mu = mean(Zall, 1);
sd = std(Zall, 1);
Zall = (Zall - mu) ./ sd;

npatt = size(Vpatt, 2);
pattResp = NaN(size(Zall, 1), npatt);
for k = 1:npatt
    %outer product of pattern k
    Pk = Vpatt(:,k) * Vpatt(:,k)';

    %removing diagonal elements of Pk to discard contributions from a
    %single neuron
    Pk(eye(size(Pk))>0) = 0;

    %computing the expression strength of pattern k as Zall(t,:) * Pk *
    %Zall(t,:)' for each t.
    Ztemp = Zall * Pk;
    pattResp(:,k) = sum(Ztemp .* Zall, 2);
end

%Detecting events of activation of each pattern as maxima when the 
%projected z-score value exceed a certain threshold.
%Initilizing the pattern response array
pattSpike = zeros(size(Zall, 1), npatt);
for k = 1:npatt
    binary = double(pattResp(:,k) > pattparams.strength_th);
    [~,loc] = findpeaks(pattResp(:,k) .* binary);
    pattSpike(loc,k) = 1;
end

%%
%Populate the output structure with results to be saved
pattparams.traintidx = traintidx;
Pat.pattparams = pattparams;

ncells_orig = size(Srep, 2);
Pat.weights = zeros(ncells_orig, size(Vpatt, 2));
Pat.cellAssemblies = false(ncells_orig, size(cellAss, 2));

Pat.weights(cellidx,:) = Vpatt;
Pat.cellAssemblies(cellidx,:) = cellAss;

Pat.Sparity = PSparity;
Pat.strength = pattResp;
Pat.activation = pattSpike;
end