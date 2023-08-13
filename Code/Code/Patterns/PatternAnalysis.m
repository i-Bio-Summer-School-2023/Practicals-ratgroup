function Pat = PatternAnalysis(Nav, Srep, pattparams)
% PatternAnalysis - Identify cell assemblies in Srep, based on independent component analysis.
%
%
% Pat = PatternAnalysis(Nav, Srep, pattparams)
%
% Inputs:
%   Nav: Structure containing navigation data (timestamps, positions, speeds, etc.).
%   Srep: Spike train data for each neuron (timepoints x neurons).
%   pattparams: Structure with pattern analysis parameters.
%
% Outputs:
%   Pat: Structure containing pattern analysis results wtih the following
%   fields:
%   - weights: Weight matrix representing the contribution of each neuron to each assembly pattern.
%   - cellAssemblies: Logical matrix indicating which neurons belong to each detected cell assemblies.
%   - Sparity: Sparsity measure of each assembly pattern (values between 0 and 1).
%   - strength: Expression strength of each assembly pattern over time.
%   - activation: events of assembly pattern activation.
%
% See Also:
%   fastICA, GaussianSmooth, findpeaks
%
% Method based on the description in the supplementary information of
% van de Ven et al., 2016 in Neuron.
%
%
% Usage:
%   Pat = PatternAnalysis(Nav, Srep, pattparams)
%
% Written by J. Fournier in 08/2023 for the iBio Summer school.

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

%Time indices over which patterns will be detected
traintidx = ismember(Nav.Condition, pattparams.condition) &...
            ismember(Nav.XDir, pattparams.dir) &...
            Nav.Spd > pattparams.spdthreshold &...
            ~isnan(Nav.Xpos);

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
%matrix)
[P,D] = eig(Corrmat);
eigval = diag(D);

%Selecting only principal components with an eigenvalue higher than 
%expected (see van de Ven et al., 2016).
eigval_th = (1 + ncells / ntimepts)^2;
Psign = P(:,eigval > eigval_th);

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