function pattparams = SetPattParams(Nav,Srep)
% SetPattParams - Define parameters for detecting cell assembly patterns.
%
%   pattparams = SetPattParams(Nav, Srep) defines a set of parameters for
%   detecting cell assembly patterns from spike response data.
%
% INPUTS:
% - Nav: A structure containing at least a field called 'sampleTimes' with
%   the sample times of the data and some additional fields with the
%   explanatory variables
% - Srep: Array of responses (ntimes x ncells) from which patterns will be
%   identified.
%
% OUTPUT:
% - pattparams: A structure containing pattern analysis parameters with the
%   following fields:
%   - subset: A structure defining conditions over the fields of Nav for
%     which patterns will be detected. Fields of subset correspond to the
%     names of the fields in Nav, and the subset specifies values and
%     operations to apply. For example:
%     pattparams.subset.Condition = [1 3];
%     pattparams.subset.Condition_op = 'ismember';
%     pattparams.subset.Spd = 5;
%     pattparams.subset.Spd_op = '>=';
%   - cellidx: A logical array indicating a subset of cells for pattern detection.
%   - nspk_th: Minimal number of spikes over the train set to consider a cell.
%   - Marcenko: If true, select principal components (PCs) based on the Marcenko-Pastur law.
%   - nShuffle: Number of shuffle controls to perform for eigenvalue distribution
%     if pattparams.Marcenko is false.
%   - variablenames: Cell array of variable names in Nav used for shuffling.
%   - binedges: Cell array of bin edges for discretizing variables.
%   - NoiseCov: If true, remove average covariance from covariance matrix.
%   - pvalshf_th: P-value threshold for selecting PCs using shuffling.
%   - strength_th: Pattern activation threshold to convert strength into "spikes".
%   - sampleRate: Sampling rate of the data (in Hz).
%   - timewin: Size of the spike count window in seconds.
%
% USAGE:
%    Nav = LoaddataNav(loadparams);
%    Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
%    Srep = Spk.spikeTrain;
%    pattparams = SetPattParams(Nav, Srep);
%
% Written by J Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"
%%
%Conditions over the fields of Nav for which patterns will be detected
%pattparams.subset should be a structure where fields have names of the 
%fields of Nav to which the condition should apply to.
pattparams.subset = [];

%For instance, for the example data set, we can define the following fields
pattparams.subset.Condition = [1 3 5];
pattparams.subset.Condition_op = 'ismember';

pattparams.subset.XDir = [-1 1];
pattparams.subset.XDir_op = 'ismember';

pattparams.subset.laptype = [-1 0 1];
pattparams.subset.laptype_op = 'ismember';

pattparams.subset.Spd =  2.5;
pattparams.subset.Spd_op = '>=';

pattparams.subset.Xpos =  0;
pattparams.subset.Xpos_op = '>=';

pattparams.subset.Xpos =  100;
pattparams.subset.Xpos_op = '<=';

%Subset of cells used for pattern detection.
pattparams.cellidx = true(1, size(Srep, 2));

%Minimal number of spikes over the train set to consider a cell for 
%pattern detection
pattparams.nspk_th = 0;

%If true, the PC will be selected according to the Marcenko-Pastur law
pattparams.Marcenko = true;

%Number of shuffle controls to perform for randomization if Marcenko-Pastur
%law is not used to select PCs
pattparams.nShuffle = 100;

%Names of varaibles in Nav that will be used to build the shuffle
%distribution of eigenvalues. Time points will be shuffled wihtin bins of
%these variables. If empty, the shuffling procedure will simply circularly 
%shift the spike counts by a random amount > 1 second.
pattparams.variablenames{1} = 'Xpos';
pattparams.variablenames{2} = 'Spd';
pattparams.variablenames{3} = 'XDir';

%Bin edges to discretize variables indicated in pattparams.variablenames.
pattparams.binedges{1} = 0 : 2: 100;
pattparams.binedges{2} = [0 : 5 : 50 inf];
pattparams.binedges{3} = [-2 0 2];

%If true, the covariance average across all shuffle controls, which provide
%an estimate of the signal covariance, will be removed from the overall
%covariance before proceeding to the PCA and ICA.
pattparams.NoiseCov = false;

%P-value to use as a threshold when selecting the PCs using the shuffling
%approach
pattparams.pvalshf_th = 0.05;

%Pattern activation threshold to convert activation strength into
%activation "spikes".
pattparams.strength_th = 5;

%Sampling rate of the data
pattparams.sampleRate = 1 / mean(diff(Nav.sampleTimes), 'omitnan');

%Size of the spike count window in seconds
pattparams.timewin = .02;
end