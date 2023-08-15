function PhsMaps = PhsMapsAnalysis(Nav, stcell, Lrep, Phsmapsparams, sampleTimes_Lrep)

%%
%If no Y variable are indicated in mapsparams, we'll just compute a 1D map
if ~isempty(Phsmapsparams.Xvariablename)
    X = Nav.(Phsmapsparams.Xvariablename);
else
    X = ones(size(Nav.sampleTimes));
    Phsmapsparams.Xbinedges = 1;
    Phsmapsparams.XsmthNbins = 0;
end

%Selecting time indices over which to perform the phase coupling analysis,
%according to parameters defined in Phsmapsparams.subset
tidx = true(size(X));
pnames = fieldnames(Phsmapsparams.subset);
fnames = fieldnames(Nav);
for i = 1:numel(pnames)
    if ismember(pnames{i},fnames)
        fn = str2func(Phsmapsparams.subset.([pnames{i} '_op']));
        tidx = tidx & fn(Nav.(pnames{i}), Phsmapsparams.subset.(pnames{i}));
    elseif ~strcmp(pnames{i}(end-2:end),'_op')
        warning('some fields of Phsmapsparams.subset are not matching fields of Nav')
    end
end

%Selecting time indices where X is in the range of bins
tidx = tidx &...
       X >= Phsmapsparams.Xbinedges(1) & X <= Phsmapsparams.Xbinedges(end) &...
       ~isnan(X);

%number of cells selected for phase coupling analysis
ncells = size(stcell, 2);

for icell = 1:ncells
    tidxspk = interp1(Nav.sampleTimes, tidx, stcell{icell}, 'nearest');
    stcell{icell}(~tidxspk) = [];
end

%Number of spikes for each cell
nspk = cellfun(@length, stcell);

%Selecting cells for which we'll run the phase coupling analysis
cellidx = find(Phsmapsparams.cellidx & nspk > Phsmapsparams.nspk_th);

%Subsetting the cell array of spike times
st = stcell(cellidx);

%number of cells selected for phase coupling analysis
ncells = size(st, 2);

%number of phase bins
nPhsbins = numel(Phsmapsparams.Phsbinedges) - 1;

%Subsetting X and Nav.sampleTimes
X = X(tidx);
sampleTimes_Nav = Nav.sampleTimes(tidx);

%%
%Filtering Lrep data in the frequency band indicated in 
%Filtering in the theta band
LsampleRate = 1 / mean(sampleTimes_Lrep, 'omitnan');
Lfilt = filterellip(Lrep, Phsmapsparams.freqrange(1), Phsmapsparams.freqrange(2), LsampleRate);

%Computing Hilbert transform of the filtered signal to get the 
%instantaneous phase of the pass-band signal.
hilbTrans = hilbert(Lfilt);
phs = angle(hilbTrans);

%subsetting phs and sampleTimes_Lrep
tidxLrep = interp1(Nav.sampleTimes, tidx, sampleTimes_Lrep, 'nearest');
phs = phs(tidxLrep);
sampleTimes_Lrep = sampleTimes_Lrep(tidxLrep);

%discretizing phases according to Phsmapsparams.Phsbinedges
phs_discrete = discretize(phs, Phsmapsparams.Phsbinedges);

%%
%Statistics based on distribution of spike theta phase.

%Initializing arrays.
scmap = NaN(ncells, nPhsbins);%spike count array
rL = NaN(ncells,1);%Resultant length array
Rpval = NaN(ncells,1);%Rayleigh test p-value
phsmean = NaN(ncells,1);%Rayleigh test p-value

occmap = ComputeMap(phs_discrete, [], ones(size(phs_discrete)), nPhsbins);

%Looping across cells to build phase-locking statistics
for icell = 1:ncells
    %Getting the phase at spike times by interpolation
    sph = interp1(sampleTimes_Lrep, phs, st{icell}, 'nearest');

    %Building histogram of spike phases.
    sph_discrete = discretize(sph_discrete, Phsmapsparams.Phsbinedges);
    scmap(icell,:) = ComputeMap(sphs_discrete, [], ones(size(sphs_discrete)), nPhsbins);

    %Resultant vector length.
    rL(icell) = circ_r(sph / 180 * pi);
    phsmean(icell) = circ_mean(sph / 180 * pi) / pi * 180;
    Rpval(icell) = circ_rtest(sph / 180 * pi);

end

%%
% Another way to compute the modulation of spiking by theta phase is by 
% computing a map (i.e. spike count map / occupancy map).

% %Retrieving the default parameters to compute maps.
mapsparams = SetMapsParams(Nav,stcell);

%Replacing values for time indices selection.
mapsparams.subset = phsmapsparams.subset;
mapsparams.Yvariablename = [];
mapsparams.Xvariablename = 'Phs';
mapsparams.XsmthNbins = phsmapsparams.PhssmthNbins;
mapsparams.Xbinedges = phsmapsparams.Phsbinedges;

%Creating a field in Nav corresponding to the phase.
Nav.ThetaPhase = thetaphs;

%STOPPED HERE:need to change the behavior of MapsAnalysis to work with
%spike times


%Now we can run MapsAnalysis with the new Nav and phsmapsparams.
maps = MapsAnalyses(Nav, Spk.spikeTrain, phsmapsparams);

%%
%First plotting position versus theta phase of spikes (after splitting
%according to the direction of travel...)
% figure;
% for icell = 1:ncells
% scatter(Nav.Xpos(tidx & Spk.spikeTrain(:,cellidx(icell))>0),Lfp.ThetaPhase(tidx & Spk.spikeTrain(:,cellidx(icell))>0),'.');
% pause
% end

%%
%Construct 2D maps of mean firing rate as a function of position and theta
%phase.

%number of position and theta phase bins
nPhsbins = numel(Phsmapsparams.Phsbinedges) - 1;
nXbins = numel(Phsmapsparams.Xbinedges) - 1;

%Discretizing theta phases
phs = Lfp.ThetaPhase(tidx);
phs_discrete = discretize(phs, Phsmapsparams.Phsbinedges);

%Discretizing positions
Xpos = Nav.Xpos(tidx);
Xpos_discrete = discretize(Xpos, Phsmapsparams.Xbinedges);

%Computing occupancy map
flat = 1/Phsmapsparams.sampleRate * ones(size(Xpos_discrete));
occmap = ComputeMap(phs_discrete, Xpos_discrete, flat, nPhsbins, nXbins);

%Removing occupancy for position x theta phase bins below the occupancy 
%threshold
occmap(occmap <= Phsmapsparams.occ_th) = NaN;

%Smoothing the occupancy map with a 2D gaussian window.
occmap = repmat(occmap, [1 3]);
occmap = GaussianSmooth(occmap, [Phsmapsparams.XsmthNbins Phsmapsparams.PhssmthNbins]);
occmap = occmap(:,(nPhsbins+1):2*nPhsbins);

%Computing and smoothing the spike count map for each cell
scmap = NaN(ncells, nXbins, nPhsbins);
for icell = 1:ncells
    scmaptemp = ComputeMap(phs_discrete, Xpos_discrete, spikeTrain(:,icell), nPhsbins, nXbins);
    scmaptemp(isnan(occmap)) = NaN;
    scmaptemp = repmat(scmaptemp, [1 3]);
    scmaptemp = GaussianSmooth(scmaptemp, [Phsmapsparams.XsmthNbins Phsmapsparams.PhssmthNbins]);
    scmap(icell,:,:) = scmaptemp(:, (nPhsbins+1):2*nPhsbins);
end

%Calculating the place field x theta phase maps by dividing scmap and 
%occmap
occmap = permute(occmap, [3 1 2]);%permuting dimension for convenience
mapXTheta = scmap ./ occmap;

%%

%Estimating the average decoding error as a function of the phase of the 
% theta oscillation

%Running the decoder (direction x position) with the default parameters
decparams = DefineDecParams(Nav, Spk);
Dec = DecodingAnalysis2D(Nav, Spk.spikeTrain, decparams);

%Calculating the decoding error 
XErrMax = (Dec.XDecMax - Dec.X) .* sign(Nav.XDir);

%Selecting time indices avoiding edges
minmax = [25 75];
tidxdec = tidx & Nav.Xpos > minmax(1) &  Nav.Xpos < minmax(2);
phs = Lfp.ThetaPhase(tidxdec);

%Discretizing theta phases
phs_discrete = discretize(phs, Phsmapsparams.Phsbinedges);

%Computing and smoothing cicularly the sum of decoding errors across theta
%phases
summap = ComputeMap(phs_discrete, [], XErrMax(tidx & Nav.Xpos > minmax(1) &  Nav.Xpos < minmax(2)), nPhsbins, []);
summap = GaussianSmooth(repmat(summap,[1 3]), [0 1]);
summap = summap((nPhsbins+1):2*nPhsbins);

%Computing the occupancy map
occmap = ComputeMap(phs_discrete, [], ones(size(phs_discrete)), nPhsbins, []);
occmap = GaussianSmooth(repmat(occmap,[1 3]), [0 1]);
occmap = occmap((nPhsbins+1):2*nPhsbins);

%Calculating the average decoding error across theta phases
ThetaXDec = summap ./ occmap;

%%
ncells_orig = size(Srep, 2);

PhsMaps.rL = NaN(ncells_orig,1);
PhsMaps.phsmean = NaN(ncells_orig,1);
PhsMaps.Rpval = NaN(ncells_orig,1);

PhsMaps.rL(cellidx) = rL;
PhsMaps.phsmean(cellidx) = phsmean;
PhsMaps.Rpval(cellidx) = Rpval;

PhsMaps.Phsmaps = maps;
PhsMaps.mapXTheta = mapXTheta;
PhsMaps.ThetaXDec = ThetaXDec;

%    %%
% %Same thing as above but after centering position on the position of max
% %firing rate and normalizing positions by the field width.
% 
% %number of position and theta phase bins
% nPhsbins = numel(Phsmapsparams.Phsbinedges) - 1;
% 
% %Discretizing theta phases
% phs = Lfp.ThetaPhase(tidx);
% phs_discrete = discretize(phs, Phsmapsparams.Phsbinedges);
% 
% %position of max firing rate
% mapX = Maps.mapX(cellidx,:);
% [~, imaxPos] = max(mapX, [], 2);
% maxPos = Maps.bincenters(imaxPos);
% 
% %Identifying limits of subfield edges.
% ifieldstart = NaN(1, ncells);
% ifieldend = NaN(1, ncells);
% field_th = 0.10;%amplitude threshold to identify field limits
% for icell = 1:ncells
%     ma = max(mapX(icell,:));
%     mi = min(mapX(icell,:));
%     startidx = find(mapX(icell,1:imaxPos(icell)) <= mi + field_th * (ma - mi), 1, 'last');
%     if isempty(startidx)
%         startidx = 1;
%     end
%     ifieldstart(icell) = startidx;
% 
%     endidx = imaxPos(icell) + find(mapX(icell,(imaxPos(icell)+1):end) <= mi + field_th * (ma - mi), 1, 'first');
%     if isempty(endidx)
%         endidx = size(mapX, 2);
%     end
%     ifieldend(icell) = endidx;
% end
% 
% fieldstart = Maps.bincenters(ifieldstart);
% fieldend = Maps.bincenters(ifieldend);
% 
% %Computing linear-circular coefficients and position x theta phase maps for
% %each cell after centering and nomralizing the position of the animal to
% %within field limits.
% Xbinedges_norm = -1:0.1:1;
% nXbins_norm = numel(Xbinedges_norm) - 1;
% XsmthNbins_norm = 2;
% mapXTheta2 = NaN(ncells, nXbins_norm, nPhsbins);
% ThetaRho = NaN(ncells, 1);
% ThetaRho_pval = NaN(ncells, 1);
% ThetaSlope = NaN(ncells, 1);
% ThetaPhi0 = NaN(ncells, 1);
% ThetaNspk = NaN(ncells, 1);
% for icell = 1:ncells
%     %Discretizing positions
%     if numel(unique(Maps.mapsparams.dir)) > 1
%         error('phase precession on linear track should be estimated from traversals in a single direction')
%     end
%     Xpos = (Nav.Xpos(tidx) - maxPos(icell)) * sign(Maps.mapsparams.dir);
%     Xpos(Xpos < 0) = Xpos(Xpos < 0) / abs(maxPos(icell) - fieldstart(icell));
%     Xpos(Xpos > 0) = Xpos(Xpos > 0) / abs(fieldend(icell) - maxPos(icell));
%     Xpos(abs(Xpos) > 1) = NaN;
%     Xpos_discrete = discretize(Xpos, Xbinedges_norm);
% 
%     %Estimating the linear-circular coefficient and its significance
%     spkidx = spikeTrain(:,icell) > 0;
%     [ThetaRho(icell) ,ThetaRho_pval(icell), ThetaSlope(icell) ,ThetaPhi0(icell)] = ...
%         circlin_regression(Xpos(spkidx), phs(spkidx) /180 *pi);
%     ThetaNspk(icell) = sum(spkidx);
% 
%     %Computing occupancy map
%     flat = 1/Phsmapsparams.sampleRate * ones(size(Xpos_discrete));
%     occmap = ComputeMap(phs_discrete, Xpos_discrete, flat, nPhsbins, nXbins_norm);
% 
%     %Removing occupancy for position x theta phase bins below the occupancy
%     %threshold
%     occmap(occmap <= Phsmapsparams.occ_th) = NaN;
% 
%     %Smoothing the occupancy map with a 2D gaussian window.
%     occmap = repmat(occmap, [1 3]);
%     occmap = GaussianSmooth(occmap, [XsmthNbins_norm Phsmapsparams.PhssmthNbins]);
%     occmap = occmap(:,(nPhsbins+1):2*nPhsbins);
% 
%     %Computing and smoothing the spike count map for each cell
% 
%     scmap = ComputeMap(phs_discrete, Xpos_discrete, spikeTrain(:,icell), nPhsbins, nXbins_norm);
%     scmap(isnan(occmap)) = NaN;
%     scmap = repmat(scmap, [1 3]);
%     scmap = GaussianSmooth(scmap, [Phsmapsparams.XsmthNbins Phsmapsparams.PhssmthNbins]);
%     scmap = scmap(:, (nPhsbins+1):2*nPhsbins);
% 
%     mapXTheta2(icell,:,:) = scmap ./ occmap;
% end
% 
% 
% %%
% %Decoding positions around ripple times
% 
% %Detecting ripple peaks
% ripFs = 1 / mean(diff(Lfp.sampleTimes));
% riptimes = find(Lfp.ripplepeak == 1);
% 
% 
% %Extracting snippets of Decoded positions around ripple peaks.
% idxwin = -round(0.1 * ripFs):round(0.1 * ripFs);
% [~, ~, lrip] = ComputeTriggeredAverage(Dec.XDecMax, riptimes, idxwin);

end