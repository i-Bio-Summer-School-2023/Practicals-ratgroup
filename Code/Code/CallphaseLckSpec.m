
% Calculate the Spike locking to slow (0-20Hz) oscillations


function CallphaseLckSpec3(Spk,Lfp,Nav)

% for pyramidal cells

% index to pyramidal cells
pyrIdx = find(Spk.PyrCell==1);

% number of pyr cells
nPyr =  length(pyrIdx);


% Select only run section -- NEED the code
idxcond = Nav.Condition == 3;

% assuming run is niot disrupted, there is definately better way to do
% this but I dotn know all paramters stored
runTstart = Nav.sampleTimes(find(idxcond,1,'first'));
runTend = Nav.sampleTimes(find(idxcond,1,'last'));

LfpRun = Lfp.LfpHpc_raw;
LfpRun = LfpRun(Lfp.sampleTimes_raw>runTstart & Lfp.sampleTimes_raw<runTend );
LfpRun_timstamps = Lfp.sampleTimes_raw(Lfp.sampleTimes_raw>runTstart & Lfp.sampleTimes_raw<runTend );


% the frequencies used for spectrum
spec_Frq = logspace(log10(1),log10(20),96);
% spec_Frq = logspace(log10(1),log10(20),24);

% LFP ampling frequency
Lfp_Fs = 1./mean(diff(LfpRun_timstamps));

% calculatig the slow osc phase and power by morelet transform
[~, Phs] = MorlPow2(LfpRun, spec_Frq, Lfp_Fs, 7, 3); % wavelet power

spec_Phs = cell2mat(Phs');



for cellNumber=1:nPyr

    fignumber = cellNumber;

    % calculate and plot spike-LFP locking spectrum
    %-------------------------

    % find spikes of currrent cell
    spk_timestamps = Spk.spikeTimes(Spk.spikeID == cellNumber);

    % limited to run session
    spk_timestamps = spk_timestamps(spk_timestamps>runTstart & spk_timestamps<runTend );

    % caclulate
    [lckSpec] = phaseLckSpec(spk_timestamps, LfpRun_timstamps, spec_Phs, spec_Frq);

    % plot 
    plotphaseLckSpec(cellNumber,lckSpec,fignumber);
   

    % also draw spike-theta phas emodulation:    
    %-------------------------
    % calc spike phase by linear interpolation (NaN for times outside eeg timing)
    spikeThPhase = interp1(Lfp.sampleTimes_theta,Lfp.ThetaPhase,spk_timestamps','linear',NaN);
    
    %removing nan values
    spikeThPhase(isnan(spikeThPhase))=[];
    
    BinSize = 30;
    BarBin=15:BinSize:345;

    % calculate theta phase modulation paramters
    [phLck]= phaseLocking(spikeThPhase, BarBin);
    
    nHist = phLck.nHist;
    pRay    = phLck.pR;
    
    % plot theta-phase modulation
    histPlot(spikeThPhase,pRay,BarBin, BinSize,nHist,fignumber)
    %-------------------------


    input('enter to go next cell ..')    
end






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sigPow, sigPhs] = MorlPow2(sig, totF, Fs, NW, WaveDur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Power calcuations made in cell (instead of matrix), because cells do not
% have to be contiguous and so Out of Memory error wont happen. Mind the
% size if use matrices instead (600*4800*160*8/1024^3 = 3.4332GB)

% In this version, the phase and power are normalized seperately for
% precise acquiring both power and phase values

% str = strcat('computing wavelet power for f:',...
%     num2str(totF(1),'%2.2f, '),' ',num2str(totF(2),' %2.2f, '),'.., ',...
%     num2str(totF(end),' %2.2f Hz'));
% fprintf('%s \n',str);

% signal as a horizontal vector
sig = sig(:)';

Nfrq = length(totF);

sigPow =cell(1,Nfrq);
sigPhs =cell(1,Nfrq);

for idx=1:Nfrq
    f0 = totF(idx);
    St = NW/(2*pi*f0);
    wt  = -WaveDur*St:1/Fs:WaveDur*St;
    
    % Wavelet Normalization
    A = (St*sqrt(pi))^(-1/2);
    
    nMor = length(wt);
    nMorHalf = floor(nMor/2);
    nSig = length(sig);
    nConv = nMor+nSig-1;
    
    nConvPow2 = pow2(nextpow2(nConv));
    wMorl = A .* exp((-wt.^2)/(2*St^2)) .* exp(2*pi*f0*1i*wt);
    
    sigFFT  = fft(sig,nConvPow2);
    MorlFFT = fft(wMorl, nConvPow2);
    
    % acquiring the signal phase --------------------------------------
    sigConv = ifft(sigFFT .* MorlFFT);
    sigConv = sigConv(1:nConv);
    sigConv = 2*sigConv(nMorHalf:nSig+nMorHalf-1);
    
    sigPhs{idx} = angle(sigConv).*(180/pi);  % phase of signal in degree (-180 to 180)
    sigPhs{idx} = sigPhs{idx} + 180;         % Shift to set Trough as 0°. Degrees between (0 to +360)
    
    % acquiring the signal power --------------------------------------
    %Normalization of wavelet power over its amplitude
    MorlFFT = MorlFFT ./ max(MorlFFT);
    
    sigConv = ifft(sigFFT .* MorlFFT);
    sigConv = sigConv(1:nConv);
    sigConv = 2*sigConv(nMorHalf:nSig+nMorHalf-1);
    
    sigPow{idx} = abs(sigConv).^2;           % amplitude of signal
end







% function phaseLckSpecPerCell(cellNumber, spk_timestamps, Lfp_timstamps, spec_Frq, spec_Phs,fignumber)
% 
% % constants
% 
%  [lckSpec] = phaseLckSpec(spk_timestamps, Lfp_timstamps, spec_Phs, spec_Frq);
%  plotphaseLckSpec(cellNumber,lckSpec,fignumber);












% creating spike to osc-phase locking spectrum
function [lckSpec] = phaseLckSpec(spk_timestamps, Lfp_timstamps, spec_Phs, spec_Frq)

% phase bins
bin = 1:2:360;

% Theta (range) locking
%--------------------------------------
% phase input
Phs = spec_Phs;

% eeg (phase) time
% eeg_timstamps = (1:size(Phs,2))./specFs;

for ii=1:size(Phs,1)
    % calc spike phase by linear interpolation (NaN for times outside eeg timing)
    tsPh = interp1(Lfp_timstamps,Phs(ii,:),spk_timestamps','linear',NaN);
    
    %removing nan values
    tsPh(isnan(tsPh))=[];
    
    % calculate phase locking paramters
    [phLck]= phaseLocking(tsPh, bin);
    
    nHist = phLck.nHist;
    
    % smoothing the spectra
    nHist = repmat(nHist,1,3);
    nHist = kernelSmooth(nHist,length(bin),'guassian');
    nHist = nHist(length(bin)+1:2*length(bin));
    
    % creating output
    lckSpec.LO.nHist(ii,:) = nHist;
    %     lckSpec.LO.ppc(ii)     = phLck.ppc;
    lckSpec.LO.rL(ii)      = phLck.rL;
    lckSpec.LO.rAng(ii)    = phLck.rAng;
    lckSpec.LO.pR(ii)      = phLck.pR;
    % lckSpec.LO.D_rL(ii)    = phLck.D_rL;
end


% % % 
% % % % Gamma (range) locking
% % % %--------------------------------------
% % % 
% % % % phase input
% % % Phs = phsHO;
% % % 
% % % % eeg (phase) time
% % % tt = (1:size(Phs,2))./specFs;
% % % 
% % % for ii=1:size(Phs,1)
% % %     % calc spike phase by linear interpolation (NaN for times outside eeg timing)
% % %     tsPh = interp1(tt,Phs(ii,:),ts','linear',NaN);
% % %     
% % %     %removing nan values
% % %     tsPh(isnan(tsPh))=[];
% % %     
% % %     % calculate phase locking paramters
% % %     [phLck]= phaseLocking(p, tsPh, bin, 0);
% % %     
% % %     nHist = phLck.nHist;
% % %     
% % %     % smoothing the spectra
% % %     nHist = repmat(nHist,1,3);
% % %     nHist = kernelSmooth(nHist,length(bin),'guassian');
% % %     nHist = nHist(length(bin)+1:2*length(bin));
% % %     
% % %     % creating output
% % %     lckSpec.HO.nHist(ii,:) = nHist;
% % %     %     lckSpec.HO.ppc(ii)     = phLck.ppc;
% % %     lckSpec.HO.rL(ii)      = phLck.rL;
% % %     lckSpec.HO.rAng(ii)    = phLck.rAng;
% % %     lckSpec.HO.pR(ii)      = phLck.pR;
% % %     lckSpec.HO.D_rL(ii)    = phLck.D_rL;
% % % end

% creating output
%--------------------------------------
lckSpec.phaseBin = bin;
lckSpec.LO.frq   = spec_Frq;
% % lckSpec.HO.frq   = frqHO;



% measures statistics on spaikes phases (phase of a rythm during spikings)
function [phLck]= phaseLocking(spikePhase, bin)

% make sure spikePhase is coloumn-vector
spikePhase = spikePhase(:);

% if no spike
if isnan(spikePhase)
    phLck.nHist = NaN.*ones(length(bin));
    phLck.ppc   = NaN;
    phLck.rL    = NaN;
    phLck.rAng  = NaN;
    phLck.pR    = NaN;
    % phLck.D_rL  = NaN;
    return;
end

% Number of spikes
numSpk = length(spikePhase);

% Remove possible NaN
spikePhase(isnan(spikePhase)) = [];

% Firing rate. Note that if gamma cycle or theta cycle not selected all
% rates represent rateAll
% rateAll = numAll/lengthAll;

% Degree --> Radian spike phase
spikeRadPhase = circ_ang2rad(spikePhase);


% Resultant vector length
if numSpk==0; rL=NaN; else rL = circ_r(spikeRadPhase); end

% Resultant vector phase [deg]
if numSpk==0; rAng=NaN; else rAng = circ_rad2ang(circ_mean(spikeRadPhase)); end

% Rayleigh test (circ_rtest)
if numSpk==0; pR=NaN; else pR = circ_rtest(spikeRadPhase); end

% spike count $ normalized
nHist = hist(spikePhase,bin);
nHist = nHist./numSpk;


phLck.nHist = nHist;
phLck.rL    = rL;
phLck.rAng  = rAng;
phLck.pR    = pR;




function alpha = circ_ang2rad(alpha)
alpha = alpha * pi /180;

function alpha = circ_rad2ang(alpha)
alpha = alpha / pi *180;

% smoothing a signal with a kernel with specific length
function [outSig]=kernelSmooth(inSig,windowLen,type)

switch type
    case 'uniform'
        kernl = ones(1, windowLen)./windowLen;
    case 'guassian'
        kernl = gausswin(windowLen)./sum(gausswin(windowLen));
end

outSig = conv(inSig, kernl, 'same');







function r = circ_r(alpha, w, d, dim)
% r = circ_r(alpha, w, d)
%   Computes mean resultant vector length for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%     [dim  compute along this dimension, default is 1]
%
%     If dim argument is specified, all other optional arguments can be
%     left empty: circ_r(alpha, [], [], dim)
%
%   Output:
%     r		mean resultant length
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%   Topics in circular statistics, S.R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if nargin < 4
    dim = 1;
end

if nargin < 2 || isempty(w)
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) ~= size(alpha,2) || size(w,1) ~= size(alpha,1)
        error('Input dimensions do not match');
    end
end

if nargin < 3 || isempty(d)
    % per default do not apply correct for binned data
    d = 0;
end

% compute weighted sum of cos and sin of angles
r = sum(w.*exp(1i*alpha),dim);

% obtain length
r = abs(r)./sum(w,dim);

% for data with known spacing, apply correction factor to correct for bias
% in the estimation of r (see Zar, p. 601, equ. 26.16)
if d ~= 0
    c = d/2/sin(d/2);
    r = c*r;
end




function [mu ul ll] = circ_mean(alpha, w)
%
% mu = circ_mean(alpha, w)
%   Computes the mean direction for circular data.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		weightings in case of binned angle data]
%
%   Output:
%     mu		mean direction
%     ul    upper 95% confidence limit
%     ll    lower 95% confidence limit
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

% check vector size
if size(alpha,2) > size(alpha,1)
    alpha = alpha';
end

if nargin<2
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) > size(w,1)
        w = w';
    end
end

% compute weighted sum of cos and sin of angles
r = w'*exp(1i*alpha);

% obtain mean by
mu = angle(r);

% confidence limits if desired
if nargout > 1
    t = circ_confmean(alpha,0.05,w);
    ul = mu + t;
    ll = mu - t;
end



function [pval z] = circ_rtest(alpha, w, d)
%
% [pval, z] = circ_rtest(alpha,w)
%   Computes Rayleigh test for non-uniformity of circular data.
%   H0: the population is uniformly distributed around the circle
%   HA: the populatoin is not distributed uniformly around the circle
%   Assumption: the distribution has maximally one mode and the data is
%   sampled from a von Mises distribution!
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%
%   Output:
%     pval  p-value of Rayleigh's test
%     z     value of the z-statistic
%
% PHB 7/6/2008
%
% References:
%   Statistical analysis of circular data, N. I. Fisher
%   Topics in circular statistics, S. R. Jammalamadaka et al.
%   Biostatistical Analysis, J. H. Zar
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html

if size(alpha,2) > size(alpha,1)
    alpha = alpha';
end

if nargin < 2
    r =  circ_r(alpha);
    n = length(alpha);
else
    if length(alpha)~=length(w)
        error('Input dimensions do not match.')
    end
    if nargin < 3
        d = 0;
    end
    r =  circ_r(alpha,w(:),d);
    n = sum(w);
end

% compute Rayleigh's R (equ. 27.1)
R = n*r;

% compute Rayleigh's z (equ. 27.2)
z = R^2 / n;

% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-R^2))-(1+2*n));







% This function plots spike-LFP phase locking spectra
function [] = plotphaseLckSpec(cellNumber,lckSpec,fignumber)

% extracting inputs
x = lckSpec.phaseBin;
nHistLO = lckSpec.LO.nHist;
frqLO   = lckSpec.LO.frq;

% nHistHO = lckSpec.HO.nHist;
% frqHO   = lckSpec.HO.frq;

% ppcLO = lckSpec.LO.ppc;
% ppcHO = lckSpec.HO.ppc;

% Normalize on phase bins
nHistLO2 = zscore(nHistLO')';
% nHistHO2 = zscore(nHistHO')';

% Plots
%--------------------------------------
nHist = {nHistLO;  nHistLO2;};

% initializing figure;
figure(fignumber); clf;
% set(gcf,'units','normalized','Position',figpos)

for ii=1:length(nHist)

    % Theta-----
    subplot(2,2,ii); hold on;
    surf([x x+360], frqLO, [nHist{ii} nHist{ii}],'EdgeColor','none');
    colormap(jet);
    colorbar('off');
    % colorbar('eastoutside');
    view(2);
    % mm = max(max(nHistTh));
    % plot3(0:1:720,6*ones(length(0:1:720)),mm.*ones(length(0:1:720)),...
    %     'Color','black','LineStyle','--','LineWidth',1);
    % plot3(0:1:720,10*ones(length(0:1:720)),mm.*ones(length(0:1:720)),...
    %     'Color','black','LineStyle','--','LineWidth',1);
    set(gca,'FontSize',10,'xlim',[0 720],'ylim',[1 20],...
        'xtick',[0 180 360 540 720],...
        'YScale', 'log','ytick',[1 6 10 20],'box','off')

    strTitle = strcat('cell# ',num2str(cellNumber,'%g'));
    if ii==1
        title(strTitle,'FontSize',10)
        xlabel('LFP phase (°)')
        ylabel('LFP frq (Hz)')
    else
        title('normalized','FontSize',10)
    end
    % colorbar
    
    % Gamma-----
    % % subplot(2,2,ii+1); hold on;
    % % surf([x x+360], frqHO, [nHist{ii+1} nHist{ii+1}],'EdgeColor','none');
    % % colormap(jet);
    % % colorbar('off');
    % % % colorbar('eastoutside');
    % % view(2);
    % % % mm = max(max(nHistGm));
    % % % plot3(0:1:720,25*ones(length(0:1:720)),mm.*ones(length(0:1:720)),...
    % % %     'Color','black','LineStyle','--','LineWidth',1);
    % % % plot3(0:1:720,50*ones(length(0:1:720)),mm.*ones(length(0:1:720)),...
    % % %     'Color','black','LineStyle','--','LineWidth',1);
    % % % plot3(0:1:720,100*ones(length(0:1:720)),mm.*ones(length(0:1:720)),...
    % % %     'Color','black','LineStyle','--','LineWidth',1);
    % % set(gca,'FontSize',10,'xlim',[0 720],'ylim',[20 140],...
    % %     'xtick',[0 180 360 540 720],...
    % %     'YScale', 'log','ytick',[25 50 100 140],'box','off')
    % % if ii==1
    % %     title('Gamma phase locking','FontSize',12)
    % % else
    % %     title('z-score normalized','FontSize',12)
    % % end
end

% keyboard;


% plotting the  resultant length of vector length for Spike-phase count

LockingResLen = lckSpec.LO.rL;

subplot(2,2,3); hold on;
plot(frqLO,LockingResLen)
xlabel('LFP frequency (Hz)');
ylabel('strength of locking (1)')

%-----------------------------------------------
% save image
% saveImage(gcf,p,session,tetrode,unit,figext)







% plot the spike-theta phase distrbution
function histPlot(spikePhase,pR,BarBins,BinSize,nHist,figNumber)

figure(figNumber);
subplot(2,2,4)
if isnan(spikePhase)
    return;
end

% Degree --> Radian spike phase
spikeRadPhase = circ_ang2rad(spikePhase);

% Histgram-----------------------------------------------------------------
MinH = 0.00;
MaxH = 0.2;
xDetail = 0:360;
lblesFsize = 10;

% Theta -------------------------------------------------------------------
if pR<0.05
    bar([BarBins BarBins+360],[nHist nHist],1,'r')
    title ('p(Rayleigh)<0.05')
else
    bar([BarBins BarBins+360],[nHist nHist],1,'b')
end
% Curve fitting
if ~isempty(spikePhase)
    % von Mises fitting of original data
    [mu kappa] = circ_vmpar(spikeRadPhase);
    % get von Mises distribution
    vm = circ_vmpdf(circ_ang2rad(xDetail),mu,kappa);
    vm = vm/sum(vm)*length(xDetail)/(360/BinSize);
    vm=vm(:);
    %plot
    hold on
    plot([xDetail xDetail+360],[vm' vm'],'k-')
    %     title('von Mises fitting','FontSize',lblesFsize)

end


% Setting
set(gca,'XTick',[0 180 360 540 720]); axis([0 720 MinH MaxH]);
xlabel('theta phase (°)','FontSize',lblesFsize);
ylabel('spike count (norm.)','FontSize',lblesFsize);
figSetting(lblesFsize-2);




function [thetahat kappa] = circ_vmpar(alpha,w,d)

% r = circ_vmpar(alpha, w, d)
%   Estimate the parameters of a von Mises distribution.
%
%   Input:
%     alpha	sample of angles in radians
%     [w		number of incidences in case of binned angle data]
%     [d    spacing of bin centers for binned data, if supplied
%           correction factor is used to correct for bias in
%           estimation of r, in radians (!)]
%
%   Output:
%     thetahat		preferred direction
%     kappa       concentration parameter
%
% PHB 3/23/2009
%
% References:
%   Statistical analysis of circular data, N.I. Fisher
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de

alpha = alpha(:);
if nargin < 2
    w = ones(size(alpha));
end
if nargin < 3
    d = 0;
end

r = circ_r(alpha,w,d);
kappa = circ_kappa(r);

thetahat = circ_mean(alpha,w);




function kappa = circ_kappa(alpha,w)
%
% kappa = circ_kappa(alpha,[w])
%   Computes an approximation to the ML estimate of the concentration
%   parameter kappa of the von Mises distribution.
%
%   Input:
%     alpha   angles in radians OR alpha is length resultant
%     [w      number of incidences in case of binned angle data]
%
%   Output:
%     kappa   estimated value of kappa
%
%   References:
%     Statistical analysis of circular data, Fisher, equation p. 88
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de - www.kyb.mpg.de/~berens/circStat.html


alpha = alpha(:);

if nargin<2
    % if no specific weighting has been specified
    % assume no binning has taken place
    w = ones(size(alpha));
else
    if size(w,2) > size(w,1)
        w = w';
    end
end

N = length(alpha);

if N>1
    R = circ_r(alpha,w);
else
    R = alpha;
end

if R < 0.53
    kappa = 2*R + R^3 + 5*R^5/6;
elseif R>=0.53 && R<0.85
    kappa = -.4 + 1.39*R + 0.43/(1-R);
else
    kappa = 1/(R^3 - 4*R^2 + 3*R);
end

if N<15 && N>1
    if kappa < 2
        kappa = max(kappa-2*(N*kappa)^-1,0);
    else
        kappa = (N-1)^3*kappa/(N^3+N);
    end
end



function [p alpha] = circ_vmpdf(alpha, thetahat, kappa)

% [p alpha] = circ_vmpdf(alpha, w, p)
%   Computes the circular von Mises pdf with preferred direction thetahat
%   and concentration kappa at each of the angles in alpha
%
%   The vmpdf is given by f(phi) =
%   (1/(2pi*I0(kappa))*exp(kappa*cos(phi-thetahat)
%
%   Input:
%     alpha     angles to evaluate pdf at, if empty alphas are chosen to
%               100 uniformly spaced points around the circle
%     [thetahat preferred direction, default is 0]
%     [kappa    concentration parameter, default is 1]
%
%   Output:
%     p         von Mises pdf evaluated at alpha
%     alpha     angles at which pdf was evaluated
%
%
%   References:
%     Statistical analysis of circular data, Fisher
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens and Marc J. Velasco, 2009
% velasco@ccs.fau.edu

% if no angles are supplied, 100 evenly spaced points around the circle are
% chosen
if nargin < 1 || isempty(alpha)
    alpha = linspace(0, 2*pi, 101)';
    alpha = alpha(1:end-1);
end
if nargin < 3
    kappa = 1;
end
if nargin < 2
    thetahat = 0;
end

alpha = alpha(:);

% evaluate pdf
C = 1/(2*pi*besseli(0,kappa));
p = C * exp(kappa*cos(alpha-thetahat));


% Seting figure parameters
function figSetting(txtSize)

set(gca,'Box','Off')
set(gca, 'TickDir', 'Out');
set(gca, 'TickLength', [0.02 0.02]);
set(gca,'FontSize',txtSize);
set(gca,'XTick',[0 180 360 540 720])