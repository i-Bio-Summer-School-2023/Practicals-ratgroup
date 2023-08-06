%This is the order of the functions to go through 
%We still need to add some functions to compare different conditions and some
%functions to plot the results

%Day1
loadparams = DefineLoadParams;
Nav = LoaddataNav(loadparams);
Spk = LoaddataSpk(loadparams, Nav.sampleTimes);
Lfp = LoaddataLfp(loadparams, Nav.sampleTimes);

%Day2
mapsparams = DefineMapsParams(Nav,Spk);%define parameters
Maps = MapsAnalyses(Nav, Spk.spikeTrain, mapsparams);%estimate place fields
Maps = MapsAnalyses2D(Nav, Spk.spikeTrain, mapsparams);%estimate place fields
glmsparams = DefineGLMsParams(Nav,Spk);%define parameters for GLMs
GLMs = GLMAnalyses(Nav, Spk.spikeTrain, glmsparams);%estimate GLMs
GLMs = GLMAnalysesCV(Nav, Spk.spikeTrain, glmsparams);%estimate GLMs with cross-validation

%Day3
decparams = DefineDecParams(Nav,Spk);%define parameters for decoding
Dec = DecodingAnalysis(Nav, Spk.spikeTrain, decparams);%compute Bayesian decoding, 1D
DecXLR = DecodingAnalysis2D(Nav, Spk.spikeTrain, decparams);%compute Bayesian decoding, 2D
pattparams = DefinePattParams(Nav,Spk);%define parameters for pattern detection
Pat = PatternAnalysis(Nav, Spk.spikeTrain, pattparams);%Estimate cell assembly patterns
%Estimate Maps based on the vectors of activation of cell assemblies.
mapspattparams = DefineMapsParams(Nav,Spk);
mapspattparams.cellidx = true(1, size(Pat.spikeTrain, 2)); 
PatMaps = MapsAnalyses(Nav, Pat.activation, mapspattparams);

%optional: same for GLM analysis
glmspattparams = DefineGLMsParams(Nav,Spk);
glmspattparams.cellidx = true(1, size(Pat.activation, 2)); 
PatGLMs = GLMAnalyses(Nav, Pat.activation, glmspattparams);
%optional: same for decoding analysis
decpattparams = DefineDecParams(Nav,Spk);
decpattparams.cellidx = true(1, size(Pat.activation, 2)); 
PatDec = DecodingAnalysis(Nav, Pat.activation, decpattparams);
PatDecXLR = DecodingAnalysis2(Nav, Pat.activation, decpattparams);

%Day4
crossparams = DefineCrossParams(Nav,Spk);%define parameters for correlation analysis
Cross = CrossCorrelationAnalyses(Nav, Spk.spikeTrain, crossparams);%Compute pair-wise correlations
TFparams = DefinelfpTFparams(Nav, Spk, Lfp);%Define parameters for time frequency analyses
TFlfp = TimeFreqAnalysis(Nav, Spk, Lfp, TFparams);%wavelet transform and coherence
thetaparams = DefineThetaParams(Nav,Spk);%Parameters related to theta modulation and precession
ThetaAna = ThetaAnalysis(Nav, Spk, Lfp, thetaparams);%Theta related analysis on spikes and decoding


