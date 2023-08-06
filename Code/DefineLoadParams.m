function loadparams = DefineLoadParams(datadirpath)
%Define a set of parameters needed to load the data.
%
%INTPUT: 
% - datadirpath: path to the directory where data are stored. Default to
%   '/ Matlab Drive / Data' if not provided.
%
%OUTPUT:
% A structure with the following fields:
% - animalname: name of the animal.
%
% - sessionID: number identiying the session to load.
%
% - Datafolder: Path to the folder where the data are located, that is
%   datadirpath argument.
%
% - catevtfilename: file name for timestamps of each subsession.
%
% - pufevtfilename: file name for timestamps at which anair puff was 
%   delivered.
%
% - rrwevtfilename: file name for timestamps at which an reward was 
%   delivered when the animal reached the right platform.
%
% - lrwevtfilename: file name for timestamps at which an reward was 
%   delivered when the animal reached the left platform.
%
% - ripevtfilename: file name for timestamps at which a ripple had been 
%   detected.
%
% - posfilename: file name for tracked positions.
%
% - laptypefilename: file name for the type of trial (+1 if the animal
%   went from the left to the right platform; -1 if it went from the right
%   to the left one and 0 if made a u-turn and returned to the same
%   platform.
%
% - accfilename: file name for the accelerometer data (3 axis data).
%
% - pix2cm: Size of the pixels of the video tracking in cm.
%
% - sampleRate: Final sampling rate at which Behavioral and spiking data
%   will be resampled.
%
% - spkfilename: file name for the spike timestamps and clusters ID.
%
% - spkinfofilename: file name for information about each cluster (shank
%   ID, Pyramidal / Interneuron, etc).
%
% - ShankList: List of shank numbers from which spikes will be loaded 
%  (defaults to all shanks located in the hippocampus).
%
% - lfpfilename: filename for raw Lfp data.
%
% - LfpChannel_Hpc: Lfp channel to load from hippocampal recording.
%
% - LfpChannel_Bla: Lfp channel  to load from amygdala recording.
%
% - sampleRate_rawLfp: smapling rate at which raw Lfp data will be 
%   resampled.
%
% - ThetaBand: frequency range of the theta band [minFq maxFq]. Default is
%   from 6 to 9 Hz.
%
% - statefilename: file name for start and end timestamps of the state of
%   the animal (awake, drowsy, REM-sleep, SW-sleep).
%
% Usage:   
% loadparams = DefineLoadParams(datadirpath)
%
% written by J.Fournier 07/2023.

%animal name and session ID
loadparams.animalname = 'Rat08';
loadparams.session = 20130713;

%directory where data are stored
if nargin == 1
    loadparams.Datafolder = datadirpath;
else
    loadparams.Datafolder = [filesep 'MATLAB Drive' filesep 'Dataset' filesep 'Session1'];
end

if contains(loadparams.Datafolder,'MATLAB Drive')
    %Making sure the code is on Matlab path if it is run from the shared
    %matlad drive folders
    addpath(genpath('/MATLAB Drive/Code'));
end

%Names of the event files to load
loadparams.catevtfilename = [loadparams.animalname '-' num2str(loadparams.session) '.cat.evt'];
loadparams.pufevtfilename = [loadparams.animalname '-' num2str(loadparams.session) '.puf.evt'];
loadparams.rrwevtfilename = [loadparams.animalname '-' num2str(loadparams.session) '.rrw.evt'];
loadparams.lrwevtfilename = [loadparams.animalname '-' num2str(loadparams.session) '.lrw.evt'];
loadparams.ripevtfilename = [loadparams.animalname '-' num2str(loadparams.session) '.rip.evt'];

%Names of the .mat files to load
loadparams.posfilename = 'Positions.mat';
loadparams.laptypefilename = 'LapType2.mat';
loadparams.spkfilename = 'HippoSpikes.mat';
loadparams.spkinfofilename = 'IndexType.mat';
loadparams.lfpfilename = 'LFP2.mat';
loadparams.accfilename = 'LFP2.mat';
loadparams.statefilename = 'States.mat';

%Final smapling rate at which Behavioral, spiking and theta data will be 
%resampled
loadparams.sampleRate = 50;

%Final smapling rate at which raw Lfp data will be resampled.
loadparams.sampleRate_rawLfp = 600;

%Size of the pixels of the video tracking in cm
loadparams.pix2cm = 0.43;

%List of shank numbers from which spikes will be loaded (corresponds to
%shanks located in the hippocampus
loadparams.ShankList = 1:4;

%LFP channel to load from hippocampal recording
loadparams.LfpChannel_Hpc = 2;

%LFP channel  to load from amygdala recording
loadparams.LfpChannel_Bla = 2;

%Frequency range to filter Lfp in the theta band
loadparams.ThetaBand = [6 9];
end