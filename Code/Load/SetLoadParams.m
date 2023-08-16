function loadparams = SetLoadParams(datadirpath)
% Define a set of parameters needed to load the data.
%
% INPUT:
% - datadirpath: path to the directory where data are stored. Default to
%   '/ Matlab Drive / School / Data' if not provided.
%
% OUTPUT:
% - loadparams: a structure whose fields contain parameters required to call
%   LoaddataNav, LoaddataSpk, and LoaddataLfp.
%
% Fields of loadparams are the following:
% - animalname: name of the animal.
% - sessionID: number identifying the session to load.
% - Datafolder: Path to the folder where the data are located, that is
%   datadirpath argument.
% - catevtfilename: file name for timestamps of each subsession.
% - pufevtfilename: file name for timestamps at which an air puff was 
%   delivered.
% - rrwevtfilename: file name for timestamps at which a reward was 
%   delivered when the animal reached the right platform.
% - lrwevtfilename: file name for timestamps at which a reward was 
%   delivered when the animal reached the left platform.
% - ripevtfilename: file name for timestamps at which a ripple had been 
%   detected.
% - posfilename: file name for tracked positions.
% - laptypefilename: file name for the type of trial (+1 if the animal
%   went from the left to the right platform; -1 if it went from the right
%   to the left one and 0 if made a u-turn and returned to the same
%   platform.
% - accfilename: file name for the accelerometer data (3-axis data).
% - pix2cm: Size of the pixels of the video tracking in cm.
% - sampleRate: Final sampling rate at which Behavioral and spiking data
%   will be resampled.
% - spkfilename: file name for the spike timestamps and clusters ID.
% - spkinfofilename: file name for information about each cluster (shank
%   ID, Pyramidal / Interneuron, etc).
% - ShankList: List of shank numbers from which spikes will be loaded 
%   (defaults to all shanks located in the hippocampus).
% - ShankList_hpc: shank numbers located in the hippocampus.
% - ShankList_blaL: shank numbers located in the left amygdala.
% - ShankList_blaR: shank numbers located in the right amygdala.
% - lfpfilename: filename for raw Lfp data.
% - LfpChannel_Hpc: Lfp channel to load from hippocampal recording.
% - LfpChannel_Bla: Lfp channel to load from amygdala recording.
% - sampleRate_raw: sampling rate at which raw Lfp data will be 
%   resampled.
% - ThetaBand: frequency range of the theta band [minFq maxFq]. Default is
%   from 6 to 9 Hz.
% - statefilename: file name for start and end timestamps of the state of
%   the animal (awake, drowsy, REM-sleep, SW-sleep).
%
% USAGE:   
% datadirpath = <path to the directory containing your data>
% loadparams = SetLoadParams(datadirpath);
%
% See also: LoaddataNav, LoaddataSpk, LoaddataLfp
%
% Written by J.Fournier in 08/2023 for the Summer school "Advanced
% computational analysis for behavioral and neurophysiological recordings"

%animal name and session ID
loadparams.animalname = 'Rat08';
loadparams.session = 20130713;

%directory where data are stored
if nargin == 1
    loadparams.Datafolder = datadirpath;
else
    loadparams.Datafolder = [filesep 'MATLAB Drive' filesep 'School' filesep 'Dataset' filesep 'Session1'];
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
loadparams.laptypefilename = 'LapType.mat';
loadparams.spkfilename = 'AllSpikes.mat';
loadparams.spkinfofilename = 'IndexType.mat';
loadparams.lfpfilename = 'LFP.mat';
loadparams.accfilename = 'LFP.mat';
loadparams.statefilename = 'States.mat';

%Final smapling rate at which Behavioral, spiking and theta data will be 
%resampled
loadparams.sampleRate = 50;

%Final smapling rate at which raw Lfp data will be resampled.
loadparams.sampleRate_raw = 600;

%Size of the pixels of the video tracking in cm
loadparams.pix2cm = 0.43;

%List of shanks located in the hippocampus
loadparams.ShankList= 1:4;

%List of shanks located in the hippocampus
loadparams.ShankList_hpc = 1:4;

%List of shanks located in the left bla
loadparams.ShankList_blaL = 5:8;

%List of shanks located in the right bla
loadparams.ShankList_blaR = 13:19;

%LFP channel to load from hippocampal recording
loadparams.LfpChannel_Hpc = 2;

%LFP channel  to load from amygdala recording
loadparams.LfpChannel_Bla = 2;

%Frequency range to filter Lfp in the theta band
loadparams.ThetaBand = [6 9];
end