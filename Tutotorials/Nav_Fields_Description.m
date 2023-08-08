% Fields of Nav are the following:
%
% - sampleTimes: time stamps of the samples for behavioral data
%
% - X: positions of the animal on the X axis
%
% - Y: position of the animal on ht eY axis
%
% - Xpos: position along the X axis in percentage of the track length.
%
% - XDir: direction of movement along the X axis (+1 for left to right, -1
% for right to left)
%
% - Spd: speed of movement
%
% - smthSpd: speed smoothed by a box car window
%
% - laptype: +1 if the animal went from the left to the right platform; 
%            -1 if it went from the right one to the left one;
%             0 if it went back to the same plateform before reaching the
%            end of the track
%
% - uturn: +1 on trial where the animal returned to the same platform; 0
% otherwise
%
% - trialID: trial number. Trial were defined as any continuous period of 
%   time where the animal was on the track
%
% - reward: +1 when a reward was delivred; 0 othewise.
%
% - airpuff: +1 when an air puff was delivred; 0 othewise.
%
% - state: +1 for awake; 0 for drowsy; -1 for REM sleep; -2 for slow wave 
%   sleep
%
% - acc: 3-axis accelerometer data (ntimes x 3 array)
%
% All fields of Nav have time along the first dimension.