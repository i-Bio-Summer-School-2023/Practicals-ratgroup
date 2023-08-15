function events = LoadEvents(filepath)
% LoadEvents - Load and parse event data from a .evt file.
%
%   events = LoadEvents(filepath) loads event data from the specified .evt
%   file located at the given filepath. The function returns a MATLAB
%   structure containing timestamps and descriptions of the events.
%
% INPUT:
%   filepath - Full path to the .evt file.
%
% OUTPUT:
%   events - MATLAB structure with the following fields:
%     * timestamps - A list of timestamps in seconds.
%     * description - A cell array of strings describing each timestamp.
%
% USAGE:
%   events = LoadEvents(filepath);
%
% SEE ALSO:
%   LoaddataNav, LoaddataSpk, LoaddataLfp, SetLoadParams
%
% Written by J. Fournier in 08/2023 for the Summer school
% "Advanced computational analysis for behavioral and neurophysiological 
% recordings"


f = fopen(filepath,'r');
events.timestamps = [];
events.description = [];
i = 0;
while ~feof(f)
    i = i + 1;
    events.timestamps(i) = fscanf(f, '%f', 1);
    line = fgetl(f);
    start = regexp(line, '[^\s]', 'once');
    events.description{i} = sscanf(line(start:end),'%c');
end
fclose(f);

%convert timestamps from ms to seconds
events.timestamps = events.timestamps / 1000;

end