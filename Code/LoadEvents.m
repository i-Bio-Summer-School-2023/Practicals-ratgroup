function events = LoadEvents(filepath)
%Load the .evt file located at filepath.
%
%INPUT:
% - filepath: full path to the .evt file.
%OUTPUT:
% A structure with the following fields:
% - timestamps: a list of timestamps in seconds
% - description: a cell array of string describing each timestamps
%
% Usage:
% events = LoadEvents(filepath);
%
%
% by J. Fournier 07/2023

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