% Fields of mapsparams are the following:
%
% condition: experimental conditions over which place fields will be 
% estimated.
%
% dir: An array representing the directions over which place fields will 
% be estimated.
% 
% laptype: An array defining lap types over which place fields will be 
% estimated.
% 
% spdthreshold: The minimum speed threshold (in units of velocity) over 
% which place fields will be computed.
% 
% cellidx: A logical array indicating a subset of cells for which place 
% fields will be computed.
% 
% sampleRate: The sampling rate of the data (in Hz).
% 
% scalingFactor: A scaling factor applied to the response data, typically 
% set to 1 / samplingRate to convert spiking data to spikes per second.
% 
% Xvariablename: The name of the independent variable used to map the 
% response along the X-axis.
% 
% Xrange: The range of X values over which place fields will be estimated.
% 
% Xbinsize: The size of the X-axis bins.
% 
% Xsmthbinsize: The size of the Gaussian window for smoothing along the 
% X-axis.
% 
% XsmthNbins: The number of bins used for smoothing place fields along the 
% X-axis.
% 
% Xbinedges: The edges of position bins used to discretize the X-axis.
% 
% Yvariablename: The name of the independent variable used to map the 
% response along the Y-axis.
% 
% Yrange: The range of Y values over which place fields will be estimated.
% 
% Ybinsize: The size of the Y-axis bins.
% 
% Ysmthbinsize: The size of the Gaussian window for smoothing place fields
% along the Y-axis.
% 
% YsmthNbins: The number of bins used for smoothing place fields along the 
% Y-axis.
% 
% Ybinedges: The edges of position bins used to discretize the Y-axis.
% 
% occ_th: An occupancy threshold in seconds above which positions are 
% included in the place field estimate.
% 
% nspk_th: The minimal number of spikes required to consider a cell.
% 
% nShuffle: The number of shuffle controls performed for randomization.
% 
% kfold: The number of folds considered for cross-validation.