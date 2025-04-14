%**************************************************************************
%*****  Resample and animate nakhla simulation output  ********************
%**************************************************************************

clear; close all

%*****  User input parameters and options  ********************************

RunID        = '2D_DEMO_d03';       % run identifier
start        =  0;              % first output frame
stop         =  115;            % last output frame 

field        = 'c_oxd';           % field name as stored in data files
fieldtitle   = 'SiO$_2$ $\bar{c}$'; % field name as displayed over animated figure
fieldoffset  =  0;         % offset to remove from field before plotting
fieldscale   =  1;              % field scale
fieldunits   =  '[wt$\%$]';  % field units
fieldlimits  =  [];    % field colorbar limits (leave empty [] for unclipped colorbar)
fieldindex   =  1;              % set to select compositional components

timescale    =  3600;           % time scale
timeunits    =  '[hr]';         % time units
lengthscale  =  1;              % length scale
lengthunits  =  '[m]';          % length units

framerate    = 30;              % frame rate for video playback [1/s]
framenumber  = 15*framerate;    % number of frames to resample to [1]


%*****  Prepare plotting parameters and options  **************************

% Load run parameter file
load(['../out/',RunID,'/',RunID,'_par.mat']);

% Set plotting macros
TX = {'Interpreter','Latex'}; FS = {'FontSize',15};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',12};

% Run init to initialise all variables, coefficients, and auxiliary fields
postprc = 1;
run('../src/init.m');

% Preallocate arrays for simulation times and fields
numstp    = stop-start;
time_all  = zeros(numstp, 1);
field_all = zeros(Nz, Nx, numstp);


%*****  Load and resample simulation output time series  ******************

for k = start:stop
    load(['../out/',RunID,'/',RunID,'_',int2str(k),'.mat']);
    time_all(max(1,k)) = time/timescale;
    run('../src/update.m');
    field_value = evalin('base',field);
    field_all(:, :, max(1,k)) = (field_value(:,:,fieldindex)-fieldoffset)/fieldscale;
end

% Create uniformly spaced time vector from the start to end simulation time
t_uniform = linspace(time_all(1),time_all(end),framenumber);

% Resample field data using linear interpolation
field_all_reshaped = reshape(field_all, Nz*Nx, numstp)'; % Size: [numstp x (Nz*Nx)]
field_resampled = interp1(time_all, field_all_reshaped, t_uniform, 'linear');  % [framenumber x (Nz*Nx)]


%*****  Plot field and write to video file  *******************************

% Set up the VideoWriter object for output
videoFile = ['../out/',RunID,'/',RunID,'_',field,num2str(fieldindex)];
v = VideoWriter(videoFile, 'MPEG-4');
v.FrameRate = framerate;
open(v);

% Create a figure (set to 'off' to avoid display during processing)
fig = figure('Visible', 'on');  colormap(colmap);

% Loop over the uniformly spaced times, rendering and recording each frame
for i = 1:length(t_uniform)
    % Reshape row vector back to the original 2D field dimensions
    field_frame = reshape(field_resampled(i, :), Nz, Nx);
    
    % Render the 2D field as scaled colormap image
    imagesc(Xc/lengthscale,Zc/lengthscale,field_frame); axis ij equal tight; box on; 
    cb = colorbar; if ~isempty(fieldlimits); clim(fieldlimits); end
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); 
    title([fieldtitle,' ',fieldunits,'; ',sprintf('t = %.2f %s', t_uniform(i),timeunits)],TX{:},FS{:}); 
    xlabel(['Width ',lengthunits],TX{:},FS{:}); 
    ylabel(['Depth ',lengthunits],TX{:},FS{:}); 
    
    drawnow;
    
    % Capture the current frame and write it to the video file
    frame = getframe(fig);
    writeVideo(v, frame);
end

% Finish and close the video file.
close(v);

fprintf('Animation saved to %s\n', videoFile);