%**************************************************************************
%*****  Postprocess nakhla simulation output  *****************************
%**************************************************************************

clear; close all

%*****  User input parameters and options  ********************************

RunID        = '2D_DEMO_d3';       % run identifier
postdir      =  '../out';       % directory for storing postprocessed output
start        =  1;              % first output frame
stop         =  120;            % last output frame 

% Load run parameter file
load(['../out/',RunID,'/',RunID,'_par.mat']);
outdir  = postdir;
postprc = 1;

% Run init to initialise all variables, coefficients, and auxiliary fields
run('../src/init.m');

% Load output variables and postprocess all coefficients and auxiliary fields

for k = start:stop
    % load raw output frame
    load([postdir,'/',RunID,'/',RunID,'_',int2str(k),'.mat']);
    
    % run coefficient update
    run('../src/update.m');
    run('../src/phseql.m');
    run('../src/corrl.m');

    % store postprocessed output frame
    save([postdir,'/',RunID,'/',RunID,'_post_',int2str(k),'.mat']);

end

% 