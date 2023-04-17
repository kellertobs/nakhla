% prepare workspace
clear; close all;

% load default parameters
run('./par_default')

% set run parameters
runID    =  '2D_andesSVZ_rchrg'; % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
N        =  100;                 % number of grid points in z-direction
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]
L        =  D;                   % chamber width (equal to h for 1-D mode) [m]

% set model timing parameters
Nt       =  1e5;                 % number of time steps to take
tend     =  1*yr;                % end time for simulation [s]
dt       =  36;                  % initial time step [s]
dtmax    =  360;                 % maximum time step [s]

% set initial thermo-chemical state
T0       =  765;                 % temperature top layer [deg C]
c0       =  [0.01,0.05,0.20,0.74,0.04]; % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  c0;                         % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [1,1,-1,-1,0]*1e-5;
dcg      =  [0,0,0,0,0];
te0      =  [1,0.3,0.1,0.01];    % trace elements top layer [wt ppm]
te1      =  te0;                 % trace elements base layer [wt ppm]
ir0      =  [-2,0.72];           % isotope ratios top layer [delta]
ir1      =  ir0;                 % isotope ratios base layer [delta]


% set thermo-chemical boundary parameters
bndmode  =  2;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
tau_T    =  8*hr;                % wall cooling/assimilation time [s]
tau_a    =  8*hr;                % wall cooling/assimilation tie [s]
Twall    =  [nan,1150,nan];      % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [nan,nan,nan,nan,nan; ...
             0.10,0.28,0.59,0.03,0.02; ...
             nan,nan,nan,nan,nan];
tewall   =  [nan,nan,nan,nan; ...
             2,1,0.3,0.03; ...
             nan,nan,nan,nan];   % [top,bot,sds] wall rock trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan; ...
             5, 0.70; ...
             nan,nan];           % [top,bot,sds] wall rock isotope ratios [delta] (nan = no assimilation)
Ptop     =  1.25e8;              % top pressure [Pa]

% set thermo-chemical material parameters
calID    =  'andesSVZ';          % phase diagram calibration
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]

% set numerical model parameters
TINT     =  'bd2si';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.25;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
cnvreg   =  1;                   % convection regularisation parameter


%*****  RUN NAKHLA MODEL  *************************************************
run('../src/main')
%**************************************************************************

