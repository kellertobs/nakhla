% set run parameters
runID    =  'default';           % run identifier
opdir    =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  1;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100;                 % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/N;                 % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  3e5;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1*yr;                % end time for simulation [s]
dt       =  10;                  % initial time step [s]
dtmax    =  1e3;                 % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  10;                  % regularisation of initial random perturbation
zlay     =  2.0;                 % layer thickness (relative to domain depth D)
dlay     =  0.0;                 % random perturbation to layer thickness (relative to grid spacing h)
wlay_T   =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1150;                % temperature top layer [deg C]
T1       =  1150;                % temperature base layer [deg C]
dTr      =  0;                   % amplitude of random noise [deg C]
dTg      =  0;                   % amplitude of centred gaussian [deg C]
c0       =  [0.14,0.49,0.28,0.09,0.02];  % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  [0.14,0.49,0.28,0.09,0.02];  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0];         % amplitude of random noise [wt SiO2]
dcg      =  [0,0,0,0,0];         % amplitude of centred gaussian [wt SiO2]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
te0      =  [1,1,1,1];           % trace elements top layer [wt ppm]
te1      =  [1,1,1,1];           % trace elements base layer [wt ppm]
dter     =  [0,0,0,0];           % trace elements random noise [wt ppm]
dteg     =  [0,0,0,0];           % trace elements centred gaussian [wt ppm]
ir0      =  [1, 1];              % isotope ratios top layer [delta]
ir1      =  [1, 1];              % isotope ratios base layer [delta]
dirr     =  [0, 0];              % isotope ratios random noise [delta]
dirg     =  [0, 0];              % isotope ratios centred gaussian [delta]

% set thermo-chemical boundary parameters
Ptop     =  125e6;               % top pressure [Pa]
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
bnd_h    =  [0,0,0];             % internal wall rock layer thickness [m]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
tau_a    =  24*hr;               % wall cooling/assimilation tie [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  [nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan; ...
             nan,nan,nan,nan,nan]; % [top,bot,sds] wall rock major component [wt SiO2] (nan = no assimilation)
tewall   =  [nan,nan,nan,nan; ...
             nan,nan,nan,nan; ...
             nan,nan,nan,nan];   % [top,bot,sds] wall rock trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan; ...
             nan,nan; ...
             nan,nan];           % [top,bot,sds] wall rock isotope ratios [delta] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'andesSVZ';          % phase diagram calibration
kT0      =  4;                   % thermal conductivity [W/m/K]
cP       =  1200;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]
tau_r    =  0;                   % reaction time scale (set to zero for quasi-equilibrium mode)

% set model buoyancy parameters
aT       =  4e-5;                % thermal expansivity [1/K]
bPx      =  1e-11;               % fluid compressibility [1/Pa]
bPm      =  3e-11;               % fluid compressibility [1/Pa]
bPf      =  1e-8;                % fluid compressibility [1/Pa]
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
rtol     =  1e-6;                % outer its relative tolerance
atol     =  1e-9;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
alpha    =  0.6;                 % iterative step size parameter
beta     =  0.1;                 % iterative damping parameter
etacntr  =  1e+6;                % maximum viscosity contrast
Delta    =  2*D/50;              % correlation length for eddy viscosity
Prt      =  10;                  % turbulent Prandtl number (ratio of momentum to heat diffusivity)
mink     =  1e-8;                % minimum diffusivity for phase, component fractions

% set various options
calibrt  =  0;                   % not in calibrate mode
TINY     =  1e-16;               % minimum cutoff phase, component fractions
bnchm    =  0;                   % not a benchmark run

