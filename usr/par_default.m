% set run parameters
runID    =  'default';           % run identifier
opdir    =  '../out';            % output directory
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  100;                 % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on to live plot results
save_op  =  0;                   % switch on to save output to file
plot_cv  =  0;                   % switch on to live plot iterative convergence

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  100 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
Nt       =  3e5;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
yr       =  24*365.25*hr;        % conversion seconds to years
tend     =  1*yr;                % end time for simulation [s]
dt       =  1;                   % initial time step [s]
dtmax    =  1;                   % maximum time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  10;                  % regularisation of initial random perturbation
zlay     =  2.0;                 % layer thickness (relative to domain depth D)
wlay_T   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  2*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1250;                % temperature top layer [deg C]
T1       =  1250;                % temperature base layer [deg C]
dTr      =  0;                   % amplitude of random noise [deg C]
dTg      =  0;                   % amplitude of centred gaussian [deg C]
c0       =  0.51;                % major component top layer [wt SiO2]
c1       =  0.51;                % major component base layer [wt SiO2]
dcr      =  0e-5;                % amplitude of random noise [wt SiO2]
dcg      =  0e-5;                % amplitude of centred gaussian [wt SiO2]
v0       =  0.00;                % volatile component top layer [wt H2O]
v1       =  0.00;                % volatile component base layer [wt H2O]
dvr      =  0e-6;                % amplitude of random noise [wt H2O]
dvg      =  0e-6;                % amplitude of centred gaussian [wt H2O]

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
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls)
bnd_w    =  h;                   % boundary layer width [m]
bnd_h    =  0;                   % internal wall rock layer thickness [m]
fin      =  1;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
tau_a    =  24*hr;               % wall cooling/assimilation tie [s]
Twall    =  300;                 % wall temperature [degC] (nan = insulating)
cwall    =  nan;                 % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  nan;                 % wall volatile component [wt H2O] (nan = no assimilation)
tewall   =  [nan,nan,nan,nan];   % wall trace elements [wt ppm] (nan = no assimilation)
irwall   =  [nan,nan];           % wall isotope ratios [delta] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'default';           % phase diagram calibration
kT0      =  5;                   % thermal conductivity [W/m/K]
cP       =  1200;                % heat capacity [J/kg/K]
Dsx      = -300;                 % entropy change of crystallisation [J/kg]
Dsf      =  400;                 % entropy change of exsolution [J/kg]
tau_r    =  0;                   % reaction time scale (set to zero for quasi-equilibrium mode)

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
theta    =  1/2;                 % time stepping mode (0 explicit Euler, 1/2 Crank-Nicolson, 1 implicit Euler)
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  50;                  % maximum outer its
lambda   =  0.25;                % iterative lag parameter equilibration
etareg   =  1e0;                 % viscosity regularisation parameter
mink     =  1e-9;                % minimum diffusivity for phase, component fractions

