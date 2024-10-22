% set run parameters
runID    =  'default';           % run identifier
srcdir   =  '../src';            % output directory
outdir   =  '../out';            % output directory
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
init_mode= 'layer';              % T initial condition mode ('layer' or 'linear')
seed     =  24;                  % random perturbation seed
smth     =  10;                  % regularisation of initial random perturbation
zlay     =  2.0;                 % layer thickness (relative to domain depth D)
dlay     =  0.0;                 % random perturbation to layer thickness (relative to grid spacing h)
wlay_T   =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
wlay_c   =  0*h/D;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  1175;                % temperature top layer [deg C]
T1       =  1175;                % temperature base layer [deg C]
dTr      =  0;                   % amplitude of random noise [deg C]
dTg      =  0;                   % amplitude of centred gaussian [deg C]
c0       =  [0.04,0.12,0.44,0.24,0.14,0.02,0.005];  % components (maj comp, H2O) top layer [wt] (will be normalised to unit sum!)
c1       =  [0.04,0.12,0.44,0.24,0.14,0.02,0.005];  % components (maj comp, H2O) bot layer [wt] (will be normalised to unit sum!)
dcr      =  [0,0,0,0,0,0,0];     % amplitude of random noise [wt SiO2]
dcg      =  [0,0,0,0,0,0,0];     % amplitude of centred gaussian [wt SiO2]

% set model trace and isotope geochemistry parameters (must match # trace elements and isotope ratios in calibration!)
trc0     =  [1,1,1,1,1,1];       % trace elements top layer [wt ppm]
trc1     =  [1,1,1,1,1,1];       % trace elements base layer [wt ppm]
dr_trc   =  [0,0,0,0,0,0];       % trace elements random noise [wt ppm]
dg_trc   =  [0,0,0,0,0,0];       % trace elements centred gaussian [wt ppm]

% set thermo-chemical boundary parameters
fractxtl =  0;                   % fractional crystallisation mode for 0-D (Nz=Nx=1)
fractmlt =  0;                   % fractional melting mode for 0-D (Nz=Nx=1)
fractres =  0.25;                % residual fraction for fractionation mode
dPdT     =  0e5;                 % decompression rate for 0D models
Ptop     =  125e6;               % top pressure [Pa]
periodic =  0;                   % set side boundaries to periodic
bndmode  =  3;                   % boundary assimilation mode (0 = none; 1 = top only; 2 = bot only; 3 = top/bot only; 4 = all walls; 5 = only sides)
bnd_w    =  h;                   % boundary layer width [m]
bnd_h    =  [0,0,0];             % internal wall rock layer thickness [m]
fin      =  0;                   % ingassing factor (0 = no ingassing; 1 = free flow ingassing)
fout     =  1;                   % outgassing factor (0 = no outgassing; 1 = free flow outgassing)
tau_T    =  12*hr;               % wall cooling/assimilation time [s]
tau_a    =  24*hr;               % wall cooling/assimilation tie [s]
Twall    =  [300,300,nan];       % [top,bot,sds] wall rock temperature [degC] (nan = insulating)
cwall    =  nan(3,7);            % [top,bot,sds] wall rock major component [wt SiO2] (nan = no assimilation)
trcwall  =  nan(3,6);            % [top,bot,sds] wall rock trace elements [wt ppm] (nan = no assimilation)

% set thermo-chemical material parameters
calID    =  'MORB';              % phase diagram calibration
aTm      =  5e-5;                % melt  thermal expansivity [1/K]
aTx      =  1e-5;                % xtal  thermal expansivity [1/K]
aTf      =  1e-4;                % fluid thermal expansivity [1/K]
kTm      =  2;                   % melt  thermal conductivity [W/m/K]
kTx      =  5;                   % xtal  thermal conductivity [W/m/K]
kTf      =  0.1;                 % fluid thermal conductivity [W/m/K]
cPm      =  1200;                % melt  heat capacity [J/kg/K]
cPx      =  1000;                % xtal  heat capacity [J/kg/K]
cPf      =  2000;                % fluid heat capacity [J/kg/K]
tau_r    =  0;                   % reaction time scale (set to zero for quasi-equilibrium mode)

% set model buoyancy and pressure parameters
bPx      =  1e-11;               % solid compressibility [1/Pa]
bPm      =  3e-11;               % melt  compressibility [1/Pa]
bPf      =  1e-9;                % fluid compressibility [1/Pa]
dm0      =  1e-3;                % melt film size [m]
dx0      =  1e-3;                % crystal size [m]
df0      =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set chamber pressure parameters
Pchmb0   =  0;                   % initial chamber pressure [Pa]
eta_wall =  1e15;                % wall rock viscosity [Pas]
mod_wall =  1e10;                % wall rock elastic modulus [Pa]

% set numerical model parameters
TINT     =  'bd2im';             % time integration scheme ('be1im','bd2im','cn2si','bd2si')
ADVN     =  'weno5';             % advection scheme ('centr','upw1','quick','fromm','weno3','weno5','tvdim')
CFL      =  1.00;                % (physical) time stepping courant number (multiplies stable step) [0,1]
maxcmp   =  0.01;                % maximum change in phase fraction due to compaction
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-8;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
alpha    =  0.75;                % iterative step size parameter
beta     =  0.00;                % iterative damping parameter
gamma    =  1e-3;                % artificial horizontal inertia parameter (only applies if periodic)
lambda1  =  0e-7;                % pressure regularisation parameter
lambda2  =  0e-7;                % pressure regularisation parameter
etacntr  =  1e+8;                % maximum viscosity contrast
Delta_cnv=  h/2;                 % correlation length for eddy diffusivity (multiple of h, 0.5-1)
Delta_sgr=  dx0*10;              % correlation length for phase fluctuation diffusivity (multiple of dx0, df0, 10-20)
Prt      =  3;                   % turbulent Prandtl number (ratio of momentum to heat diffusivity)
Sct      =  3;                   % turbulent Schmidt number (ratio of momentum to mass diffusivity)
etamin   =  0.01;                % minimum viscosity
kmin     =  1e-9;                % minimum diffusivity
kmax     =  1e+9;                % maximum diffusivity
Pcouple  =  0;                   % coupling phase equilibria and material properties to dynamic pressure

% set various options
calibrt  =  0;                   % not in calibrate mode
bnchm    =  0;                   % not a benchmark run

