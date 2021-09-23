clear; close all;

% set run parameters
runID    =  'demo';             % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  1;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to live plot results
save_op  =  0;                   % switch on (1) to save output to file
plot_cv  =  0;                   % switch on (1) to live plot iterative convergence
isotherm =  0;
isochem  =  0;
diseq    =  0;                   % disequilibrium phase evolution

% set model domain parameters
D        =  5;                   % chamber depth [m]
L        =  5;                   % chamber width [m]
N        =  100 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]
dw       =  0.05;                % boundary layer thickness for cooling/outgassing/assimilation [m]

% set model timing parameters
M        =  1e5;                 % number of time steps to take
tend     =  3600*24;             % end time for simulation [s]
dt       =  1;                   % initial time step [s]
hr       =  3600;                % conversion seconds to hours

% set model vesicularity parameters
ftop     =  nan;                 % top outgassing vesicularity [vol] (nan = no outgassing)
tau_f    =  1*hr;                % top outgassing time [s]
kf       =  1e-9;                % volume diffusivity [m^2/s]
smth     =  (N/25)^2;            % regularisation of initial random perturbation

% set model thermo-chemical parameters
c0       =  0.57;                % initial magma composition [wt% SiO2]
c1       =  0;                   % amplitude of random noise [wt% SiO2]
v0       =  0;                   % initial magma water content [wt% H2O]
v1       =  0;                   % amplitude of random noise [wt% H2O]
T0       =  1275;                % initial magma temperature [degC]
T1       =  1;                   % amplitude of random noise [degC]
Ptop     =  1e8;                 % top pressure [Pa]
ctop     =  nan;                 % top composition [wt% SiO2] (nan = no assimilation)
tau_c    =  48*hr;               % chamber wall assimilation time [s]
Twall    =  500;                 % wall temperature [degC] (nan = insulating)
coolmode =  3;                   % mode of wall cooling (0 = no cooling; 1 = top only; 2 = top/bot only; 3 = all walls)
tau_T    =  12*hr;               % chamber wall cooling time [s]
kc       =  1e-7;                % chemical diffusivity [m^2/s]
kc       =  1e-7;                % chemical diffusivity [m^2/s]
kTm      =  4;                   % melt thermal diffusivity [m2/s]
kTx      =  1;                   % xtal thermal diffusivity [m2/s]
kTf      =  0.02;                % mvp  thermal diffusivity [m2/s]
Cpm      =  1400;                % melt heat capacity [J/kg/K]
Cpx      =  1000;                % xtal heat capacity [J/kg/K]
Cpf      =  2000;                % mvp  heat capacity [J/kg/K]

% set model phase equilibrium parameters
cphs0    =  0.35;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.70;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  5.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.60;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.525;               % peritectic solidus  composition [wt SiO2]
perT     =  1050;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  40;                  % solidus shift from water content [degC/wt]
tau_r    =  20;                  % crystallisation time [s]
DLx      = -300e3;               % latent heat [J/kg]
DLf      =  400e3;               % latent heat [J/kg]

% set model rheology parameters
eta0     =  1e3;                 % background melt viscosity [Pas]
A        = -1.0;                 % bubble weakening exponent
B        =  2.0;                 % crystal stiffening exponent

% set model buoyancy parameters
rhom     =  2400;                % melt phase density [kg/m3]
rhox     =  3000;                % crystal phase density [kg/m3] 
rhof     =  250;                 % bubble phase density [kg/m3]
dx       =  0.000;               % crystal size [m]
df       =  0.000;               % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.25;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'UPW2';              % advection scheme ('UPW2', 'UPW3', or 'FRM')
theta    =  1.00;                % time-stepping scheme selector (1=BE, 1/2=CN, 0=FE)
rtol     =  1e-4;                % outer its relative tolerance
atol     =  1e-7;                % outer its absolute tolerance
maxit    =  10;                  % maximum outer its
alpha    =  0.0;                 % iterative lag parameter phase equilibrium
beta     =  0.5;                 % iterative lag parameter reaction rates
gamma    =  0.0;                 % numerical compressibility [0,1]
delta    =  0;                   % regularisation of viscosity
etactr   =  1e8;                 % minimum viscosity for regularisation
TINY     =  0;                   % minimum cutoff phase, component fractions

% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
    parfile = ['../out/',runID,'/',runID,'_par'];
    save(parfile);
end

% run code
run('../src/main')
