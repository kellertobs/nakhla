clear; close all;

% set run parameters
runID    =  'rchrg2';             % run identifier
restart  =  0;                   % restart from file (0: new run; <1: restart from last; >1: restart from specified frame)
nop      =  20;                  % output frame plotted/saved every 'nop' time steps
plot_op  =  1;                   % switch on (1) to live plot results
save_op  =  1;                   % switch on (1) to save output to file
plot_cv  =  1;                   % switch on (1) to live plot iterative convergence
isotherm =  0;                   % switch on (1) isothermal mode
isochem  =  0;                   % switch on (1) isochemical mode
diseq    =  0;                   % disequilibrium phase evolution

% set model domain parameters
D        =  10;                  % chamber depth [m]
L        =  10;                  % chamber width [m]
N        =  120 + 2;             % number of grid points in z-direction (incl. 2 ghosts)
h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

% set model timing parameters
M        =  2e4;                 % number of time steps to take
hr       =  3600;                % conversion seconds to hours
tend     =  hr*480;              % end time for simulation [s]
dt       =  1;                   % initial time step [s]

% set initial thermo-chemical state
seed     =  15;                  % random perturbation seed
smth     =  (N/25)^2;            % regularisation of initial random perturbation
zlay     =  0.80;                % layer thickness (relative to domain depth D)
wlay     =  1e-16;               % thickness of smooth layer boundary (relative to domain depth D)
T0       =  700;                 % temperature top layer [deg C]
T1       =  1150;                % temperature base layer [deg C]
dT       =  2;                   % amplitude of random noise [deg C]
c0       =  0.66;                % major component top layer [wt SiO2]
c1       =  0.48;                % major component base layer [wt SiO2]
dc       =  0.00;                % amplitude of random noise [wt SiO2]
v0       =  0.02;                % volatile component top layer [wt H2O]
v1       =  0.02;                % volatile component base layer [wt H2O]
dv       =  0.00;                % amplitude of random noise [wt H2O]

% set model trace and isotope geochemistry parameters
it0      =  10;                  % incompatible tracer top layer [wt ppm]
it1      =  1;                   % incompatible tracer base layer [wt ppm]
dit      =  0.0;                 % incompatible tracer random noise [wt ppm]
KIT      =  1e-2;                % incompatible tracer partition coefficient
ct0      =  1;                   % compatible tracer top layer [wt ppm]
ct1      =  10;                  % compatible tracer base layer [wt ppm]
dct      =  -0.0;                % compatible tracer random noise [wt ppm]
KCT      =  1e2;                 % compatible tracer partition coefficient
si0      =  -5;                  % stable isotope ratio top layer [delta]
si1      =  5;                   % stable isotope ratio base layer [delta]
dsi      =  0;                   % stable isotope ratio random noise [delta]

% set thermo-chemical boundary parameters
Ptop     =  1e8;                 % top pressure [Pa]
bndmode  =  0;                   % mode of wall cooling/outgassing/assimilation (0 = none; 1 = top only; 2 = top/bot only; 3 = all walls)
Twall    =  500;                 % wall temperature [degC] (nan = insulating)
dw       =  h;                   % boundary layer thickness for cooling/outgassing/assimilation [m]
tau_T    =  5*hr;                % chamber wall cooling time [s]
fwall    =  0.02;                % wall outgassing vesicularity [wt] (nan = no outgassing)
tau_f    =  1/10*hr;             % wall outgassing time [s]
cwall    =  nan;                 % wall major component [wt SiO2] (nan = no assimilation)
vwall    =  nan;                 % wall volatile component [wt H2O] (nan = no assimilation)
itwall   =  nan;                 % wall incomp. tracer [wt SiO2] (nan = no assimilation)
ctwall   =  nan;                 % wall comp. tracer [wt SiO2] (nan = no assimilation)
siwall   =  nan;                 % wall stable isotope [wt SiO2] (nan = no assimilation)
tau_c    =  48*hr;               % wall assimilation time [s]

% set thermo-chemical material parameters
kc       =  1e-8;                % chemical diffusivity [m^2/s]
kTm      =  4;                   % melt thermal diffusivity [m2/s]
kTx      =  1;                   % xtal thermal diffusivity [m2/s]
kTf      =  0.02;                % mvp  thermal diffusivity [m2/s]
Cpm      =  1400;                % melt heat capacity [J/kg/K]
Cpx      =  1000;                % xtal heat capacity [J/kg/K]
Cpf      =  2000;                % mvp  heat capacity [J/kg/K]

% set model phase equilibrium parameters
cphs0    =  0.36;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.72;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  4.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.54;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.50;                % peritectic solidus  composition [wt SiO2]
perT     =  1050;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  1300;                % solidus shift from water content [degC/wt^0.75]
tau_r    =  60;                  % crystallisation time [s]
DLx      = -400e3;               % latent heat [J/kg]
DLf      =  500e3;               % latent heat [J/kg]

% set model rheology parameters
etam0    =  1e3;                 % melt viscosity [Pas]
etaf0    =  1e-1;                % fluid viscosity [Pas]
etax0    =  1e15;                % crystal viscosity [Pas]
phic     =  0.5;                 % close packing bubble fraction [vol]
chic     =  0.5;                 % close packing crystal fraction [vol]
A        = -1.0;                 % bubble weakening exponent
B        =  2.0;                 % crystal stiffening exponent
lambda   =  30;                  % melt-weakening of crystal viscosity [1]
Fmc      =  1e+5;                % major component weakening factor of melt viscosity [1]
Fmv      =  0.3;                 % volatile component weakening factor of melt viscosity [1]
Fxc      =  1e-5;                % major component weakening factor of crystal viscosity [1]
Em       =  150e3;               % activation energy melt viscosity [J/mol]
Ex       =  250e3;               % activation energy crystal viscosity [J/mol]

% set model buoyancy parameters
rhom0    =  3000;                % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhox0    =  3300;                % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhof0    =  500;                 % bubble phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aTm      =  3e-5;                % melt thermal expansivity [1/K]
aTx      =  1e-5;                % xtal thermal expansivity [1/K]
aTf      =  1e-4;                % mvp  thermal expansivity [1/K]
gCm      =  0.5;                 % melt compositional expansion [1/wt]
gCx      =  0.6;                 % xtal compositional expansion [1/wt]
bPf      =  1e-8;                % mvp compressibility [1/Pa]
dx       =  0.001;               % crystal size [m]
df       =  0.001;               % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set numerical model parameters
CFL      =  0.50;                % (physical) time stepping courant number (multiplies stable step) [0,1]
ADVN     =  'FRM';               % advection scheme ('UPW2', 'UPW3', or 'FRM')
rtol     =  1e-3;                % outer its relative tolerance
atol     =  1e-6;                % outer its absolute tolerance
maxit    =  20;                  % maximum outer its
alpha    =  0.85;                % iterative lag parameter equilibration
delta    =  0;                   % regularisation of settling speed
etamin   =  1e2;                 % minimum viscosity for stabilisation
etamax   =  1e8;                 % maximum viscosity for stabilisation
TINY     =  1e-16;               % minimum cutoff phase, component fractions

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
