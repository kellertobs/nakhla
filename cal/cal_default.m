% specify melting model composition parameters
cal.nc     =  8;
cal.oxdStr =  {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O'};
cal.oxd    = [  42.7    0.0    0.0    0.0   57.3    0.0    0.0    0.0      % for
                30.0    0.6    0.0   69.4    0.0    0.0    0.0    0.0      % fay
                54.9    0.2    7.5    3.0   33.9    0.4    0.1    0.0      % opx
                49.2    1.0    1.0   30.0    1.7   15.0    1.5    0.6      % cpx
                44.4    0.0   35.8    0.0    0.0   19.2    0.6    0.0      % ant
                67.3    0.0   20.2    0.0    0.0    0.8   11.0    0.7      % alb
                64.8    0.0   18.3    0.0    0.0    0.0    0.9   16.0      % kfs
                99.99   0.01   0.0    0.0    0.0    0.0    0.0    0.0];    % qtz

cal.cmpStr = {'for','fay','opx','cpx','ant','alb','kfs','qtz'};
cal.cmp    = [  97.5    2.5    0.0    0.0    0.0    0.0    0.0    0.0      % cphs0 => dunite
                 8.0    1.5   17.0   28.0   37.5    8.0    0.0    0.0      % perCx => pyroxenite
                 3.0    1.0   11.0   37.0   27.0   19.5    1.0    0.5      % perCm => basalt
                 0.0    0.0    0.0   13.0    6.0   12.0   31.0   38.0];    % cphs1 => rhyolite

cal.oxd = cal.oxd./sum(cal.oxd,2)*100;
cal.cmp = cal.cmp./sum(cal.cmp,2)*100;

cal.cmp_oxd = cal.cmp*cal.oxd./100;

% specify melting model phase diagram parameters
cal.cphs0    =  cal.cmp_oxd(1,1)/100; % phase diagram lower bound composition [wt SiO2]
cal.cphs1    =  cal.cmp_oxd(4,1)/100; % phase diagram upper bound composition [wt SiO2]
cal.Tphs0    =  850;                  % phase diagram lower bound temperature [degC]
cal.Tphs1    =  1850;                 % phase diagram upper bound temperature [degC]
cal.PhDg     =  [8.0,3.5,1.3,1.3];    % phase diagram curvature factor (> 1)
cal.perCm    =  cal.cmp_oxd(3,1)/100; % peritectic liquidus composition [wt SiO2]
cal.perCx    =  cal.cmp_oxd(2,1)/100; % peritectic solidus  composition [wt SiO2]
cal.perT     =  1125;                 % peritectic temperature [degC]
cal.clap     =  1e-7;                 % Clapeyron slope for P-dependence of melting T [degC/Pa]
cal.dTH2O    =  [1400,1300,1200];     % solidus shift from water content [degC/wt^0.75]

% specify geochemical model parameters
cal.nte      =  4;           % number of trace elements
cal.nir      =  2;           % number of isotope ratios
cal.Kte_cmp  =  [0.01  0.01  0.01  0.01  0.01  0.01  0.01  0.01  % component-wise trace element partition coefficients
                 0.10  0.10  0.10  0.10  0.10  0.10  0.10  0.10
                 3.00  3.00  3.00  3.00  3.00  3.00  3.00  3.00
                 10.0  10.0  10.0  10.0  10.0  10.0  10.0  10.0];

% specify density parameters
cal.rhox0 = [3270,4390,3000,3250,2730,2620,2580,2650];                     % solid  component reference densities [kg/m3]
cal.rhom0 = [2710,3580,2580,2850,2530,2310,2290,2360];                     % liquid component reference densities [kg/m3]
cal.rhof0 =  500;                                                          % fluid reference density [kg/m3]
cal.aT    =  4e-5;                                                         % thermal expansivity [1/K]
cal.gH    =  0.75;                                                         % hydrous melt expansivity [1/(wt H2O)]
cal.bP    =  1e-8;                                                         % fluid compressibility [1/Pa]

% specify mixture viscosity parameters (Costa et al., 2009)
cal.Bphi    = 2.0;              % Einstein-Roscoe powerlaw coefficient bubbles
cal.Bchi    = 2.0;              % Einstein-Roscoe powerlaw coefficient crystals
cal.chi_pck = 0.60;             % rheologically critical crystal fraction
cal.gamma   = 2.50;             % step-function steepness coefficient
cal.delta   = 27;               % solid viscosity melt-weakening slope
cal.xi      = 4.5e-4;           % solid viscosity level
cal.etaf0   = 0.1;              % fluid viscosity constant

% specify segregation coefficient parameters
cal.bm      = 100;              % melt permeability geometric factor (k0 = dx^2/bm)
cal.cm      = 0.001;            % melt percolation threshold
cal.nm      = 3;                % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.mm      = 2;                % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.bf      = 100;              % fluid permeability geometric factor (k0 = dx^2/bm)
cal.cf      = 0.1;              % fluid percolation threshold
cal.nf      = 4;                % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
cal.mf      = 2;                % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
