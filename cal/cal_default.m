% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 11;
cal.nmsy   = 7;
cal.ncmp   = 4;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O'};
cal.memStr = {'for','fay','ens','hyp','aug','pig','ant','alb','spn','kfs','qtz'};
cal.msyStr = {'olv','opx','cpx','plg','spn','kfs','qtz'};
cal.cmpStr = {'dun','pxn','bas','rhy'};

% oxide composition of mineral end-members
cal.mem_oxd    = [  42.7    0.0    0.0    0.0   57.3    0.0    0.0    0.0      % forsterite (for)
                    30.0    0.6    0.0   69.4    0.0    0.0    0.0    0.0      % fayalite (fay)
                    56.0    0.08   4.49   4.33  33.92   1.11   0.07   0.0      % enstatite (ens)
                    50.0    0.15   8.85  22.17  16.08   2.70   0.05   0.0      % hypersthene (hyp)
                    52.0    0.08   7.16   4.88  18.05  16.03   0.20   0.6      % augite (aug)
                    48.0    1.68   2.05  33.27   4.19  10.08   0.13   0.6      % pigeonite (pig)
                    44.4    0.0   35.8    0.0    0.0   19.2    0.6    0.0      % anorthite (ant)
                    67.3    0.0   20.2    0.0    0.0    0.8   11.0    0.7      % albite (alb)
                    0.0    23.0    1.0   67.0    9.0    0.0    0.0    0.0      % spinel (spn)
                    64.8    0.0   18.3    0.0    0.0    0.0    0.9   16.0      % k-feldspar (kfs)
                    99.99   0.01   0.0    0.0    0.0    0.0    0.0    0.0];    % quartz (qtz)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  1  1  0  0  0    % plagioclase (plg)
               0  0  0  0  0  0  0  0  1  0  0    % spinel (spn)
               0  0  0  0  0  0  0  0  0  1  0    % k-feldspar (kfs)
               0  0  0  0  0  0  0  0  0  0  1];  % quartz (qtz)

% mineral end-member composition of melting model components
cal.cmp_mem =[  97.5    2.5    0.0    0.0   0.0    0.0    0.0    0.0    0.0    0.0   0.0      % cphs0 => dunite (dun)
                14.0   16.0   43.0    0.0  19.0    0.0    7.0    0.0    1.0    0.0   0.0      % perCx => pyroxenite (pxn)
                 0.1    2.0    1.0   11.0  21.0   13.0   25.0   23.0    2.0    1.0   0.9      % perCm => basalt (bas)
                 0.0    0.0    0.0    1.0   0.0   11.0    6.0   12.0    0.5   31.0  38.5];    % cphs1 => rhyolite (rhy)
cal.cmp_mem = cal.cmp_mem./sum(cal.cmp_mem,2)*100;

% mineral systems composition of melting model components
cal.cmp_msy = cal.cmp_mem*cal.msy_mem.';

% oxide composition of melting model components
cal.cmp_oxd = cal.cmp_mem*cal.mem_oxd./100;

% oxide composition of mineral systems in melting model components
for i=1:cal.ncmp
    for j=1:cal.nmsy
        cal.cmp_msy_oxd(i,j,:) = cal.cmp_mem(i,cal.msy_mem(j,:)==1)*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(cal.cmp_mem(i,cal.msy_mem(j,:)==1)+1e-32);
    end
end

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
cal.Kte_mem  =  [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
cal.rhox0 = [3270,4390,3000,3150,3250,2730,2630,2620,2580,4300,2650];      % solid  component reference densities [kg/m3]
cal.rhom0 = [2710,3580,2580,2750,2850,2530,2430,2310,2290,3700,2360];      % liquid component reference densities [kg/m3]
cal.rhof0 =  1000;                                                         % fluid reference density [kg/m3]
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
cal.bm      = 50;               % melt permeability geometric factor (k0 = dx^2/bm)
cal.cm      = 0.001;            % melt percolation threshold
cal.nm      = 3;                % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.mm      = 2;                % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.bf      = 50;               % fluid permeability geometric factor (k0 = dx^2/bm)
cal.cf      = 0.05;             % fluid percolation threshold
cal.nf      = 4;                % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
cal.mf      = 2;                % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
