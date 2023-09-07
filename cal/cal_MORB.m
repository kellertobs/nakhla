% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 7;
cal.nmem   = 9;
cal.nmsy   = 5;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','dps','aug','mgt','ant','alb','qtz','wat'};
cal.msyStr = {'olv','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'dun','ogb','fbs','dac','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1          3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   Al2O3     FeO     MgO     CaO    Na2O     H2O
cal.mem_oxd    = [  42.71    0.0     0.0    57.29    0.00    0.0     0.0     % forsterite (for)
                    29.49    0.0    70.51    0.0     0.00    0.0     0.0     % fayalite (fay)

                    55.02    1.77    2.84   21.48   18.55    0.34    0.0     % diopside (dps)
                    49.68    0.83   32.40    0.27   14.34    2.48    0.0     % augite (aug)

                     0.0     1.55   98.04    0.41    0.0     0.0     0.0     % magnetite (mgt)

                    44.00   36.10    0.0     0.0    19.52    0.38    0.0     % anorthite (ant)
                    66.63   20.99    0.0     0.0     1.77   10.61    0.0     % albite (alb)

                   100.0     0.0     0.0     0.0     0.0     0.0     0.0     % quartz (qtz)
                     0.0     0.0     0.0     0.0     0.0     0.0   100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  1  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                    for       fay       dps       aug       mgt       ant       alb       qtz       wat
cal.cmp_mem =   [81.3484   18.6516         0         0         0         0         0         0         0
                 12.0019    1.8816   27.3737         0         0   41.1869   17.5558         0         0
                  0.0107    2.3401    0.2122   45.1515   10.6537    9.5873   32.0444         0         0
                       0    1.8042    0.2264   18.7293    2.8218    4.4565   12.2755   59.6864         0
                       0    0.0142         0         0    0.6394         0   41.1559   58.1904         0
                       0         0         0         0         0         0         0         0  100.0000];
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

% set pure component melting points T_m^i at P=0
cal.T0  = [1850.0  1140.0 1030.0 900.0 830.0];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./300;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = 0*cal.A + 1;

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [20.0  5.0  10.0  10.0  10.0];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.nte     = 4;                    % number of trace elements
cal.nir     = 2;                    % number of isotope ratios
cal.Kte_mem = [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
cal.rhox0   = [3270,4390,3220,3520,4930,2680,2570,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 1000;                 % fluid reference density [kg/m3]

% specify mixture viscosity parameters (Costa et al., 2009)
cal.Bphi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient bubbles
cal.Bchi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient crystals
cal.chi_pck = 0.60;                 % rheologically critical crystal fraction
cal.gamma   = 2.50;                 % step-function steepness coefficient
cal.delta   = 27;                   % solid viscosity melt-weakening slope
cal.xi      = 4.5e-4;               % solid viscosity level
cal.etaf0   = 0.1;                  % fluid viscosity constant

% specify segregation coefficient parameters
cal.bm      = 50;                   % melt permeability geometric factor (k0 = dx^2/bm)
cal.cm      = 0.001;                % melt percolation threshold
cal.nm      = 3;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.mm      = 2;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.bf      = 50;                   % fluid permeability geometric factor (k0 = dx^2/bm)
cal.cf      = 0.05;                 % fluid percolation threshold
cal.nf      = 4;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
cal.mf      = 2;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
