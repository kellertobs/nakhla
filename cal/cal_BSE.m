% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 7;
cal.nmem   = 12;
cal.nmsy   = 6;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','hyp','fsl','mau','fau','mgt','ant','alb','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'cmp1','cmp2','cmp3','cmp4','cmp5','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1          3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   Al2O3     FeO     MgO     CaO    Na2O     H2O
cal.mem_oxd    = [  42.70    0.0     0.0    57.3     0.0     0.0     0.0     % forsterite (for)
                    29.50    0.0    70.5     0.0     0.0     0.0     0.0     % fayalite (fay)

                    56.16    2.06    8.27   33.36    0.15    0.0     0.0     % enstatite (ens)
                    52.10    3.78   17.85   24.79    1.48    0.0     0.0     % hypersthene (hyp)
                    48.14    2.29   34.46   11.47    3.64    0.0     0.0     % ferrosilite (hyp)

                    51.47    1.63   14.73   12.60   19.20    0.37    0.0     % Mg-augite (mau)
                    49.54    0.78   32.01    3.02   13.31    1.34    0.0     % Fe-augite (fau)

                     0.0     1.84   97.40    0.76    0.0     0.0     0.0     % magnetite (mgt)

                    44.4    35.8     0.0     0.0    19.5     0.3     0.0     % anorthite (ant)
                    68.0    20.0     0.0     0.0     0.5    11.5     0.0     % albite (alb)

                    100.0    0.0    0.0     0.0     0.0      0.0     0.0     % quartz (qtz)

                     0.0     0.0     0.0     0.0     0.0     0.0   100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                    for       fay       ens       hyp       fsl       mau       fau       mgt       ant       alb       qtz       wat
cal.cmp_mem =   [  100.0000         0         0         0         0         0         0         0         0         0         0         0
   68.7993   31.2007         0         0         0         0         0         0         0         0         0         0
         0   14.5954   38.5037         0         0         0         0         0   46.9009         0         0         0
         0         0         0    2.0009         0   79.2930         0         0    6.9138   11.7924         0         0
         0    0.5056         0         0    0.4245         0   27.1232    0.9959   14.4074    5.5379   51.0055         0
         0         0         0         0         0         0         0         0         0         0         0  100.0000];
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
cal.T0  = [1890.0; 1150; 1050; 950; 850];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./300;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = 0*cal.A + 1;

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [70; 25; 20; 18; 25];

% specify melting point dependence on H2O
cal.dTH2O   = 1500;                 % solidus shift from water content [degC/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [degC/wt^pH2O]

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
