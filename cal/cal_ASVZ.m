% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 12;
cal.nmsy   = 5;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'hyp','fsl','dps','mau','fau','msp','tsp','ant','alb','san','qtz','wat'};
cal.msyStr = {'opx','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'ano','spa','gbn','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   TiO2    Al2O3     FeO     MgO     CaO   Na2O     K2O    H2O
cal.mem_oxd    = [ 51.59    0.0     3.42   19.21   22.98    2.81     0.0    0.0    0.0     % hypersthene (hyp)
                   48.83    0.0     0.38   40.90    8.88    1.01     0.0    0.0    0.0     % ferrosillite (fsl)

                   52.57    0.12    2.06    9.01   16.93   19.15    0.14    0.02   0.0     % diopside (dps)
                   52.48    0.18    0.71   17.19    8.78   18.79    1.75    0.12   0.0     % Mg-augite (mau)
                   50.39    1.18    0.63   32.49    2.75    8.92    3.31    0.33   0.0     % Fe-augite (fau)

                    0.0     2.78    3.71   82.94   10.57    0.0     0.0     0.0    0.0     % Mg-spinel (msp)
                    0.0    29.42    0.0    70.18    0.40    0.0     0.0     0.0    0.0     % Ti-spinel (tsp)

                   44.4     0.0    35.83    0.0     0.0    19.20    0.57    0.0    0.0     % anorthite (ant)
                   68.0     0.0    19.96    0.0     0.0     0.59   11.45    0.0    0.0     % albite (alb)
                   66.91    0.0    18.59    0.0     0.0     0.0     2.12   12.38   0.0     % sanidine (san)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  1  1  1  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  1  1  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   hyp     fsl     dps     mau     fau     msp     tsp     ant     alb     san     qtz     wat
cal.cmp_mem =   [     0       0       0       0       0       0       0 100.000       0       0       0       0
                      0       0       0       0       0  21.931   3.779  65.030   9.260       0       0       0
                 27.366       0  29.188       0       0   7.510   4.269  24.440   7.227       0       0       0
                      0   2.999       0   5.126       0   0.537   1.653       0  81.653   8.032       0       0
                      0   1.650       0       0   3.876       0   0.085       0  23.716  30.807  39.866       0
                      0       0       0       0       0       0       0       0       0       0       0 100.000];
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
cal.T0 =  [1553.0  1183.2  1112.8  992.6  708.4];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A  =  (cal.T0+273.15)./300;

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  0*cal.A + 1;

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [33.31  4.00  8.10  16.40  20.4100];

% specify melting point dependence on H2O
cal.dTH2O   = 1500;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.nte     = 4;                    % number of trace elements
cal.nir     = 2;                    % number of isotope ratios
cal.Kte_mem = [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
cal.rhox0   = [3270,3600,3250,3320,3450,4550,4950,2690,2570,2520,2650,1000]; % mineral end-member reference densities [kg/m3]
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
