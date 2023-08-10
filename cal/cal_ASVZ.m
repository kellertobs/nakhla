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
cal.memStr = {'hyp','fsl','mau','fau','pig','msp','tsp','ant','alb','san','qtz','wat'};
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
cal.mem_oxd    = [ 52.75    0.0     3.55   14.47   27.01    2.22    0.0     0.0    0.0     % hypersthene (hyp)
                   50.59    0.0     0.21   34.70   13.39    1.11    0.0     0.0    0.0     % ferrosillite (fsl)

                   52.63    0.17    2.26    7.51   16.45   20.40    0.52    0.06   0.0     % Mg-augite (mau)
                   53.56    0.16    0.70   11.76   13.59   18.95    1.22    0.06   0.0     % Fe-augite (fau)
                   52.14    0.82    0.73   21.12    6.28   15.42    3.17    0.32   0.0     % pigeonite (pig)

                    0.0     0.0     4.03   84.63   11.34    0.0     0.0     0.0    0.0     % Mg-spinel (msp)
                    0.0    39.72    0.0    59.70    0.58    0.0     0.0     0.0    0.0     % Ti-spinel (tsp)

                   44.4     0.0    35.8     0.0     0.0    19.3     0.5     0.0   0.0     % anorthite (ant)
                   68.0     0.0    20.0     0.0     0.0     0.5    11.5     0.0   0.0     % albite (alb)
                   64.8     0.0    18.3     0.0     0.0     0.0     0.5    16.4   0.0     % sanidine (san)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0   0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0 100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  1  1  1  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  1  1  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                hyp       fsl       mau       fau       pig       msp       tsp       ant       alb       san       qtz     wat
cal.cmp_mem =   [         0         0         0         0         0         0         0  100.0000         0         0         0         0
                          0         0         0         0         0   19.8667    4.1042   61.7253   14.3037         0         0         0
                    19.5097    6.9877   20.2093    9.6940         0    5.5143    4.2463   25.3282    8.5105         0         0         0
                          0    6.7887         0    0.0163    6.5783         0    1.3157    3.5683   75.1501    6.5828         0         0
                          0    1.1032         0    1.2930    1.2323         0    0.2699         0   31.4133   25.3432   39.3451         0
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
cal.T0 =  [1553.0  1156.0  1106.8  990.9  748.2];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A  =  (cal.T0+273.15)./300;

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  0*cal.A + 1;

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [29.68  0.90  6.13  13.20  28.28];

% specify melting point dependence on H2O
cal.dTH2O   = 1500;                 % solidus shift from water content [degC/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [degC/wt^pH2O]

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
