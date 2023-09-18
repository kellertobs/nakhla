% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 15;
cal.nmsy   = 6;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','hyp','fsl','dps','mau','fau','ulv','mgt','ilm','ant','alb','san','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'ano','tro','gbn','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   TiO2   Al2O3     FeO     MgO     CaO    Na2O     K2O    H2O
cal.mem_oxd    = [ 42.71    0.0     0.0     0.0    57.29    0.00    0.0     0.0    0.0     % forsterite (for)
                   29.49    0.0     0.0    70.51    0.0     0.00    0.0     0.0    0.0     % fayalite (fay)

                   51.37    0.0     3.42   20.47   21.89    2.85    0.0     0.0    0.0     % hypersthene (hyp)
                   48.90    0.0     0.11   41.55    8.60    0.84    0.0     0.0    0.0     % ferrosillite (fsl)

                   52.44    0.22    2.21    8.98   15.33   20.35    0.41    0.06   0.0     % diopside (dps)
                   52.57    0.05    0.81   16.15   11.30   18.01    1.08    0.03   0.0     % Mg-augite (mau)
                   51.71    0.94    0.80   24.09    3.95   14.96    3.23    0.32   0.0     % Fe-augite (fau)

                    0.0    38.62    2.61   41.09   17.69    0.0     0.0     0.0    0.0     % ulvospinel (msp)
                    0.0     0.0     0.44   99.56    0.0     0.0     0.0     0.0    0.0     % magnetite (mgt)
                    0.0    52.00    0.00   48.00    0.0     0.0     0.0     0.0    0.0     % ilmenite (ilm)

                   44.00    0.0    36.10    0.0     0.0    19.52    0.38    0.0    0.0     % anorthite (ant)
                   66.93    0.0    20.69    0.0     0.0     1.72   10.61    0.05   0.0     % albite (alb)
                   68.11    0.0    18.19    0.0     0.0     0.0     5.52    8.18   0.0     % sanidine (san)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  1  0  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                  for     fay     hyp     fsl     dps     mau     fau     ulv     mgt     ilm     ant     alb     san     qtz     wat
cal.cmp_mem =   [    0       0       0       0       0       0       0       0       0       0   89.47   10.53       0       0       0
                 20.06   12.91       0       0       0       0       0   10.02       0       0   52.73    4.28       0       0       0
                  0.68   12.63   18.95       0   31.48       0       0    3.04    2.22       0   12.60   18.41       0       0       0
                     0    1.10       0    6.36       0    1.70    5.56       0    0.36    0.85    0.01   75.30    8.77       0       0
                     0       0       0    0.33       0       0    5.36       0       0    0.03       0    1.56   54.15   38.57       0
                     0       0       0       0       0       0       0       0       0       0       0       0       0       0  100.00];
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
cal.T0 =  [1530  1118  1070  950  749];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A  =  (cal.T0+273.15)./300;

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  0*cal.A + 1;

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [35.7  2.0  8.1  10.5  16.2];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.nte     = 4;                    % number of trace elements
cal.nir     = 2;                    % number of isotope ratios
cal.Kte_mem = [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
%              for  fay  hyp  fsl  dps  may  fay  ulv  mgt  ilm  ant  alb  san  qtz  wat
cal.rhox0   = [3270,4390,3270,3600,3250,3320,3450,4550,4980,4930,2690,2570,2520,2650,1000]; % mineral end-member reference densities [kg/m3]
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
