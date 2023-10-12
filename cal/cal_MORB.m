% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 11;
cal.nmsy   = 5;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ant','alb','dps','aug','tma','tim','ilm','qtz','wat'};
cal.msyStr = {'olv','fsp','pxn','oxs','qtz'};
cal.cmpStr = {'dun','tro','ogb','gbn','fbs','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2    TiO2    Al2O3    FeO     MgO     CaO    Na2O     H2O
cal.mem_oxd    = [  42.78    0.0     0.0     0.07   57.15    0.0     0.0     0.0     % forsterite (for)
                    35.79    0.0     0.0    37.15   27.06    0.00    0.0     0.0     % fayalite (fay)

                    44.18    0.0    35.97    0.0     0.0    19.37    0.48    0.0     % anorthite (ant)
                    67.05    0.0    20.62    0.0     0.0     1.37   10.96    0.0     % albite (alb)

                    53.92    0.0     2.86    2.00   20.01   20.94    0.27    0.0     % diopside (dps)
                    50.10    0.0     0.56   34.22    8.43    5.42    1.27    0.0     % augite (aug)

                     0.0    37.96    2.18   33.29   26.57    0.0     0.0     0.0     % Ti-Mg-Al-bearing spinel (tma)
                     0.0    20.15    1.08   76.39    2.38    0.0     0.0     0.0     % Ti-magnetite (tim)
                     0.0    50.53    0.0    49.47    0.0     0.0     0.0     0.0     % ilmenite (ilm)

                   100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     % quartz (qtz)
                     0.0     0.0     0.0     0.0     0.0     0.0     0.0   100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  1  1  0  0  0  0  0    % pyroxene (pxn)
               0  0  0  0  0  0  1  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                  for     fay     ant     alb     dps     aug     tma     tim     ilm     qtz     wat
cal.cmp_mem =   [83.78   16.22       0   16.22       0       0       0       0       0       0       0
                 27.10    7.42       0    7.42       0       0       0       0       0       0       0
                  2.83    2.73   40.34    2.73   40.34   40.34    3.56       0       0       0       0
                  2.83    2.73   40.34    2.73   40.34   40.34    3.56       0       0       0       0
                  2.80    9.48    0.58    9.48    0.58    0.58   40.86    3.92    0.36       0       0
                     0    1.38    3.39    1.38    3.39    3.39   12.44    2.41    0.13   52.52       0
                     0       0       0       0       0       0       0       0       0       0  100.00];
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
cal.T0  = [1850  1250  1180  1140  1051  799];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [9 5 4 3.5 3 2.5];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [26.2  2.0  5.9  10.0  12.9  5.4];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;0.9;3.00;10.0;1.1].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  dps  aug  mgt  ant  alb  qtz  wat
cal.rhox0   = [3270,4390,3220,3520,4850,4950,2680,2570,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
cal.etax0 = 1e18;                   % solid reference viscosity [Pas]
cal.etaf0 = 1e-1;                   % vapour reference viscosity [Pas]
cal.Eax   = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA    =[ 0.65, 0.25, 0.35; ...  % permission slopes
             0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
             0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB    =[ 0.55, 0.18, 0.27; ...  % permission step locations
             0.64,0.012,0.348; ...  % each row sums to 1
             0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC    =[[0.30, 0.30, 0.40]*0.7; ... % permission step widths
            [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
            [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

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
