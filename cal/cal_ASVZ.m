% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 15;
cal.nmsy   = 6;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'ant','alb','san','for','fay','tms','mgt','dps','mau','fau','hyp','fsl','ilm','qtz','wat'};
cal.msyStr = {'fsp','olv','oxs','cxp','opx','qtz'};
cal.cmpStr = {'ano','str','osg','fbs','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   TiO2   Al2O3     FeO     MgO     CaO    Na2O     K2O    H2O
cal.mem_oxd    = [ 44.40    0.0    35.84    0.0     0.0    19.20    0.56    0.0    0.0     % anorthite (ant)
                   67.39    0.0    20.35    0.0     0.0     1.07   11.19    0.0    0.0     % albite (alb)
                   65.71    0.0    18.59    0.0     0.0     0.0     2.82   12.88   0.0     % sanidine (san)

                   38.83    0.0     0.0    20.97   40.21    0.0     0.0     0.0    0.0     % forsterite (for)
                   32.28    0.0     0.0    55.74   11.97    0.0     0.0     0.0    0.0     % fayalite (fay)

                    0.0    41.15    0.0    40.02   18.83    0.0     0.0     0.0    0.0     % Ti-Mg-bearing spinel (tms)
                    0.0     0.0     0.0   100.00    0.0     0.0     0.0     0.0    0.0     % magnetite (mgt)

                   52.40    0.0     2.45    7.99   16.19   20.81    0.16    0.0    0.0     % diopside (dps)
                   52.98    0.0     0.89   13.90   12.82   18.60    0.81    0.0    0.0     % Mg-augite (mau)
                   52.15    0.0     0.44   25.27    3.68   15.29    3.17    0.0    0.0     % Fe-augite (fau)

                   51.51    0.0     3.29   20.35   22.08    2.77    0.0     0.0    0.0     % hypersthene (hyp)
                   49.06    0.0     0.31   40.11    9.63    0.89    0.0     0.0    0.0     % ferrosillite (fsl)

                    0.0    51.59    0.0    48.41    0.0     0.0     0.0     0.0    0.0     % ilmenite (ilm)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  1  0  0  0  0  0  0  0  0  0  0  0  0    % feldspar (fsp) 
               0  0  0  1  1  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  0  0  0  1  1  0  0  0  0  0  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  0  0  0  1  1  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                ant     alb     san     for     fay     tms     mgt     dps     mau     fau     hyp     fsl     ilm     qtz     wat
cal.cmp_mem = [91.54    8.46       0       0       0       0       0       0       0       0       0       0       0       0       0
               69.85    6.87       0   17.33       0    5.96       0       0       0       0       0       0       0       0       0
               36.45    0.99       0   19.20    7.02    5.63    5.44    2.99   22.28       0       0       0       0       0       0
               23.61    5.52    5.10       0    2.73    2.33    1.95    0.73    1.96   24.82   31.25       0       0       0       0
                6.01   72.49    1.76       0       0       0    0.03       0    2.59    5.07    1.25    9.66    1.13       0       0
                   0   13.31   38.03       0       0       0       0       0       0       0       0    3.27       0   45.39       0
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
cal.T0 =  [1530  1170  1130  1099  972  780];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8 5 4.4 4 3 2.5];

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [28.5  2.1  2.2  8.0  15.4  8.5];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.00;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              ant  alb  san  for  fay  ulv  mgt  dps  mau  fau  hyp  fsl  ilm  qtz  wat
cal.rhox0   = [2690,2590,2530,3400,3930,4120,4970,3270,3400,3470,3390,3660,4700,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
cal.etax0 = 1e18;                   % solid reference viscosity [Pas]
cal.etaf0 = 1e-1;                   % vapour reference viscosity [Pas]
cal.Eax   = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA  =  [ 0.65, 0.25, 0.35; ...  % permission slopes
             0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
             0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB  =  [ 0.55, 0.18, 0.27; ...  % permission step locations
             0.64,0.012,0.348; ...  % each row sums to 1
             0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC  =  [[0.30, 0.30, 0.40]*0.7; ... % permission step widths
            [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
            [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

% % specify mixture viscosity parameters (Costa et al., 2009)
% cal.Bphi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient bubbles
% cal.Bchi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient crystals
% cal.chi_pck = 0.60;                 % rheologically critical crystal fraction
% cal.gamma   = 2.50;                 % step-function steepness coefficient
% cal.delta   = 27;                   % solid viscosity melt-weakening slope
% cal.xi      = 4.5e-4;               % solid viscosity level
% cal.etaf0   = 0.1;                  % fluid viscosity constant
% 
% % specify segregation coefficient parameters
% cal.bm      = 50;                   % melt permeability geometric factor (k0 = dx^2/bm)
% cal.cm      = 0.001;                % melt percolation threshold
% cal.nm      = 3;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
% cal.mm      = 2;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
% cal.bf      = 50;                   % fluid permeability geometric factor (k0 = dx^2/bm)
% cal.cf      = 0.05;                 % fluid percolation threshold
% cal.nf      = 4;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
% cal.mf      = 2;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
