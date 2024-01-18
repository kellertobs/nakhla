% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 13;
cal.nmsy   = 6;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','hyp','tdp','dps','pig','ant','alb','ams','ulv','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','fsp','oxs','qtz'};
cal.cmpStr = {'dun','hrz','opt','gbr','bas','fbs','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   TiO2   Al2O3     FeO     MgO     CaO    Na2O    H2O
cal.mem_oxd    = [ 42.03    0.0     0.0     0.45   57.35    0.17    0.0    0.0     % forsterite (for)
                   34.21    0.0     0.0    44.59   20.61    0.59    0.0    0.0     % fayalitic olivine (fay)

                   55.01    0.0     5.28    5.77   32.54    1.40    0.0    0.0     % enstatite (ens)
                   51.18    0.0     8.08   10.79   27.53    2.42    0.0    0.0     % hypersthene (hyp)

                   49.60    0.28    8.81    7.77   17.79   15.57    0.18   0.0     % tschermakitic diopside (tds)
                   50.94    0.53    4.93    9.95   18.51   14.94    0.20   0.0     % diopside (dps)
                   47.81    1.42    2.65   28.85    8.98   10.15    0.14   0.0     % pigeonite (pig)

                   47.60    0.0    33.60    0.0     0.0    17.00    1.80   0.0     % anorthite (ant)
                   55.34    0.0    28.62    0.0     0.0    11.93    4.12   0.0     % albitic plagioclase (alb)

                    0.0    18.04    6.20   72.39    3.37    0.0     0.0    0.0     % Al-Mg-bearing spinel (ams)
                    0.0    21.62    2.72   74.79    0.87    0.0     0.0    0.0     % ulvospinel (ulv)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  1  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                 for     fay     ens     hyp     tdp     dps     pig     ant     alb     ams     ulv     qtz     wat
cal.cmp_mem = [100.00       0       0       0       0       0       0       0       0       0       0       0       0   % dunite (dun)
                22.87   37.68   39.45       0       0       0       0       0       0       0       0       0       0   % harzburgite (hrz)
                 1.25    1.00   65.40   31.35    1.01       0       0       0       0       0       0       0       0   % orthopyroxenite (opt)
                    0    6.82       0    1.01   45.13    1.96       0   45.08       0       0       0       0       0   % gabbro (gbr)
                    0       0       0       0       0   32.59   16.36   41.23    2.97    6.85       0       0       0   % basalt (bas)
                    0       0       0       0       0       0   47.47       0   30.22       0   17.57    4.74       0   % ferro-basalt (fbs)
                    0       0       0       0       0       0       0       0       0       0       0       0  100.00]; % volatiles (vol)
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
cal.T0 =  [1890  1540  1475  1177  1073  1014];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  [8.37  6.82  5.27  2.38  2.00  1.85];

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [26.2  55.4  15.2  9.6  6.3  9.6];

% initial composition used in calibration
cal.c0 = [0.30  0.33  0.09  0.21  0.04  0.03  0.00];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.00;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  ens  hyp  tdp  dps  pig  ant  alb  ams  ulv  qtz  wat
cal.rhox0   = [3270,4390,3270,3600,3200,3250,3450,2690,2570,4750,4950,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
cal.etax0 = 1e19;                   % solid reference viscosity [Pas]
cal.etaf0 = 1e-1;                   % vapour reference viscosity [Pas]
cal.Eax   = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA  =  [ 0.65, 0.25, 0.35; ...  % permission slopes
             0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
             0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB  =  [ 0.50, 0.20, 0.30; ...  % permission step locations
             0.64,0.012,0.348; ...  % each row sums to 1
             0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC  =  [[0.30, 0.30, 0.40]*0.7; ... % permission step widths
            [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
            [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions
