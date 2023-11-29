% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 12;
cal.nmsy   = 6;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','hyp','dps','pig','ant','alb','ams','ulv','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','fsp','oxs','qtz'};
cal.cmpStr = {'dun','pxn','cp3','cp4','cp5','cp6','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   TiO2   Al2O3     FeO     MgO     CaO    Na2O    H2O
cal.mem_oxd    = [ 41.87    0.0     0.0     1.85   56.28    0.00    0.0    0.0     % forsterite (for)
                   34.37    0.0     0.0    44.95   20.68    0.00    0.0    0.0     % fayalite (fay)

                   54.99    0.0     5.30    5.76   32.55    1.40    0.0    0.0     % enstatite (ens)
                   51.17    0.0     8.08   10.80   27.53    2.42    0.0    0.0     % hypersthene (hyp)

                   50.32    0.33    7.36    7.61   18.66   15.53    0.19   0.0     % diopside (dps)
                   48.12    1.36    2.63   27.46    9.66   10.62    0.15   0.0     % pigeonite (mau)

                   47.56    0.0    33.65    0.0     0.0    16.97    1.82   0.0     % anorthite (ant)
                   55.38    0.0    28.56    0.0     0.0    11.98    4.08   0.0     % albite (alb)

                    0.0    18.19    6.21   72.22    3.38    0.0     0.0    0.0     % Al-Mg-bearing spinel (ams)
                    0.0    21.40    2.74   74.97    0.89    0.0     0.0    0.0     % ulvospinel (ulv)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  1  1  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   for     fay     ens     hyp     dps     pig     ant     alb     ams     ulv     qtz     wat
cal.cmp_mem =   [100.00       0       0       0       0       0       0       0       0       0       0       0
                  11.19   17.15   71.66       0       0       0       0       0       0       0       0       0
                   9.09   11.67       0    0.01   31.43       0   47.79       0       0       0       0       0
                      0       0       0       0       0   40.64   18.02   11.94       0   20.23    9.18       0
                      0       0       0       0       0   40.64   18.02   11.94       0   20.23    9.18       0
                      0       0       0       0       0   40.64   18.02   11.94       0   20.23    9.18       0
                      0       0       0       0       0       0       0       0       0       0       0  100.00];
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
cal.T0 =  [1890    1537    1470    1263    1201    1069];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  [9.7000    7.5000    4.7000    3.6000    2.7000    2.5000];

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [32.0198   30.2777   25.1996   20.0875   14.5849    6.6382];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;0.9;3.00;10.0;1.1].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  ens  hyp  dps  pig  ant  alb  ams  ulv  qtz  wat
cal.rhox0   = [3270,4390,3270,3600,3250,3450,2690,2570,4750,4950,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
cal.etax0 = 1e18;                   % solid reference viscosity [Pas]
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
