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
cal.memStr = {'for','fay','ens','hyp','tds','dps','pig','ant','alb','ams','ulv','qtz','wat'};
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
cal.mem_oxd    = [ 41.8500         0         0    1.8300   56.3200         0         0         0   % forsterite (for)
                   34.3900         0         0   45.0000   20.6100         0         0         0   % fayalitic olivine (fay)

                   55.0100         0    5.2800    5.7700   32.5400    1.4000         0         0   % enstatite (ens)
                   51.1800         0    8.0800   10.7900   27.5300    2.4200         0         0   % hypersthene (hyp)

                   49.7000    0.3000    8.8600    7.4500   17.7200   15.9700         0         0   % tchermakitic diopside (tds)
                   50.8400    0.4100    5.8400    9.2000   18.8700   14.8400         0         0   % diopside (dps)
                   47.2200    1.8000    0.8700   35.4700    6.0000    8.6400         0         0   % pigeonite (pig)

                   46.9700         0   33.4700    0.8100         0   17.0400    1.7100         0   % anorthite (ant)
                   54.1700         0   28.1300    1.9500         0   11.7300    4.0200         0   % albitic plg (alb)

                         0   18.0400    6.2000   72.3900    3.3700         0         0         0   % Al-Mg-bearing ulvospinel (ams)
                         0   21.6300    2.7200   74.7800    0.8700         0         0         0   % ulvospinel (ulv)

                  100.0000         0         0         0         0         0         0         0   % quartz (qtz)
                         0         0         0         0         0         0         0  100.0000 ];% water (wat)

cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  1  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                 for     fay     ens     hyp     tds     dps     pig     ant     alb     ams     ulv     qtz     wat
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
cal.B  =  [8.37  6.52  5.09  2.47  2.11  1.95];

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx =  350;
cal.Dsf =  450;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [26.2  55.4  15.2  9.6  6.3  9.6];

% initial composition used in calibration
cal.c0 = [0.30  0.33  0.09  0.21  0.04  0.03  0.00];

% specify melting point dependence on H2O
cal.dTH2O   = [1400,1400,1400,1400,1400,1400];  % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                             % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.00;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  ens  hyp  tdp  dps  pig  ant  alb  ams  ulv  qtz  wat
cal.rhox0   = [3270,4390,3270,3600,3200,3250,3450,2690,2570,4750,4950,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
%             for  fay  ens  hyp  tdp  dps  pig  ant  alb  ams  ulv  qtz  wat
cal.etax0 = [1e18,1e18,1e19,1e19,1e19,1e20,1e20,1e17,1e17,1e16,1e16,1e17,1e0]; % mem ref viscosities [Pas]
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

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;
