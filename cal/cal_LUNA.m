% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 10;
cal.nmsy   = 5;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','dps','pig','ant','alb','ulv','qtz','wat'};
cal.msyStr = {'olv','pxn','fsp','oxs','qtz'};
cal.cmpStr = {'dun','opx','web','gbr','bas','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7       9]; % oxdie indices for viscosity, density functions

%           plg  olv  cpx  qtz
cal.imsy = [ 3    1    2    5];  % mineral system indices for plotting basalt tetrahedron

% oxide composition of mineral end-members
%                   SiO2        TiO2       Al2O3     FeO     MgO          CaO       Na2O    H2O
cal.mem_oxd    = [ 41.8500         0         0    1.8300   56.3200         0         0         0   % forsterite (for)
                   34.3900         0         0   45.0000   20.6100         0         0         0   % fayalitic olivine (fay)

                   53.7000    0.0100    5.7400    6.5600   33.9900         0         0         0   % enstatite
                   49.5700    0.1200    9.6100    4.7000   18.5600   17.4400         0         0   % diopside (dps)
                   48.0400    1.3700    2.6300   26.1600    8.2000   13.6000         0         0   % pigeonite (pig)

                   46.4600         0   33.8500    0.7300         0   17.4200    1.5400         0   % anorthite (ant)
                   54.6700         0   27.7500    2.0300         0   11.3600    4.1900         0   % albitic plg (alb)

                         0   19.4900    4.7900   73.3600    2.3600         0         0         0   % ulvospinel (ulv)

                  100.0000         0         0         0         0         0         0         0   % quartz (qtz)
                         0         0         0         0         0         0         0  100.0000 ];% water (wat)

cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0    % pyroxene (pxn)
               0  0  0  0  0  1  1  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                  for       fay       ens       dps       pig       ant       alb       ulv       qtz       wat
cal.cmp_mem = [85.0000   15.0000         0         0         0         0         0         0         0         0   % dunite
                1.3000    6.6000   92.1000         0         0         0         0         0         0         0   % orthopyroxenite
                9.4000   17.2000   19.8000   51.0000         0    2.6000         0         0         0         0   % olivine websterite
                     0   23.1000         0   36.1000    0.2000   34.3000    6.3000    0.1000         0         0   % olivine gabbro
                     0         0         0         0   44.9000         0   30.1000   19.9000    5.1000         0   % basal
                     0         0         0         0         0         0         0         0         0  100.0000];
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
cal.T0 =  [1800  1534  1165  1074  1000];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [6.2  4.9  3.3  2.8  2.5];

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  [6.0  5.0  2.7  2.6  2.4];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [21.8  20.6  20.1  5.2  14.1];

% initial composition used in calibration
cal.c0 = [0.386  0.175  0.158  0.220  0.061  0.000];

% specify melting point dependence on H2O
cal.dTH2O   = [933  1084  1344  1527  1680];  % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                             % solidus shift from water content [K/wt^pH2O]

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx =  350;
cal.Dsf =  450;

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.00;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  ens  dps  pig  ant  alb   ulv  qtz  wat
cal.rhox0   = [3270,4390,3470,3250,3450,2690,2570,4850,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
%             for  fay  ens  dps  pig  ant  alb  ulv  qtz  wat
cal.etax0 = [1e18,1e18,1e19,1e20,1e20,1e17,1e17,1e16,1e17,1e0]; % mem ref viscosities [Pas]
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
cal.tol     = 1e-12;
cal.alpha   = 0.75;
