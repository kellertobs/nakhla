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
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O + K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ant','alb','dps','aug','ulv','mgt','qtz','wat'};
cal.msyStr = {'olv','fsp','cxp','spn','qtz'};
cal.cmpStr = {'dun','gbr','bas','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                SiO2    TiO2   Al2O3     FeO     MgO     CaO    Na2O     H2O
cal.mem_oxd = [41.3700         0         0    7.5700   51.0600         0         0         0
               29.7900         0         0   69.0900    1.1200         0         0         0

               46.5900         0   34.4300         0         0   17.5300    1.4500         0
               68.9400         0   18.6800         0         0         0   12.3800         0

               53.3000    0.0300    2.7400    5.5400   19.4700   18.9200         0         0
               50.0000    1.1800    0.3400   30.8200         0   14.0100    3.6500         0

                     0   40.2200    2.6000   31.6000   25.5800         0         0         0
                     0   14.5700    1.2200   84.2100         0         0         0         0

              100.0000         0         0         0         0         0         0         0
                     0         0         0         0         0         0         0  100.0000]; % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  1  1  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%               for   fay   ant   alb   dps   aug   ulv   mgt   qtz    wat
cal.cmp_mem = [94.3   5.7     0     0     0     0     0     0     0     0   % dun
                2.6   2.6  66.6     0  28.2     0     0     0     0     0   % gbr
                6.5   5.8  30.6   4.2  29.5  19.2   4.2     0     0     0   % bas
                  0   5.7   8.4  73.5     0  11.5   0.1   0.8     0     0   % tra
                  0     0     0  52.7     0     0   2.2   0.1  45.0     0   % rhy
                  0     0     0     0     0     0     0     0     0 100.0]; % vol
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
cal.T0  = [1850  1188  1107  1012  850];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [7.5  4.6  3.5  2.4  2.4];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [7.5  4.6  3.5  2.4  2.4];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [37.0  3.0  4.9  9.7  3.1];

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx  = 350;
cal.Dsf  = 450;

% specify melting point dependence on H2O
cal.dTH2O = [908  1414  1518  1661  1976];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;                           % solidus shift from water content exponent

% primary and evolved end-member compositions used in calibration
cal.c0     = [0.090  0.297  0.414  0.169  0.030  0.005];

cal.c0_oxd = [50.00  1.26  15.08  9.10  10.58  11.31  2.67  0.50];

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  san  dps  aug  ulv  mgt  qtz  wat
cal.rhox0   = [3210,4200,2680,2550,3210,3520,3930,4730, 2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%              for  fay  ant  alb   dps  aug  ulv  mgt  qtz  wat
cal.etax0   = [1e19,1e19,1e17,1e17,1e19,1e19,1e17,1e17,1e19,1e0]; % mem ref viscosities [Pas]
cal.etaf0   = 0.1;                    % fluid viscosity constant [Pas]
cal.Eax     = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA      =[ 0.65, 0.25, 0.35; ...  % permission slopes
               0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
               0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB      =[ 0.55, 0.18, 0.27; ...  % permission step locations
               0.64,0.012,0.348; ...  % each row sums to 1
               0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC      =[[0.30, 0.30, 0.40]*0.7; ... % permission step widths
              [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
              [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;
