% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 14;
cal.nmsy   = 6;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ant','alb','san','dps','aug','ulv','mgt','ilm','hyp','fsl','qtz','wat'};
cal.msyStr = {'olv','fsp','cxp','oxs','opx','qtz'};
cal.cmpStr = {'dun','tro','gbr','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                 SiO2      TiO2     Al2O3       FeO       MgO       CaO      Na2O       K2O       H2O
cal.mem_oxd = [ 41.2800         0         0    8.0100   50.7100         0         0         0         0
                31.7400         0         0   58.6200    9.6400         0         0         0         0
                44.7800         0   35.4600         0         0   18.8000    0.9600         0         0
                67.6100         0   20.2000         0         0    0.8900   11.2700    0.0300         0
                67.5800         0   19.4600         0         0    0.2800    6.2700    6.4100         0
                53.6400         0    2.5500    4.7300   19.5900   19.4500    0.0400         0         0
                51.7500         0    0.3700   24.5700    5.1000   15.6200    2.5900         0         0
                      0   38.9300    2.6100   28.4400   30.0100         0         0         0         0
                      0   12.7200    1.0200   86.2700         0         0         0         0         0
                      0   53.0300         0   46.9700         0         0         0         0         0
                51.1700         0    2.7200   23.3700   20.1300    2.6100         0         0         0
                49.0700         0    0.3500   39.5600   10.0700    0.9500         0         0         0
               100.0000         0         0         0         0         0         0         0         0
                      0         0         0         0         0         0         0         0  100.0000];
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  1  1  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  1  1  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%              for   fay   ant   alb   san   dps   aug   pig   ulv   mgt   ilm   hyp   fsl   qtz   wat
cal.cmp_mem = [92.9000    7.1000         0         0         0         0         0         0         0         0         0         0         0         0
               43.8000    4.4000   51.9000         0         0         0         0         0         0         0         0         0         0         0
                0.1000    0.1000   34.0000   10.4000         0   52.3000         0    3.1000         0         0         0         0         0         0
                     0    7.4000   20.0000   24.8000    0.7000   13.4000   25.5000         0    3.6000    0.9000    3.7000         0         0         0
                     0         0    5.4000   64.9000    6.9000         0   11.0000         0    0.1000    1.3000         0   10.5000         0         0
                     0         0    2.1000         0   44.1000         0    3.5000         0         0         0         0    0.3000   50.0000         0
                     0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000];
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

% primary and evolved end-member compositions used in calibration
cal.c0     = [0.0090    0.2010    0.3080    0.3670    0.1150         0    0.0050];
cal.c1     = [     0         0         0         0    0.3160    0.6840    0.0240];

cal.c0_oxd = [49.20  1.01  15.11  9.58  11.55  11.31  2.12  0.12  0.30];
cal.c1_oxd = [75.26  0.24  11.51  3.40   0.77   2.01  4.61  2.20  2.40];

% set pure component melting points T_m^i at P=0
cal.T0  = [1600  1193  1158  1078  980  819];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [6.3  5.0  3.8  2.4  1.8  0.8];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [6.3  5.2  3.8  2.6  2.3  2.1];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [12.3  2.2  3.0  7.2  9.9  5.3];

% specify melting point dependence on H2O
cal.dTH2O   = [1050  1410  1460  1580  1730  2060];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O    = 0.75;                                  % solidus shift from water content exponent

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  hyp  fsl  qtz  wat
cal.rhox0   = [3200,4050,2680,2600,2570,3220,3460,3930,4760,4720,3310,3660,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  hyp  fsl  qtz  wat
cal.etax0   = [1e19,1e19,1e17,1e17,1e17,1e20,1e20,1e16,1e16,1e16,1e20,1e20,1e17,1e0]; % mem ref viscosities [Pas]
cal.etaf0   = 0.1;                  % fluid viscosity constant [Pas]
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
