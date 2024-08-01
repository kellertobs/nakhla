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
cal.mem_oxd = [41.3600         0         0    7.5900   51.0500         0         0         0         0
               31.1900         0         0   61.5900    7.2200         0         0         0         0
               44.1900         0   35.9100         0         0   19.3100    0.5900         0         0
               64.6800         0   22.1600         0         0    2.6900   10.4000    0.0700         0
               68.0700         0   19.5400         0         0    0.2300    5.5900    6.5700         0
               53.3300         0    3.1800    4.9000   20.3500   18.2400         0         0         0
               51.5800         0    0.3200   24.9400    4.9500   15.8600    2.3500         0         0
                     0   37.8600    2.8500   33.0300   26.2600         0         0         0         0
                     0    9.9500    1.2500   88.8000         0         0         0         0         0
                     0   51.7500         0   48.2500         0         0         0         0         0
               52.0800         0    3.9800   16.1200   24.7600    3.0600         0         0         0
               48.5100         0    0.5000   41.0100    8.8300    1.1500         0         0         0
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
%              for   fay   ant   alb   san   dps   aug   ulv   mgt   ilm   hyp   fsl   qtz   wat
cal.cmp_mem = [ 97     3     0     0     0     0     0     0     0     0     0     0     0     0
                40    25    32     0     0     3     0     0     0     0     0     0     0     0
                 0     0    42     7     0    37     3     3     0     0     8     0     0     0
                 0    11     4    35     1     2    33   1.5   7.3   1.2     1     3     0     0
                 0     0     0    71     4     1     9     0   0.2   0.8     2    12     0     0
                 0     0     0     0    45     0     6     0     0   0.2     0   0.8    48     0
                 0     0     0     0     0     0     0     0     0     0     0     0     0   100];
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
cal.c0     = [0.0490    0.1360    0.5200    0.1680    0.1000    0.0270    0.0030];
cal.c1     = [0.0010    0.0010    0.0010    0.0210    0.3240    0.6520    0.0240];

cal.c0_oxd = [49.20  1.01  15.11  9.58  11.55  11.31  2.12  0.12  0.30];
cal.c1_oxd = [75.26  0.24  11.51  3.40   0.77   2.01  4.61  2.20  2.40];

% set pure component melting points T_m^i at P=0
cal.T0  = [1600  1240  1150  1090  990  835];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [6.0  5.5  3.5  2.6  1.8  1.0];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [7.0  6.0  3.3  3.0  2.7  2.6];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [18  7  5  7  10  12];

% specify melting point dependence on H2O
cal.dTH2O   = [1100  1400  1500  1600  1800  2100];  % solidus shift from water content prefactor [K/wt^pH2O]
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
cal.etax0 = [1e19,1e19,1e17,1e17,1e17,1e20,1e20,1e16,1e16,1e16,1e20,1e20,1e17,1e0]; % mem ref viscosities [Pas]
cal.etaf0 = 0.1;                  % fluid viscosity constant [Pas]
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

% convergence tolerance
cal.tol     = 1e-9;
