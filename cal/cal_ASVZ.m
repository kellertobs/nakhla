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
cal.memStr = {'ant','alb','san','for','fay','ulv','ilm','dps','aug','pig','hyp','fsl','qtz','wat'};
cal.msyStr = {'fsp','olv','oxs','cxp','opx','qtz'};
cal.cmpStr = {'ano','str','sgn','and','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                     SiO2      TiO2     Al2O3       FeO       MgO       CaO      Na2O       K2O       H2O
cal.mem_oxd    = [ 43.4800         0   36.1900         0         0   19.7100    0.6200         0         0
                   66.2900         0   21.0300         0         0    1.8800   10.8000         0         0
                   67.0400         0   18.2800         0         0         0    6.0600    8.6200         0
                  
                   39.4400         0         0   17.6200   42.9400         0         0         0         0
                   37.6700         0         0   27.0900   35.2400         0         0         0         0
                 
                         0   11.5300    2.8900   73.5200   12.0600         0         0         0         0
                         0   40.6000         0   59.4000         0         0         0         0         0
                 
                   52.1400         0    2.8600    7.1200   16.5400   21.2000    0.1400         0         0
                   53.7600         0    0.7300   11.7700   13.7600   18.8500    1.1300         0         0
                   53.4000         0    0.2800   17.7700    8.5800   17.3800    2.5900         0         0
                 
                   52.4200         0    4.0800   14.4400   26.6900    2.3700         0         0         0
                   51.7200         0    0.2700   29.7200   17.2800    1.0100         0         0         0
                
                  100.0000         0         0         0         0         0         0         0         0
                         0         0         0         0         0         0         0         0  100.0000 ];
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  1  0  0  0  0  0  0  0  0  0  0  0    % feldspar (fsp) 
               0  0  0  1  1  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  0  0  0  1  1  0  0  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  0  0  0  1  1  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                 ant     alb     san     for     fay     ulv     ilm     dps     aug     pig     hyp     fsl     qtz     wat
cal.cmp_mem = [ 92     8     0     0     0     0     0     0     0     0     0     0     0     0
                51    19     0    22     0     7     0     0     0     0     0     0     0     0
                25     2     0     0     2    13     5    21     2     0    30     0     0     0
                 9    56     0     0     0     0     2     0    10     6     2    15     0     0
                 1    30    64     0     0     0     1     0     0     3     0     1     0     0
                 3     0    53     0     0     0     0     2     0     3     0     1    38     0
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

% set pure component melting points T_m^i at P=0
cal.T0 =  [1553  1171  1134  1092  990  745];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8.0  5.0  4.5  4.0  3.0  2.5];

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [21.0  5.1  2.0  9.5  12.4  7.4];

% specify melting point dependence on H2O
cal.dTH2O   = [1100  1400  1500  1600  1800  2100];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O    = 0.75;                                  % solidus shift from water content exponent

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               ant  alb  san  for  fay  ulv  ilm  dps  aug  pig  hyp  fsl  qtz  wat
cal.rhox0   = [2530,2190,2220,3340,3470,4330,4790,3250,3320,3360,3300,3500,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             ant  alb  san  for  fay  ulv  ilm  dps  aug  pig  hyp  fsl  qtz  wat
cal.etax0 = [1e17,1e17,1e17,1e18,1e18,1e16,1e16,1e19,1e19,1e19,1e19,1e19,1e0]; % mem ref viscosities [Pas]
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
