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
cal.memStr = {'ant','alb','san','for','fay','ulv','mgt','ilm','dps','aug','pig','hyp','fsl','qtz','wat'};
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
                   67.0400         0   18.2800         0         0         0    4.0600   10.6200         0

                   39.5600         0         0   17.0300   43.4100         0         0         0         0
                   37.5500         0         0   27.6900   34.7600         0         0         0         0

                         0   10.1400    6.0400   64.8200   19.0000         0         0         0         0
                         0   14.1100         0   81.2600    4.6300         0         0         0         0
                         0   42.9900         0   57.0100         0         0         0         0         0

                   52.1400         0    2.8600    7.1200   16.5400   21.2000    0.1400         0         0
                   53.7600         0    0.7300   11.7700   13.7600   18.8500    1.1300         0         0
                   53.4000         0    0.2800   17.7700    8.5800   17.3800    2.5900         0         0
                   
                   52.4600         0    4.3200   13.4800   27.2800    2.4600         0         0         0
                   51.6800         0    0.0400   30.6700   16.6900    0.9200         0         0         0

                  100.0000         0         0         0         0         0         0         0         0
                         0         0         0         0         0         0         0         0  100.0000 ];
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  1  0  0  0  0  0  0  0  0  0  0  0  0    % feldspar (fsp) 
               0  0  0  1  1  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  0  0  0  1  1  1  0  0  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  1  1  1  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  0  0  0  0  1  1  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                 ant     alb     san     for     fay     ulv     ilm     dps     aug     pig     hyp     fsl     qtz     wat
cal.cmp_mem = [90.4000    9.6000         0         0         0         0         0         0         0         0         0         0         0         0         0
               61.0000    9.0000         0   17.9000    2.5000    5.3000    4.3000         0         0         0         0         0         0         0         0
               34.0000    4.5000         0         0    0.5000    3.5000    5.1000    2.7000   16.2000    6.2000         0   27.2000         0         0         0
               13.7000   55.4000    0.5000         0         0         0    0.7000    1.5000         0    7.5000    4.9000    1.5000   14.2000         0         0
                0.5000   46.0000   48.0000         0         0         0         0    1.4000         0         0    2.9000         0    1.2000         0         0
                3.3000    1.3000   43.7000         0         0         0         0         0    1.7000    0.6000    2.5000         0    1.1000   45.9000         0
                     0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000];
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
cal.c0     = [0.1200    0.1210    0.3240    0.3130    0.0680    0.0530    0.0300];
cal.c1     = [0.0010    0.0010    0.0020    0.0430    0.3690    0.5850    0.0200];

% set pure component melting points T_m^i at P=0
cal.T0 =  [ 1500   1170   1130   1070    940    815];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [ 5.40   5.10   4.50   3.00   1.90   1.00];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [ 5.50   5.30   4.50   3.00   2.80   3.00];

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [28.00   3.00   3.00  10.00  10.00   5.00];

% specify melting point dependence on H2O
cal.dTH2O   = [1100  1400  1500  1600  1800  2100];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O    = 0.75;                                  % solidus shift from water content exponent

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               ant  alb  san  for  fay  ulv  mgt  ilm  dps  aug  pig  hyp  fsl  qtz  wat
cal.rhox0   = [2530,2190,2220,3340,3470,4330,4540,4790,3250,3320,3360,3300,3500,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             ant  alb  san  for  fay  ulv  mgt  ilm  dps  aug  pig  hyp  fsl  qtz  wat
cal.etax0 = [1e17,1e17,1e17,1e18,1e18,1e16,1e16,1e16,1e19,1e19,1e19,1e19,1e19,1e19,1e0]; % mem ref viscosities [Pas]
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
