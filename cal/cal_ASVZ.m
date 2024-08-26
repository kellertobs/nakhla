% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 14;
cal.nmsy   = 5;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'ant','alb','san','ens','hyp','fsl','ulv','mgt','ilm','dps','aug','pig','qtz','wat'};
cal.msyStr = {'fsp','opx','oxs','cxp','qtz'};
cal.cmpStr = {'ano','gnr','gbr','and','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                     SiO2      TiO2     Al2O3       FeO       MgO       CaO      Na2O       K2O       H2O
cal.mem_oxd    = [ 43.7200         0   36.0500         0         0   19.5400    0.6900         0         0
                   66.9000         0   20.6500         0         0    1.4200   11.0300         0         0
                   67.3600         0   18.2500         0         0         0    5.9000    8.4900         0

                   51.8100         0    5.8700   12.5500   28.5400    1.2300         0         0         0
                   52.2000         0    3.7600   15.9500   25.3500    2.7400         0         0         0
                   50.0500         0         0   37.4200   11.6900    0.8400         0         0         0

                         0   20.7400    4.7500   54.6500   19.8600         0         0         0         0
                         0    9.4000    1.6700   83.8300    5.1000         0         0         0         0
                         0   46.9200         0   53.0800         0         0         0         0         0

                   52.3200         0    2.7200    7.5800   16.5700   20.5700    0.2400         0         0
                   53.3600         0    0.8000   13.3500   12.8600   18.5300    1.1000         0         0
                   52.4900         0    0.2600   23.0300    5.0600   15.9900    3.1700         0         0

                  100.0000         0         0         0         0         0         0         0         0
                         0         0         0         0         0         0         0         0  100.0000 ];
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
%              ant alb san ens hyp fsl ulv mgt ilm dps aug pig qtz wat
cal.msy_mem = [ 1   1   1   0   0   0   0   0   0   0   0   0   0   0    % feldspar (fsp) 
                0   0   0   1   1   1   0   0   0   0   0   0   0   0    % orthopyroxene (opx)
                0   0   0   0   0   0   1   1   1   0   0   0   0   0    % oxides (oxs)
                0   0   0   0   0   0   0   0   0   1   1   1   0   0    % clinopyroxene (cpx)
                0   0   0   0   0   0   0   0   0   0   0   0   1   0];  % quartz (qtz)

% mineral end-member composition of melting model components
%               ant    alb    san    ens    hyp    fsl    ulv    mgt    ilm    dps    aug    pig    qtz    wat
cal.cmp_mem = [89.2000   10.8000         0         0         0         0         0         0         0         0         0         0         0         0
   31.3000    0.4000         0   47.3000         0         0   21.0000         0         0         0         0         0         0         0
   38.3000    0.1000         0         0   28.0000         0         0   10.9000    2.6000   20.1000         0         0         0         0
    9.4000   61.4000    0.1000         0         0   12.2000         0         0    1.2000         0   11.5000    4.2000         0         0
    1.1000   31.3000   63.8000         0         0    0.9000         0         0    0.7000         0         0    2.2000         0         0
    1.9000         0   49.8000         0         0         0         0         0    1.1000    4.1000         0    1.6000   41.4000         0
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
cal.c0     = [0.1990    0.1110    0.2810    0.2520    0.0830    0.0740    0.0200];
cal.c1     = [0.0010    0.0010    0.0010    0.0300    0.3220    0.6450    0.0250];

% set pure component melting points T_m^i at P=0
cal.T0  = [1550        1129        1116        1070         922         810];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [8.4000    3.3000    3.3000    2.7000    1.8000    0.6000];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8.4000    3.6000    3.5000    2.8000    1.7000    2.6000];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r   = [35.8000    4.4000    3.9000   13.1000   11.5000    3.0000];

% specify melting point dependence on H2O
cal.dTH2O   = [1084        1474        1500        1570        1816        2074];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O    = 0.75;                                  % solidus shift from water content exponent

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               ant  alb  san  ens  hyp  fsl  ulv  mgt  ilm  dps  aug  pig  qtz  wat
cal.rhox0   = [2690,2580,2540,3280,3350,3580,4220,4600,4770,3250,3330,3420,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             ant  alb  san  ens  hyp  fsl  ulv  mgt  ilm  dps  aug  pig  qtz  wat
cal.etax0 = [1e17,1e17,1e17,1e18,1e18,1e18,1e16,1e16,1e16,1e19,1e19,1e19,1e18,1e0]; % mem ref viscosities [Pas]
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
cal.alpha   = 0.5;
