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
cal.cmpStr = {'ano','snr','gbn','and','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2    TiO2   Al2O3     FeO     MgO     CaO    Na2O     K2O     H2O
cal.mem_oxd    = [ 43.72       0   36.05       0       0   19.54    0.69       0       0    % anorthite (ant)
                   66.90       0   20.65       0       0    1.42   11.03       0       0    % albite (alb)
                   67.36       0   18.25       0       0       0    5.90    8.49       0    % sanidine (san)

                   51.81       0    5.87   12.55   28.54    1.23       0       0       0    % enstatite (ens)
                   52.20       0    3.76   15.95   25.35    2.74       0       0       0    % hypersthene (hyp)
                   50.05       0       0   37.42   11.69    0.84       0       0       0    % ferrosillite (fsl)

                       0   20.74    4.75   54.65   19.86       0       0       0       0    % ulvospinel (ulv)
                       0    9.40    1.67   83.83    5.10       0       0       0       0    % magnetite (mgt)
                       0   46.92       0   53.08       0       0       0       0       0    % ilmenite (ilm)

                   52.32       0    2.72    7.58   16.57   20.57    0.24       0       0    % diopside (dps)
                   53.36       0    0.80   13.35   12.86   18.53    1.10       0       0    % augite (aug)
                   52.49       0    0.26   23.03    5.06   15.99    3.17       0       0    % pigeonite (pig)

                  100.00       0       0       0       0       0       0       0       0    % quartz (qtz)
                       0       0       0       0       0       0       0       0  100.00 ]; % water (wat)
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
cal.cmp_mem = [89.2   10.8      0      0      0      0      0      0      0      0      0      0      0      0   % anorthosite (ano)
               31.3    0.4      0   47.3      0      0   21.0      0      0      0      0      0      0      0   % spinel-norite (snr)
               38.3    0.1      0      0   28.0      0      0   10.9    2.6   20.1      0      0      0      0   % gabbro-norite (gbn)
                9.4   61.4    0.1      0      0   12.2      0      0    1.2      0   11.5    4.2      0      0   % andesite (and)
                1.1   31.3   63.8      0      0    0.9      0      0    0.7      0      0    2.2      0      0   % trachydacite (trd)
                1.9      0   49.8      0      0      0      0      0    1.1    4.1      0    1.6   41.4      0   % rhyolite (rhy)
                  0      0      0      0      0      0      0      0      0      0      0      0      0  100.0]; % volatile (vol)
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
cal.T0  = [1550  1129  1116  1070  922  810];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [8.4  3.3  3.3  2.7  1.8  0.6];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8.4  3.6  3.5  2.8  1.7  2.6];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r   = [35.8  4.4  3.9  13.1  11.5  3.0];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% specify melting point dependence on H2O
cal.dTH2O = [1084  1474  1500  1570  1816  2074];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;                                  % solidus shift from water content exponent

% primary and evolved end-member compositions used in calibration
cal.c0     = [0.1990    0.1110    0.2810    0.2520    0.0830    0.0740    0.0200];
cal.c1     = [0.0010    0.0010    0.0010    0.0300    0.3220    0.6450    0.0250];

cal.c0_oxd = [51.61  1.32  19.26  8.71  5.91  9.33  3.10  0.76  2.00];
cal.c1_oxd = [74.03  0.46  12.90  1.47  0.68  1.48  4.51  4.47  2.50];

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               ant  alb  san  ens  hyp  fsl  ulv  mgt  ilm  dps  aug  pig  qtz  wat
cal.rhox0   = [2690,2580,2540,3280,3350,3580,4220,4600,4770,3250,3330,3420,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             ant  alb  san  ens  hyp  fsl  ulv  mgt  ilm  dps  aug  pig  qtz  wat
cal.etax0  = [1e17,1e17,1e17,1e19,1e19,1e19,1e17,1e17,1e17,1e19,1e19,1e19,1e19,1e0]; % mem ref viscosities [Pas]
cal.etaf0  = 0.1;                  % fluid viscosity constant [Pas]
cal.Eax    = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA     =[ 0.65, 0.25, 0.35; ...  % permission slopes
              0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
              0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB     =[ 0.55, 0.18, 0.27; ...  % permission step locations
              0.64,0.012,0.348; ...  % each row sums to 1
              0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC     =[[0.30, 0.30, 0.40]*0.7; ... % permission step widths
             [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
             [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;
