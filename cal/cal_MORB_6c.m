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
cal.msyStr = {'olv','fsp','cxp','spn','opx','qtz'};
cal.cmpStr = {'dun','tro','gbr','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                 SiO2      TiO2     Al2O3       FeO       MgO       CaO      Na2O       K2O       H2O
cal.mem_oxd = [41.3200         0         0    7.8300   50.8500         0         0         0         0
               31.0700         0         0   62.2400    6.6900         0         0         0         0
               44.1200         0   35.9400         0         0   19.3600    0.5800         0         0
               65.8000         0   21.4100         0         0    2.3100   10.4200    0.0600         0
               67.8200         0   19.5900         0         0    0.3400    5.1400    7.1100         0
               53.3300         0    3.2600    4.4300   20.6300   18.3500         0         0         0
               51.5500         0    0.3400   25.0200    4.8900   15.8300    2.3700         0         0
                     0   38.3500    2.9600   31.0500   27.6400         0         0         0         0
                     0   10.2900    1.2700   88.4400         0         0         0         0         0
                     0   52.4900         0   47.5100         0         0         0         0         0
               52.2000         0    4.0600   15.4800   25.1900    3.0700         0         0         0
               48.4400         0    0.5000   41.2700    8.6600    1.1300         0         0         0
              100.0000         0         0         0         0         0         0         0         0
                     0         0         0         0         0         0         0         0  100.0000];  % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  1  1  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0    % spinel (spn)
               0  0  0  0  0  0  0  0  0  0  1  1  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
cal.cmp_mem = [100.0000         0         0         0         0         0         0         0         0         0         0         0         0         0
                17.3000   64.7000   18.0000         0         0         0         0         0         0         0         0         0         0         0
                      0         0   42.4000   12.4000         0   41.2000         0    0.2000         0         0    3.8000         0         0         0
                      0    6.4000   21.5000   20.6000         0   10.2000   27.1000    2.4000    5.6000         0    3.7000    2.5000         0         0
                      0         0         0   74.7000    3.6000         0   10.4000         0         0    0.7000    0.8000    9.8000         0         0
                      0         0         0         0   45.6000         0    4.8000         0         0         0         0    0.1000   49.5000         0
                      0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000]; % volatiles (vol)
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
cal.T0  = [1875  1199  1159  1092  1003  837];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [7.4200    2.9100    3.1400    2.8100    2.4300    1.0000];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [7.0400    3.3400    3.9000    3.1700    2.5300    2.5300];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [37.8000    9.1000    4.3000    8.7000   12.8000   10.7000];

% initial and final composition used in calibration
cal.c0 = [0.0900    0.0550    0.4710    0.2860    0.0950    0.0030    0.0030];
cal.c1 = [0.0020    0.0010    0.0010    0.0080    0.2440    0.7440    0.0300];

% specify melting point dependence on H2O
cal.dTH2O   = [1300        1450        1500        1650        1900        2400];                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  ens  fsl  qtz  wat
cal.rhox0   = [3210,4060,2680,2600,2560,3210,3460,3930,4760,4720,3300,3660,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  ens  fsl  qtz  wat
cal.etax0 = [1e19,1e18,1e17,1e17,1e17,1e20,1e20,1e16,1e16,1e16,1e20,1e20,1e17,1e0]; % mem ref viscosities [Pas]
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


% specify segregation coefficient parameters
cal.bm      = 50;                   % melt permeability geometric factor (k0 = dx^2/bm)
cal.cm      = 0.001;                % melt percolation threshold
cal.nm      = 3;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.mm      = 2;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.bf      = 50;                   % fluid permeability geometric factor (k0 = dx^2/bm)
cal.cf      = 0.05;                 % fluid percolation threshold
cal.nf      = 4;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
cal.mf      = 2;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)

% convergence tolerance
cal.tol     = 1e-9;
