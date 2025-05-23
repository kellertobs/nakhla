% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 14;
cal.nmsy   = 6;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','dps','aug','ant','alb','san','ulv','mgt','ilm','hyp','fsl','qtz','wat'};
cal.msyStr = {'olv','cxp','fsp','spn','opx','qtz'};
cal.cmpStr = {'dun','ogb','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

%           plg  olv  cpx  qtz
cal.imsy = [ 3    1    2    6];  % mineral system indices for plotting basalt tetrahedron

% oxide composition of mineral end-members
%                 SiO2      TiO2     Al2O3       FeO       MgO       CaO      Na2O       K2O       H2O
cal.mem_oxd = [41.3200         0         0    7.8800   50.8000         0         0         0         0
               31.6400         0         0   59.1500    9.2100         0         0         0         0
               52.6400         0    3.7100    3.3400   20.1500   20.1600         0         0         0
               51.7800         0    0.1200   24.9300    4.8500   15.7400    2.5800         0         0
               43.7700         0   36.2000         0         0   19.6500    0.3800         0         0
               66.0700         0   21.2400         0         0    2.1000   10.5700    0.0200         0
               68.1300         0   19.4700         0         0    0.1700    6.5400    5.6900         0
                     0   35.6800    3.2800   35.9700   25.0700         0         0         0         0
                     0    5.3800    1.8200   92.7200    0.0800         0         0         0         0
                     0   53.2000         0   46.8000         0         0         0         0         0
               51.1900         0    4.0200   19.0900   22.6600    3.0400         0         0         0
               49.1100         0    0.3600   39.3500   10.1500    1.0300         0         0         0
              100.0000         0         0         0         0         0         0         0         0
                     0         0         0         0         0         0         0         0  100.0000];
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  1  1  1  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0    % spinel (spn)
               0  0  0  0  0  0  0  0  0  0  1  1  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)


% mineral end-member composition of melting model components
%cal.cmp_mem = ones(cal.ncmp,cal.nmem); %uncomment when restarting process
cal.cmp_mem =[100.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.00  0.00  0.00   0.00  0.00   0.00
                9.60   1.84  62.84   0.00  25.65   0.00   0.00   0.07  0.00  0.00  0.00   0.00  0.00   0.00
                0.00   8.30   9.53  22.95  27.94  18.98   0.00   3.34  3.78  1.10  4.09   0.00  0.00   0.00
                0.00   0.00   0.15   7.98   0.00  75.26   2.84   0.00  0.01  0.36  1.50  11.90  0.00   0.00
                0.00   0.00   0.00   7.57   0.00   0.00  49.40   0.00  0.00  0.01  0.00   0.01 43.01   0.00
                0.00   0.00   0.00   0.00   0.00   0.00   0.00   0.00  0.00  0.00  0.00   0.00  0.00 100.00]; 


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
cal.T0  = [1732  1167  1110  1034  819];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [5.82  2.67  3.42  3.43  1.00];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [9.52  3.04  2.70  1.91  2.17];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [23.50 4.10  8.00  10.50 14.30];

% initial and final composition used in calibration
cal.c0 = [0.034  0.413  0.412  0.130  0.010  0.005];
cal.c1 = [0.001  0.001  0.020  0.277  0.700  0.023];

% specify melting point dependence on H2O                
cal.dTH2O   = [1300   1500  1650  1900  2200];   % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                              % solidus shift from water content [K/wt^pH2O]

% set entropy gain of crystallisation and evaporation DeltaS [J/K]
cal.Dsx = 350;
cal.Dsf = 500;

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  dps  aug  ant  alb  san  ulv  mgt  ilm  ens  fsl  qtz  wat
cal.rhox0   = [3210,4000,3200,3460,2690,2600,2560,4000,4740,4710,3360,3640,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%             for  fay  dps  aug  ant  alb  san  ulv  mgt  ilm  ens  fsl  qtz  wat
cal.etax0 = [1e19,1e19,1e20,1e20,1e17,1e17,1e17,1e16,1e16,1e16,1e20,1e20,1e17,1e0]; % mem ref viscosities [Pas]
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
