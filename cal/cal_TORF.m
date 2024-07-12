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
cal.cmpStr = {'dun','ogb','fbs','tra','rhy','vol'}; % need changes? 

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

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
cal.cmp_mem = [100.00   0.00   0.00   0.00   0.00   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.00 
                33.75  61.95   4.30   0.00   0.00   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00   0.00 
                 1.98   0.01  54.79   8.61  32.66   0.00  0.00  1.95  0.00  0.00  0.00  0.00  0.00   0.00 
                 0.00   7.51   4.30  25.71  28.58  19.84  0.00  2.98  5.16  0.00  5.93  0.00  0.00   0.00 
                 0.00   0.00   1.22   6.01   0.00  78.73  3.53  0.00  0.01  0.77  0.43  9.30  0.00   0.00 
                 0.00   0.00   0.00   6.89   0.00   0.00 47.02  0.00  0.00  0.01  0.00  0.01 46.08   0.00 
                 0.00   0.00   0.00   0.00   0.00   0.00  0.00  0.00  0.00  0.00  0.00  0.00  0.00 100.00]; 


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
cal.T0  = [1880  1276  1136  1108  1025  826];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [3.30  2.82  1.96  3.42  2.63  1.04];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8.64  3.16  2.63  3.45  2.80  2.33];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [32.40  32.70  5.80  9.60  15.10 14.60];

% initial and final composition used in calibration
cal.c0 = [0.044   0.047   0.427   0.352   0.119    0.010    0.005];
cal.c1 = [0.001   0.001   0.001   0.022   0.269    0.706    0.005];

% specify melting point dependence on H2O
%cal.dTH2O   = 1600;                
cal.dTH2O   = [1300    1450    1500    1650    1900    2400];   % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                                                                % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%             'for','fay','dps','aug','ant','alb','san','ulv','mgt','ilm','hyp','fsl','qtz','wat'
cal.rhox0   = [3150,4055,3200,3475,2680,2580,2555,3925,4730,4875,3240,3665,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
cal.etax0 = 1e18;                   % solid reference viscosity [Pas]
cal.etaf0 = 1e-1;                   % vapour reference viscosity [Pas]
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

% specify mixture viscosity parameters (Costa et al., 2009)
cal.Bphi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient bubbles
cal.Bchi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient crystals
cal.chi_pck = 0.60;                 % rheologically critical crystal fraction
cal.gamma   = 2.50;                 % step-function steepness coefficient
cal.delta   = 27;                   % solid viscosity melt-weakening slope
cal.xi      = 4.5e-4;               % solid viscosity level
cal.etaf0   = 0.1;                  % fluid viscosity constant

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
