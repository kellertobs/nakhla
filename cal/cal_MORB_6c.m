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
                 0.0200   25.4000   54.5800         0         0   19.9900         0         0         0         0         0         0         0         0
                 0.0200    0.0200   28.5500   13.9000         0   29.2400   16.0700    3.6000         0         0    8.6000         0         0         0
                      0   13.0000   26.2700   14.7300         0    9.3000   26.3100    1.5000    6.5200         0    1.3700    1.0100         0         0
                      0         0         0   73.8100    3.3100    0.0400   10.7700         0    0.0100    0.8400    1.4500    9.7700         0         0
                      0         0         0         0   43.9300         0    4.4300         0         0    0.0100         0    0.0100   51.6200         0
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
cal.T0  = [1875        1170        1134        1075        1015         839];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [10.2100    2.2700    2.6100    2.7500    2.5900    1.0100];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [7.6600    4.5900    3.8600    3.3700    2.7600    2.6100];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [37.8000    4.9000    4.6000    5.0000   16.3000   10.8000];

% initial and final composition used in calibration
cal.c0 = [0.0920    0.2280    0.5290    0.0180    0.1290    0.0040    0.0030];
cal.c1 = [0.0010    0.0010    0.0010    0.0070    0.2640    0.7250    0.0030];

% specify melting point dependence on H2O
cal.dTH2O   = [1300        1450        1500        1650        1900        2400];                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  san  dps  aug  ulv  mgt  ilm  ens  fsl  qtz  wat
cal.rhox0   = [3200,4040,2680,2590,2560,3200,3465,3925,4760,4730,3300,3665,2540,1000]; % mem ref densities [kg/m3]
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
