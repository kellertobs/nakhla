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
cal.memStr = {'for','fay','dps','mau','fau','ant','alb','san','ens','fsl','ulv','mgt','ilm','qtz','wat'};
cal.msyStr = {'olv','cxp','fsp','opx','oxs','qtz'};
cal.cmpStr = {'dun','ogb','gbn','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2     TiO2    Al2O3    FeO     MgO     CaO    Na2O     K2O     H2O
cal.mem_oxd    = [ 41.3800         0         0    6.9500   51.4400    0.2300         0         0         0   % for
                   30.6400         0         0   64.2800    5.0000    0.0800         0         0         0   % fay

                   52.7500    0.5000    3.9400    1.4800   21.2700   20.0600         0         0         0   % dps
                   52.8300    0.1600    0.8900   16.7100   11.8000   16.2100    1.4000         0         0   % mau
                   50.8900    0.7200    0.5400   26.8400    2.2800   14.8700    3.6000    0.2600         0   % fau
   
                   44.7200         0   35.5400         0         0   18.8900    0.8500         0         0   % ant
                   66.9900         0   20.6100         0         0    1.3700   10.9800    0.0500         0   % alb
                   68.5600         0   18.0700         0         0         0    7.6500    5.7200         0   % san
   
                   53.2600    0.2900    5.0400    7.3300   30.9000    3.1200    0.0600         0         0   % ens
                   48.9100         0    0.2700   40.2700    9.4700    1.0000    0.0800         0         0   % fsl
         
                         0   36.3200    3.2400   32.0400   28.4000         0         0         0         0   % ulv
                         0    0.2000    1.7600   98.0400         0         0         0         0         0   % mgt
                         0   53.7300         0   46.2700         0         0         0         0         0   % ilm

                  100.0000         0         0         0         0         0         0         0         0   % quartz (qtz)
 
                         0         0         0         0         0         0         0         0  100.000];  % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  1  1  1  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  1  1  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  1  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                for     fay     dps     mau     fau     ant     alb     san     ens     fsl     ulv     tmg     ilm     qtz     wat
cal.cmp_mem = [ 100.0000         0         0         0         0         0         0         0         0         0         0         0         0         0         0
                 25.4942    5.3840   33.0092         0         0   36.1126         0         0         0         0         0         0         0         0         0
                  2.2827    6.7650   11.9621   19.2151         0   28.2181   15.9709         0   12.7563         0    2.8299         0         0         0         0
                  0.1042    1.1218   13.2949    0.5175   10.9337   14.7233   37.8258    3.2719    0.1009   11.0633    1.4858    5.5568         0         0         0
                       0    0.1450         0    1.9989    7.8547         0   57.3760   31.0775    0.8468    0.4540         0    0.1080    0.1393         0         0
                       0         0         0         0    1.9731         0         0   45.4803         0    2.1190         0         0    0.1017   50.3259         0
                       0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000]; % volatiles (vol)
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
cal.T0  = [1880  1203  1162  1045  950  854];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [9 5 4.3 3.5 3.0 2.5];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [36.7  1.6  3.3  10.0 11.6  9.9];

% initial composition used in calibration
cal.c0 = [0.10  0.13  0.28  0.35  0.08  0.005];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  ant  alb  dps  aug  fau,  tms   mgt  ilm  hyp  fsl  qtz  wat
cal.rhox0   = [3160,3400,2680,2560,3200,3480,3520,3890,4980,4700,3280,3660,2550,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

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
