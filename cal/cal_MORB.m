% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 14;
cal.nmsy   = 7;
cal.ncmp   = 8;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ant','alb','san','dps','aug','ulv','mgt','hyp','fsl','ilm','qtz','wat'};
cal.msyStr = {'olv','fsp','cxp','spn','opx','ilt','qtz'};
cal.cmpStr = {'frt','fyt','ogb','gbn','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                 SiO2      TiO2     Al2O3       FeO       MgO       CaO      Na2O       K2O       H2O
cal.mem_oxd = [42.4300         0         0    1.9000   55.6700         0         0         0         0
               30.8400         0         0   63.4400    5.7200         0         0         0         0
               44.4400         0   35.7100         0         0   19.0900    0.7600         0         0
               66.5700         0   20.8900         0         0    1.7000   10.7900    0.0500         0
               67.8900         0   19.4200         0         0    0.1800    5.8700    7.6400         0
               53.4300    0.5300    3.0800    2.7600   21.3200   18.8800         0         0         0
               51.5200    0.3100    0.2300   25.2300    4.6600   15.4900    2.5600         0         0
                     0   43.8500    2.7700   24.9500   28.4300         0         0         0         0
                     0    5.9800    0.3200   93.7000         0         0         0         0         0
               52.7700         0    3.9700   13.4500   26.7500    3.0600         0         0         0
               48.5300         0    0.2100   41.9000    8.3600    1.0000         0         0         0
                     0   51.2400         0   48.7600         0         0         0         0         0
              100.0000         0         0         0         0         0         0         0         0    % quartz (qtz)
                     0         0         0         0         0         0         0         0   100.000];  % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  1  1  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  0  0  0  0  0    % spinel (spn)
               0  0  0  0  0  0  0  0  0  1  1  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  0  0  0  0  1  0  0    % ilmenite (ilt)
               0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
cal.cmp_mem = [ 100.0    0       0       0       0       0       0       0       0       0       0       0       0       0
                 0    100.0    8.05    1.83   25.22   21.42   22.30    0.74    6.06    3.01    4.97       0       0       0
                 3.12    3.15    8.05    1.83   25.22   21.42   22.30    0.74    6.06    3.01    4.97       0       0       0
                 3.12    3.15    8.05    1.83   25.22   21.42   22.30    0.74    6.06    3.01    4.97       0       0       0
                    0    0.66       0    2.29    4.78       0   70.72    0.97    3.59       0    0.10    0.10       0       0
                    0    0.66       0    2.29    4.78       0   70.72    0.97    3.59       0    0.10    0.10       0       0
                    0       0       0       0    5.38       0       0       0    0.10       0       0    0.10   41.52       0
                    0       0       0       0       0       0       0       0       0       0       0       0       0  100.00]; % volatiles (vol)
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
cal.T0  = [1880  1250  1177  1083 1000 975  783];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8.71 4.0 4.36  3.97 3.20 2.97  2.52];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [42.20 25.0 2.10  10.50 13.5 12.40  11.30];

% initial composition used in calibration
cal.c0 = [0.17  0.05  0.23  0.08  0.02 0.02 0.01  0.003];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  dps  mau  fau  ant  alb  san  ens  fsl  ulv  mgt  ilm  qtz  wat
cal.rhox0   = [3150,4055,3200,3385,3475,2680,2580,2555,3240,3665,3925,4730,4875,2540,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 1000;                 % fluid reference density [kg/m3]

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
