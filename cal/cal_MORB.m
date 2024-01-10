% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 13;
cal.nmsy   = 6;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ant','alb','dps','aug','tms','mgt','ilm','hyp','fsl','qtz','wat'};
cal.msyStr = {'olv','fsp','cxp','oxs','opx','qtz'};
cal.cmpStr = {'dun','tro','osg','fbs','ban','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2    TiO2    Al2O3    FeO     MgO     CaO    Na2O     H2O
cal.mem_oxd    = [  41.57    0.0     0.0     6.45   51.98    0.0     0.0     0.0     % forsterite (for)
                    38.34    0.0     0.0    23.66   38.00    0.0     0.0     0.0     % fayalite (fay)

                    44.71    0.0    35.62    0.0     0.0    18.96    0.71    0.0     % anorthite (ant)
                    67.30    0.0    20.42    0.0     0.0     1.15   11.13    0.0     % albite (alb)

                    53.92    0.0     2.64    2.36   21.27   19.81    0.0     0.0     % diopside (dps)
                    51.46    0.0     0.46   25.88    4.44   15.24    2.52    0.0     % augite (aug)

                     0.0    48.58    2.19   21.48   27.75    0.0     0.0     0.0     % Ti-Mg-bearing spinel (tms)
                     0.0     0.0     0.48   99.52    0.0     0.0     0.0     0.0     % magnetite (mgt)
                     0.0    51.38    0.0    48.62    0.0     0.0     0.0     0.0     % ilmenite (ilm)

                    53.17    0.0     3.58   13.04   27.09    3.12    0.0     0.0     % hypersthene (hyp)
                    47.75    0.0     0.83   43.11    7.07    1.24    0.0     0.0     % ferrosillite (fsl)

                   100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0     % quartz (qtz)
                     0.0     0.0     0.0     0.0     0.0     0.0     0.0   100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  1  1  0  0  0  0  0  0  0    % clinopyroxene (cxn)
               0  0  0  0  0  0  1  1  1  0  0  0  0    % oxides (oxs) 
               0  0  0  0  0  0  0  0  0  1  1  0  0    % orthopyroxene (oxn)
               0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                for     fay     ant     alb     dps     aug     tms     mgt     ilm     hyp     fsl     qtz     wat
cal.cmp_mem = [84.22   15.78       0       0       0       0       0       0       0       0       0       0       0
               15.59   13.12   71.30       0       0       0       0       0       0       0       0       0       0
                   0    4.83   37.27    9.62   45.86       0    2.42       0       0       0       0       0       0
                   0       0   20.93   23.23   30.64    6.04    2.97   12.51       0    3.67       0       0       0
                   0       0   12.59   34.80       0   26.40       0       0    1.47    3.26   21.48       0       0
                   0       0       0   49.32       0       0       0       0       0       0    3.11   47.57       0
                   0       0       0       0       0       0       0       0       0       0       0       0  100.00];
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
cal.T0  = [1880  1203  1162  1102  1045  854];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [9 5 4.3 3.5 3 2.5];

% set entropy gain of fusion DeltaS [J/K]
cal.dS  = 350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [36.7  1.6  3.3  11.6  16.2  9.9];

% initial composition used in calibration
cal.c0 = [0.10  0.13  0.28  0.35  0.06  0.08  0.005];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;0.9;3.00;10.0;1.1].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  ant  alb  dps  aug  tms   mgt  ilm  hyp  fsl  qtz  wat
cal.rhox0   = [3160,3400,2680,2560,3200,3480,3890,4980,4700,3280,3660,2550,1000]; % mineral end-member reference densities [kg/m3]
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
