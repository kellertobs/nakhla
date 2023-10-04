% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 15;
cal.nmsy   = 6;
cal.ncmp   = 6;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','hyp','fsl','dps','mau','fau','ulv','mgt','ilm','ant','alb','san','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'ano','tro','gbn','trd','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   TiO2   Al2O3     FeO     MgO     CaO    Na2O     K2O    H2O
cal.mem_oxd    = [ 42.71    0.0     0.0     0.0    57.29    0.00    0.0     0.0    0.0     % forsterite (for)
                   29.49    0.0     0.0    70.51    0.0     0.00    0.0     0.0    0.0     % fayalite (fay)

                   51.22    0.0     3.24   21.65   21.15    2.74    0.0     0.0    0.0     % hypersthene (hyp)
                   49.04    0.0     0.30   40.38    9.33    0.95    0.0     0.0    0.0     % ferrosillite (fsl)

                   52.44    0.24    2.17    9.40   15.00   20.20    0.49    0.06   0.0     % diopside (dps)
                   52.62    0.0     0.72   16.22   11.39   18.01    1.02    0.01   0.0     % Mg-augite (mau)
                   51.76    0.90    0.89   23.06    4.73   15.33    3.03    0.30   0.0     % Fe-augite (fau)

                    0.0    36.63    2.50   44.14   16.74    0.0     0.0     0.0    0.0     % ulvospinel (msp)
                    0.0     0.80    0.56   98.64    0.0     0.0     0.0     0.0    0.0     % magnetite (mgt)
                    0.0    52.00    0.00   48.00    0.0     0.0     0.0     0.0    0.0     % ilmenite (ilm)

                   44.00    0.0    36.10    0.0     0.0    19.52    0.38    0.0    0.0     % anorthite (ant)
                   66.93    0.0    20.69    0.0     0.0     1.72   10.61    0.05   0.0     % albite (alb)
                   68.91    0.0    17.79    0.0     0.0     0.0     6.12    7.18   0.0     % sanidine (san)

                  100.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                    0.0     0.0     0.0     0.0     0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  1  0  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                  for     fay     hyp     fsl     dps     mau     fau     ulv     mgt     ilm     ant     alb     san     qtz     wat
cal.cmp_mem =   [    0       0       0       0       0       0       0       0       0       0   90.85    9.15       0       0       0
                 22.76   14.60       0       0       0       0       0    8.58       0       0   50.98    3.09       0       0       0
                  0.05   12.00   19.12    0.05   25.87    1.91    1.00    3.09    1.80       0   21.39   13.72       0       0       0
                     0    0.44    1.38    3.79    1.48    1.93    0.97    1.05    1.97    0.50    0.70   72.63   13.16       0       0
                     0       0       0       0       0    4.14    1.59       0       0       0       0    0.13   55.22   38.92       0
                     0       0       0       0       0       0       0       0       0       0       0       0       0       0  100.00];
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
cal.T0 =  [1530  1143  1109  978  741];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = (cal.T0+273.15)./350;

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [8 4.5 4 3 2.5];

% set entropy gain of fusion DeltaS [J/K]
cal.dS =  350;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  =  [27.8  2.0  8.3  13.2  16.9];

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [K/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [K/wt^pH2O]

% specify geochemical model parameters
cal.ntrc    = 6;                    % number of trace elements
cal.trcStr  = {'K 0.01','K 0.10','K 0.9','K 3.00','K 10.0','K 1.1'};
cal.Ktrc_mem = [0.01;0.10;0.9;3.00;10.0;1.1].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%              for  fay  hyp  fsl  dps  may  fay  ulv  mgt  ilm  ant  alb  san  qtz  wat
cal.rhox0   = [3270,4390,3270,3600,3250,3320,3450,4550,4980,4930,2690,2570,2520,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 500;                  % fluid reference density [kg/m3]

% specify three-phase coefficient model parameters
cal.etax0 = 1e18;                   % solid reference viscosity [Pas]
cal.etaf0 = 1e-1;                   % vapour reference viscosity [Pas]
cal.Eax   = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA  =  [ 0.65, 0.25, 0.35; ...  % permission slopes
             0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
             0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB  =  [ 0.55, 0.18, 0.27; ...  % permission step locations
             0.64,0.012,0.348; ...  % each row sums to 1
             0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC  =  [[0.30, 0.30, 0.40]*0.7; ... % permission step widths
            [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
            [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

% % specify mixture viscosity parameters (Costa et al., 2009)
% cal.Bphi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient bubbles
% cal.Bchi    = 2.0;                  % Einstein-Roscoe powerlaw coefficient crystals
% cal.chi_pck = 0.60;                 % rheologically critical crystal fraction
% cal.gamma   = 2.50;                 % step-function steepness coefficient
% cal.delta   = 27;                   % solid viscosity melt-weakening slope
% cal.xi      = 4.5e-4;               % solid viscosity level
% cal.etaf0   = 0.1;                  % fluid viscosity constant
% 
% % specify segregation coefficient parameters
% cal.bm      = 50;                   % melt permeability geometric factor (k0 = dx^2/bm)
% cal.cm      = 0.001;                % melt percolation threshold
% cal.nm      = 3;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
% cal.mm      = 2;                    % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
% cal.bf      = 50;                   % fluid permeability geometric factor (k0 = dx^2/bm)
% cal.cf      = 0.05;                 % fluid percolation threshold
% cal.nf      = 4;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
% cal.mf      = 2;                    % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
