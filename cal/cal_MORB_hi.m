% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 13;
cal.nmsy   = 5;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ant','alb','san','dps','aug','pig','ulv','mgt','ilm','qtz','wat'};
cal.msyStr = {'olv','fsp','cxp','oxs','qtz'};
cal.cmpStr = {'dun','tro','gbr','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7   8   9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                SiO2    TiO2   Al2O3     FeO     MgO     CaO    Na2O     K2O     H2O
cal.mem_oxd = [ 41.37       0       0    7.57   51.06       0       0       0       0   % forsterite (for)
                29.79       0       0   69.09    1.12       0       0       0       0   % fayalite (fay)

                44.42       0   35.72       0       0   19.10    0.76       0       0   % anorthite (ant)
                68.78       0   19.33       0       0       0   11.83    0.06       0   % albite (alb)
                68.64       0   18.31       0       0       0    6.38    6.67       0   % sanidine (san)

                53.38    0.49    3.12    2.85   21.17   18.99       0       0       0   % diopside (dps)
                51.04       0    0.55   26.04    4.31   15.79    2.27       0       0   % augite (aug)
                49.90    2.39    0.73   29.63    0.07   13.40    3.88       0       0   % pigeonite (pig)

                    0   38.80    2.82   30.50   27.88       0       0       0       0   % ulvospinel (ulv)
                    0    5.82    1.47   92.71       0       0       0       0       0   % magnetite (mgt)
                    0   54.83       0   45.17       0       0       0       0       0   % ilmenite (ilm)

               100.00       0       0       0       0       0       0       0       0   % quartz (qtz)
                    0       0       0       0       0       0       0       0  100.00]; % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100; 

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0  0    % feldspar (fsp)
               0  0  0  0  0  1  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  0  1  1  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%               for   fay   ant   alb   san   dps   aug   pig   ulv   mgt   ilm   qtz   wat
cal.cmp_mem = [94.2   5.8     0     0     0     0     0     0     0     0     0     0     0   % dun
                3.6   2.3  67.1  27.1     0     0     0     0     0     0     0     0     0   % tro
                1.0   0.1  34.1   6.2   0.7  44.0   8.3     0   5.7     0     0     0     0   % gbr
                  0  15.4  18.3  30.6   0.1   4.1  26.3   2.3   0.7   1.4   0.7     0     0   % fbs 
                  0     0   4.0  80.8   8.5   0.1   0.4   5.8   0.1   0.1   0.2   0.1     0   % tra
                  0     0     0   5.4  45.1     0   4.5   1.3     0   0.1   0.2  43.4     0   % rhy
                  0     0     0     0     0     0     0     0     0     0     0     0 100.0]; % vol
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
cal.T0  = [1850  1190  1177  1071  993  830];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [6.6  5.0  4.9  3.4  2.8  1.9];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [6.6  5.0  4.9  3.4  2.8  1.9];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r  = [40.5  3.0  3.0  8.5  8.3  3.4];

% set entropy gain of fusion DeltaS [J/K]
cal.Dsx  = 350;
cal.Dsf  = 450;

% specify melting point dependence on H2O
cal.dTH2O = [908  1412  1428  1569  1692  2024];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;                                  % solidus shift from water content exponent

% primary and evolved end-member compositions used in calibration
cal.c0     = [0.104  0.185  0.388  0.270  0.008  0.044  0.005];

cal.c0_oxd = [50.00  1.26  15.08  9.10  10.58  11.31  2.49  0.18  0.50];

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  san  dps  aug  pig  ulv  mgt  ilm  qtz  wat
cal.rhox0   = [3210,4200,2680,2600,2550,3210,3470,3520,3930,4750,4700,2540,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%              for  fay  ant  alb  san  dps  aug  pig  ulv  mgt  ilm  qtz  wat
cal.etax0   = [1e19,1e19,1e17,1e17,1e17,1e19,1e19,1e19,1e17,1e17,1e17,1e19,1e0]; % mem ref viscosities [Pas]
cal.etaf0   = 0.1;                    % fluid viscosity constant [Pas]
cal.Eax     = 300e3;                  % solid viscosity activation energy [J/mol]
cal.AA      =[ 0.65, 0.25, 0.35; ...  % permission slopes
               0.20, 0.20, 0.20; ...  % generally numbers between 0 and 1
               0.20, 0.20, 0.20; ];   % increases permission slopes away from step function 

cal.BB      =[ 0.55, 0.18, 0.27; ...  % permission step locations
               0.64,0.012,0.348; ...  % each row sums to 1
               0.80, 0.12, 0.08; ];   % sets midpoint of step functions

cal.CC      =[[0.30, 0.30, 0.40]*0.7; ... % permission step widths
              [0.52, 0.40, 0.08]*1.1; ... % square brackets sum to 1, sets angle of step functions
              [0.15, 0.25, 0.60]*0.7; ];  % factor increases width of step functions

% convergence tolerance
cal.tol     = 1e-9;
cal.alpha   = 0.5;
