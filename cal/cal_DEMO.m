% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 11;
cal.nmsy   = 5;
cal.ncmp   = 7;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ant','alb','san','dps','aug','ulv','mgt','qtz','wat'};
cal.msyStr = {'olv','fsp','cpx','spn','qtz'};
cal.cmpStr = {'dun','tro','ogb','fbs','tra','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1   2    3    4   5   6   7    8   9]; % oxide indices for viscosity, density functions

%           plg  olv  cpx  qtz
cal.imsy = [ 2    1    3    5];  % mineral system indices for plotting basalt tetrahedron

% oxide composition of mineral end-members
%                SiO2 TiO2 Al2O3 FeO  MgO  CaO Na2O  K2O  H2O
cal.mem_oxd    = [ 43    0    0    0   57    0    0    0    0      % forsterite (ant) 
                   30    0    0   70    0    0    0    0    0      % fayalite (ant) 

                   44    0   36    0    0   20    0    0    0      % anorthite (ant)
                   69    0   19    0    0    0   12    0    0      % albite (ant)
                   67    0   18    0    0    0    3   12    0      % sanidine (san)

                   54    0    3    3   21   19    0    0    0      % diopside (dps)
                   49    3    1   29    0   14    4    0    0      % augite (aug)

                    0   38    3   31   28    0    0    0    0      % ulvospinel (ulv)
                    0    6    1   93    0    0    0    0    0      % magnetite (mgt)

                  100    0    0    0    0    0    0    0    0      % quartz (qtz)
                    0    0    0    0    0    0    0    0  100 ];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
%              for fay ant alb san dps aug ulv mgt qtz wat
cal.msy_mem = [ 1   1   0   0   0   0   0   0   0   0   0    % olivine (olv) 
                0   0   1   1   1   0   0   0   0   0   0    % feldspar (fsp)
                0   0   0   0   0   1   1   0   0   0   0    % clinopyroxene (cpx)
                0   0   0   0   0   0   0   1   1   0   0    % spinel (spn)
                0   0   0   0   0   0   0   0   0   1   0];  % quartz (qtz)

% mineral end-member composition of melting model components
%              for fay  ant  alb  san  dps  aug  ulv  mgt  qtz  wat
cal.cmp_mem = [ 95   5    0    0    0    0    0    0    0    0    0    % dunite (dun)
                15   5   60   20    0    0    0    0    0    0    0    % troctolite (tro)
                 2   2   36    7    0   45    8    0    0    0    0    % olivine-gabbro (ogb)
                 1  12   18   25    0    8   27    8    1    0    0    % ferro-basalt (fbs)
                 0   0    8   72    3    2    8    2    3    2    0    % trachy-andesite (tra)
                 0   0    1    8   40    0    1    0    2   48    0    % rhyolite (rhy)
                 0   0    0    0    0    0    0    0    0    0  100];  % volatile (vol)
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
cal.T0  = [1850  1225  1175  1070  965  850];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [7.0  5.0  4.0  3.0  2.5  2.0];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [7.0  5.0  4.0  3.0  2.5  2.0];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r   = [40  5  6  9  8  5];

% set entropy gain of fusion and evaporation DeltaS [J/K]
cal.Dsx  = 350;
cal.Dsf  = 450;

% specify melting point dependence on H2O
cal.dTH2O = [910  1370  1430  1570  1740  1980];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;               % solidus shift from water content exponent

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  san  dps  aug  ulv  mgt  qtz  wat
cal.rhox0   = [3200,4400,2700,2600,2650,3300,3500,3930,4750,2500,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%              for  fay  ant  alb  san  dps  aug  ulv  mgt  qtz  wat
cal.etax0  = [1e18,1e18,1e17,1e17,1e17,1e19,1e19,1e17,1e17,1e19,1e0]; % mem ref viscosities [Pas]
cal.etaf0  = 0.1;                    % fluid viscosity constant [Pas]
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
cal.tol     = 1e-12;
cal.alpha   = 0.5;
