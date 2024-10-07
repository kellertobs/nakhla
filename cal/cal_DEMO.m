% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 7;
cal.nmem   = 8;
cal.nmsy   = 4;
cal.ncmp   = 4;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ant','alb','dps','aug','qtz','wat'};
cal.msyStr = {'olv','fsp','cpx','qtz'};
cal.cmpStr = {'dun','bas','rhy','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1        3    4   5   6   7        9]; % oxdie indices for viscosity, density functions

% oxide composition of mineral end-members
%                   SiO2   Al2O3     FeO     MgO     CaO    Na2O     H2O
cal.mem_oxd    = [ 43.00       0       0   57.00       0       0       0      % forsterite (ant) 
                   30.00       0   70.00       0       0       0       0      % fayalite (ant) 
                   44.00   36.00       0       0   20.00       0       0      % anorthite (ant)
                   67.00   21.00       0       0       0   12.00       0      % albite (ant)
                   54.00       0    8.00   18.00   20.00       0       0      % diopside (dps)
                   48.00       0   26.00    8.00   18.00       0       0      % augite (aug)
                  100.00       0       0       0       0       0       0      % quartz (qtz)
                       0       0       0       0       0       0  100.00 ];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
%              for fay ant alb dps aug qtz wat
cal.msy_mem = [ 1   1   0   0   0   0   0   0    % olivine (olv) 
                0   0   1   1   0   0   0   0    % feldspar (fsp)
                0   0   0   0   1   1   1   0    % pyroxene (pxn)
                0   0   0   0   0   0   1   0];  % quartz (qtz)

% mineral end-member composition of melting model components
%               for    fay    ant    alb    dps    aug    qtz    wat
cal.cmp_mem = [100.0   0.0    0.0    0.0    0.0    0.0    0.0    0.0    % dunite (dun)
                 0.0  10.0   35.0    0.0   45.0    0.0    0.0    0.0    % gabbro (gbr)
                 0.0   0.0    0.0   50.0    0.0    5.0   45.0    0.0    % rhyolite (rhy)
                 0.0   0.0    0.0    0.0    0.0    0.0    0.0  100.0];  % volatile (vol)
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
cal.T0  = [1890  1150  850];

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A   = [6  4  2];

% set second coeff. for P-dependence of T_m^i [1]
cal.B   = [6  4  2];

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r   = [40 20 10];

% set entropy gain of fusion and evaporation DeltaS [J/K]
cal.Dsx  = 350;
cal.Dsf  = 450;

% specify melting point dependence on H2O
cal.dTH2O = [900  1500  2000];  % solidus shift from water content prefactor [K/wt^pH2O]
cal.pH2O  = 0.75;               % solidus shift from water content exponent

% specify geochemical model parameters
cal.ntrc     = 6;                    % number of trace elements
cal.trcStr   = {'K 0.01','K 0.10','K 1.0','K 3.00','K 10.0','K 1.0'};
cal.Ktrc_mem = [0.01;0.10;1.0;3.0;10.0;1.0].*ones(cal.ntrc,cal.nmem);

% specify density parameters
%               for  fay  ant  alb  dps  aug  qtz  wat
cal.rhox0   = [3200,4400,2700,2600,3300,3500,2500,1000]; % mem ref densities [kg/m3]
cal.rhof0   = 1000;                 % fluid ref density [kg/m3]

% specify three-phase coefficient model parameters
%              for  fay  ant  alb  dps  aug  qtz  wat
cal.etax0  = [1e18,1e18,1e17,1e17,1e19,1e19,1e19,1e0]; % mem ref viscosities [Pas]
cal.etaf0  = 0.1;                  % fluid viscosity constant [Pas]
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
cal.tol     = 1e-9;
cal.alpha   = 0.5;
