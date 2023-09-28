% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 11;
cal.nmsy   = 6;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','hyp','aug','pig','ant','alb','ilm','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','fsp','oxs','qtz'};
cal.cmpStr = {'dun','pxn','bas','eut','vol'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end

%           SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O
cal.ioxd = [   1    2     3   4   5   6    7       9]; % oxdie indices for viscosity, density functions

% set mineral end-member oxide compositions
%                 SiO2   TiO2   Al2O3   FeO     MgO   CaO   Na2O   H2O
cal.mem_oxd = [  42.7    0.05    0.0    0.0   57.25   0.0    0.0   0.0    % forsterite
                 30.0    0.20    0.0   69.0    0.0    0.8    0.0   0.0    % fayalite

                 56.0    0.08    4.49   4.33  33.92   1.11	 0.07  0.0    % enstatite
                 50.0    0.15    8.85  12.17  26.08   2.70   0.05  0.0    % hypersthene

                 51.0    0.08    7.16   2.88  22.25  16.43   0.20  0.0    % Mg-augite
                 47.6    1.68    2.05  33.27   5.19  10.08   0.13  0.0    % Fe-pigeonite

                 44.2    0.08   35.52   0.4    0.2   19.1    0.5   0.0    % anorthite
                 67.4    0.20   20.4    1.1    0.3    2.1    8.5   0.0    % albite

                 0.32   20.56    3.07  74.63   1.14   0.28   0.0   0.0    % ilmenite

                99.95    0.05    0.0    0.0    0.0    0.0    0.0   0.0    % quartz
                
                 0.0     0.0     0.0    0.0    0.0    0.0    0.0 100.0];  % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [ 1  1  0  0  0  0  0  0  0  0  0     % olivine (olv)
                0  0  1  1  0  0  0  0  0  0  0     % orthopyroxene (opx)
                0  0  0  0  1  1  0  0  0  0  0     % clinopyroxene (cpx)
                0  0  0  0  0  0  1  1  0  0  0     % feldspar (fsp)
                0  0  0  0  0  0  0  0  1  0  0     % oxide (oxd)
                0  0  0  0  0  0  0  0  0  1  0 ];  % quartz (qtz)

% set mineral proportions in eutectic component
cal.cmp_mem = [ 1.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000    0.0      % dunite
                0.1600	0.1800	0.5400	0.1200	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000    0.0      % pyroxenite
                0.0200	0.0700	0.0000	0.0500	0.2700	0.1000	0.3300	0.1000	0.0600	0.0000    0.0      % basalt
                0.0000	0.0000	0.0000	0.0000	0.0000	0.4500	0.1300	0.1400	0.1800	0.1000    0.0      % eutectic
                     0       0       0       0       0       0       0       0       0       0  100.0 ];   % hydrous fluid (fld)
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
cal.T0(cal.for)  =  1890;
cal.T0(cal.pxn)  =  1380;
cal.T0(cal.bas)  =  1090;
cal.T0(cal.eut)  =  1021;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.for)   =   6.1;
cal.A(cal.pxn)   =   4.7;
cal.A(cal.bas)   =   2.85;
cal.A(cal.eut)   =   2.7;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.for)   =  8.9;
cal.B(cal.pxn)   =  3.3;
cal.B(cal.bas)   =  2.5;
cal.B(cal.eut)   =  2.5;

% set entropy gain of fusion DeltaS [J/K]
cal.dS           =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.for)   =  22.0;
cal.r(cal.pxn)   =  18.0;
cal.r(cal.bas)   =   8.0;
cal.r(cal.eut)   =   6.0;

% specify melting point dependence on H2O
cal.dTH2O   = 1500;                 % solidus shift from water content [degC/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [degC/wt^pH2O]

% specify geochemical model parameters
cal.nte          =  4;  % number of trace elements
cal.nir          =  2;  % number of isotope ratios
cal.Kte_mem      =  [0.01,0.10,3.0,10.0].'.*ones(cal.nte,cal.nmem); % mineral trace element partition coefficients

% specify density parameters
cal.rhox0   = [3270,4390,3000,3100,3250,3250,2730,3260,4400,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 1000;                 % fluid reference density [kg/m3]

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
% cal.xi      = 3e-7;                 % solid viscosity level
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
