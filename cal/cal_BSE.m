% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 7;
cal.nmem   = 11;
cal.nmsy   = 5;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','hyp','fsl','aug','pig','ant','alb','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','fsp','qtz'};
cal.cmpStr = {'dun','pxn','bas','eut','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
cal.mem_oxd    = [  42.7    0.0     0.0    57.3     0.0     0.0    0.0     % forsterite (for)
                    29.5    0.0    70.5     0.0     0.0     0.0    0.0     % fayalite (fay)
                    55.03   3.47    9.71   30.67    1.12    0.0    0.0     % enstatite (ens)
                    53.77   2.04   16.71   24.70    2.77    0.0    0.0     % hypersthene (hyp)
                    48.98   1.79   35.43   10.90    2.91    0.0    0.0     % ferrosillite (fsl)

                   53.49    1.36   10.34   17.06   17.60    0.15   0.0     % augite (aug)
                   51.18    0.95   22.11    9.56   15.77    0.41   0.0     % pigeonite (pig)

                    44.4    35.8    0.0     0.0    19.3     0.5    0.0     % anorthite (ant)
                    67.3    20.2    0.0     0.0     0.5    11.5    0.0     % albite (alb)
                   100.0     0.0    0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                     0.0     0.0    0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  1  1  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   for           fay        ens       hyp      fsl       aug       pig       ant       alb       qtz       wat
cal.cmp_mem =   [   100.0000         0         0         0         0         0         0         0         0         0         0
                      0.4319   97.5873    1.1889    0.7919    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000         0
                      0.0000    0.0000    0.7777    4.3752    0.5792   15.6742    7.2374   69.4720    1.8845    0.0000         0
                      0.0000    0.0000    0.0000    0.0986    0.6145    0.1292   30.4438   10.6900    7.6699   50.3539         0
                           0         0         0         0         0         0         0         0         0         0  100.0000 ];
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
cal.T0(cal.dun) =  1890;
cal.T0(cal.pxn) =  1250;
cal.T0(cal.bas) =  1080;
cal.T0(cal.eut) =   950;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.dun)  =   6.1;
cal.A(cal.pxn)  =   4.7;
cal.A(cal.bas)  =   2.85;
cal.A(cal.eut)  =   2.7;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.dun)  =  8.9;
cal.B(cal.pxn)  =  3.3;
cal.B(cal.bas)  =  2.5;
cal.B(cal.eut)  =  2.5;

% set entropy gain of fusion DeltaS [J/K]
cal.dS          =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.dun)  =  20.0;
cal.r(cal.pxn)  =  10.1;
cal.r(cal.bas)  =  19.0;
cal.r(cal.eut)  =  9.5;

% specify melting point dependence on H2O
cal.dTH2O   = 1500;                 % solidus shift from water content [degC/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [degC/wt^pH2O]

% specify geochemical model parameters
cal.nte     = 4;                    % number of trace elements
cal.nir     = 2;                    % number of isotope ratios
cal.Kte_mem = [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
cal.rhox0   = [3270,4390,4300,7500,5600,3200,3500,3250,3300,3350,2730,2620,2520,2650,1000]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 1000;                 % fluid reference density [kg/m3]

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
