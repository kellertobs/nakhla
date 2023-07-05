% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)
clear cal;

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 15;
cal.nmsy   = 6;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','hyp','fsl','mau','aug','fau','ulv','mgt','tim','ant','alb','san','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'dun','pxn','bas','eut','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
cal.mem_oxd    = [  42.7    0.0    0.0     0.0    57.3     0.0     0.0    0.0    0.0     % forsterite (for)
                    29.5    0.0    0.0    70.5     0.0     0.0     0.0    0.0    0.0     % fayalite (fay)
                    52.45   0.0    3.48   16.29   24.93    2.84    0.0    0.0    0.0     % hypersthene (hyp)
                    49.80   0.0    2.51   28.84   16.02    2.83    0.0    0.0    0.0     % ferrosillite (fsl)

                   53.87    0.64   3.13       0   20.29   21.53    0.54   0.0    0.0     % Mg-augite (mau)
                   51.58    0.22   0.46   27.41   10.37    8.45    1.51   0.0    0.0     %    augite (aug)
                   48.33    0.21   1.25   37.34       0   11.05    1.82   0.0    0.0     % Fe-augite (fau)

                    0.0    39.82   2.68   29.29  28.21     0.0     0.0    0.0    0.0     % ulvospinel (ulv)
                    0.0     4.00   6.06   87.99   1.98     0.0     0.0    0.0    0.0     % magnetite (mgt)
                    0.0    12.33   1.58   85.72   0.37     0.0     0.0    0.0    0.0     % titano-magnetite (tim)

                    44.4    0.0    35.8    0.0     0.0    19.3     0.5    0.0    0.0     % anorthite (ant)
                    67.3    0.0    20.2    0.0     0.0     0.5    11.5    0.0    0.0     % albite (alb)
                    64.8    0.0    18.3    0.0     0.0     0.0     0.5   16.4    0.0     % sanidine (san)
                   100.0    0.0     0.0    0.0     0.0     0.0     0.0    0.0    0.0     % quartz (qtz)
                     0.0    0.0     0.0    0.0     0.0     0.0     0.0    0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  1  0  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   for           fay       hyp       fsl       mau       aug       fau      ulv       mgt       tim       ant       alb       san       qtz       wat
cal.cmp_mem =   [         0         0         0         0         0         0         0         0         0         0  100.0000         0         0         0         0
                     4.3759    1.4908    0.7864    0.0511   30.6516    0.9666    5.4939    4.6147    0.0217    0.3825   37.4736   13.6911    0.0000    0.0000         0
                     0.4517    2.0412    0.2839    3.3047    1.5665   29.3598    5.8278    0.0459    5.9759    1.2507   18.0104   31.5913    0.2903    0.0000         0
                     0.0031    0.0533    0.0000    0.0000    0.5037    0.0828   26.1125    0.0000    0.0391    0.1463    4.1985   14.6359    8.8559   45.3689         0
                          0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000 ];
% cal.cmp_mem =   [         0         0         0         0         0         0         0         0         0         0  100.0000         0         0         0         0
%                     50.1313   11.3399    0.0000    0.0000    2.2913    0.0000    0.0000   10.1192    0.0002    0.0000   26.1009    0.0172    0.0000    0.0000         0
%                      0.0334   80.9127    1.6821    0.0430    0.0000    0.0916    0.4225    0.1461    7.3995    0.0681    8.0691    0.0518    1.0800    0.0000         0
%                      0.0000    0.0053    0.0014    0.0220    0.1681    0.0003    0.0255    0.0000    0.0004    0.0023    0.0010    0.0300    0.0135   99.7302         0
%                           0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000];
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
cal.T0(cal.dun) =  1553;
cal.T0(cal.pxn) =  1140;
cal.T0(cal.bas) =  1040;
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
