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
cal.memStr = {'for','fay','hyp','fsl','mau','aug','fau','ulv','mgt','tim','ilm','ant','alb','san','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','oxs','fsp','qtz'};
cal.cmpStr = {'ano','gbr','bas','rhy','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
%                   SiO2   TiO2   Al2O3   FeO     MgO      CaO    Na2O    K2O    H2O
cal.mem_oxd    = [  42.7    0.0    0.0     0.0    57.3     0.0     0.0    0.0    0.0     % forsterite (for)
                    29.5    0.0    0.0    70.5     0.0     0.0     0.0    0.0    0.0     % fayalite (fay)

                    52.21   0.0    3.75   15.75   25.03    3.26    0.0    0.0    0.0     % hypersthene (hyp)
                    47.04   0.0    1.34   44.77    5.69    1.16    0.0    0.0    0.0     % ferrosillite (fsl)

                    52.87   0.54   3.08    5.57   18.06   19.33    0.55   0.0    0.0     % Mg-augite (mau)
                    52.00   0.16   1.00   19.42   10.56   15.51    1.35   0.0    0.0     %    augite (aug)
                    49.86   0.52   2.05   27.21    2.25   15.36    2.75   0.0    0.0     % Fe-augite (fau)

                     0.0   40.33   2.16   28.96   28.55    0.0     0.0    0.0    0.0     % ulvospinel (ulv)
                     0.0    8.87   8.45   81.71    0.97    0.0     0.0    0.0    0.0     % magnetite (mgt)
                     0.0   15.14   1.83   81.14    1.89    0.0     0.0    0.0    0.0     % titano-magnetite (tim)
                     0.0   51.76   0.0    48.24    0.0     0.0     0.0    0.0    0.0     % ilmenite

                    44.4    0.0    35.8    0.0     0.0    19.3     0.5    0.0    0.0     % anorthite (ant)
                    67.3    0.0    20.2    0.0     0.0     0.5    11.5    0.0    0.0     % albite (alb)

                   100.0    0.0     0.0    0.0     0.0     0.0     0.0    0.0    0.0     % quartz (qtz)
                     0.0    0.0     0.0    0.0     0.0     0.0     0.0    0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  1  1  1  0  0  0  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  1  1  0  0  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  0  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   for           fay       hyp       fsl       mau       aug       fau      ulv       mgt       tim       ilm        ant       alb       qtz       wat
cal.cmp_mem =   [         0         0         0         0         0         0         0         0         0         0         0  100.0000         0         0         0
                     9.7374    2.4667    1.4334    0.0476   30.3819    1.1632    4.7503    6.2928    0.0193    0.9280    0.0000   29.0877   13.6917    0.0000         0
                     1.2091   10.8779    0.2359    5.0620    1.2086   16.3261   20.0191    0.0332    2.1568    1.0551    0.0000   13.7767   28.0394    0.0000         0
                     0.0039    0.1105    0.0000    0.5602    0.0842    0.1166    1.6404    0.0000    0.0407    0.3756    0.2103    0.9875   47.4679   48.4022         0
                          0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000];
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
cal.T0(cal.ano)  =  1553;
cal.T0(cal.gbr)  =  1118;
cal.T0(cal.bas)  =  1048;
cal.T0(cal.rhy)  =   850;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A  =  (cal.T0+273.15)./300;

% set second coeff. for P-dependence of T_m^i [1]
cal.B  =  0*cal.A + 1;

% set entropy gain of fusion DeltaS [J/K]
cal.dS          =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.ano)  =  20.4;
cal.r(cal.gbr)  =  10.3;
cal.r(cal.bas)  =  18.5;
cal.r(cal.rhy)  =  9.5;

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
