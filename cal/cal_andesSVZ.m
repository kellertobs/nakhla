% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 15;
cal.nmsy   = 6;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','ulv','mgt','ilm','ens','hyp','cp1','aug','pig','ant','alb','san','qtz','wat'};
cal.msyStr = {'olv','oxs','opx','cpx','fsp','qtz'};
cal.cmpStr = {'ano','bas','and','rhy','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
cal.mem_oxd    = [  42.7    0.0    0.0    0.0   57.3    0.0    0.0    0.0    0.0     % forsterite (for)
                    29.5    0.0    0.0   70.5    0.0    0.0    0.0    0.0    0.0     % fayalite (fay)
                    0.0    32.0    2.6   46.2   19.2    0.0    0.0    0.0    0.0     % ulvospinel (ulv)
                    0.0    20.6    0.3   77.0    2.1    0.0    0.0    0.0    0.0     % magnetite (mgt)
                    0.0    42.6    5.0   52.3    0.1    0.0    0.0    0.0    0.0     % ilmenite (ilm)
                    51.5    0.1    3.4   20.0   22.5    2.5    0.0    0.0    0.0     % enstatite (ens)
                    49.4    0.0    1.0   36.5   12.1    1.0    0.0    0.0    0.0     % hypersthene (hyp)
                    52.2    0.2    2.1   10.5   14.8   19.67   0.5    0.03   0.0     % Al-augite (cp1)
                    53.0    0.2    0.8   13.1   13.4   18.7    0.9    0.08   0.0     % augite (aug)
                    51.3    0.5    1.32  21.9    5.9   16.0    2.8    0.38   0.0     % pigeonite (pig)
                    44.4    0.0   35.8    0.0    0.0   19.3    0.5    0.0    0.0     % anorthite (ant)
                    67.3    0.0   20.2    0.0    0.0    0.5   11.5    0.5    0.0     % albite (alb)
                    64.8    0.0   18.3    0.0    0.0    0.0    0.5   16.4    0.0     % sanidine (san)
                    99.99   0.01   0.0    0.0    0.0    0.0    0.0    0.0    0.0     % quartz (qtz)
                     0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0  0  0  0    % oxide (oxd)
               0  0  0  0  0  1  1  0  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  0  1  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   for       fay       ulv       mgt       ilm       ens       hyp       cp1       aug       pig       ant       alb       kfs       qtz       wat
cal.cmp_mem =   [      0         0         0         0         0         0         0         0         0         0  100.0000         0         0         0         0
                  6.7542    8.9546    3.9788    1.8484    0.2356   11.8743    2.2291   16.3188    2.2568    0.5971   34.6419   10.3104    0.0000    0.0000         0
                  0.5298    4.0581    0.3463    2.0716    0.1853    0.3080    6.7981    0.0808    5.2206    1.8054   10.6887   65.0505    1.8161    1.0405         0
                  0.0000    0.0000    0.0000    0.0222    0.1358    0.0000    0.0000    0.0000    0.0000    1.1779    0.0000   27.5962   28.1076   42.9603         0
                       0         0         0         0         0         0         0         0         0         0         0         0         0         0  100.0000 ];  % hydrous fluid (fld)
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
cal.T0(cal.ano) =  1450;
cal.T0(cal.bas) =  1180;
cal.T0(cal.and) =  980;
cal.T0(cal.rhy) =  860;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.ano)  =   6.1;
cal.A(cal.bas)  =   4.7;
cal.A(cal.and)  =   2.85;
cal.A(cal.rhy)  =   2.7;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.ano)  =  8.9;
cal.B(cal.bas)  =  3.3;
cal.B(cal.and)  =  2.5;
cal.B(cal.rhy)  =  2.5;

% set entropy gain of fusion DeltaS [J/K]
cal.dS          =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.ano)  =  20.0;
cal.r(cal.bas)  =  18.0;
cal.r(cal.and)  =  16.0;
cal.r(cal.rhy)  =  14.0;

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