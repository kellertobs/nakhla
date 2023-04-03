% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 9;
cal.nmem   = 13;
cal.nmsy   = 6;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O','H$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K','H'};
cal.memStr = {'for','fay','mgt','ulv','ens','hyp','aug','pig','ant','alb','san','qtz','wat'};
cal.msyStr = {'olv','spn','opx','cpx','fsp','qtz'};
cal.cmpStr = {'ano','bas','and','rhy','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
cal.mem_oxd    = [  42.7    0.0    0.0    0.0   57.3    0.0    0.0    0.0    0.0     % forsterite (for)
                    29.5    0.0    0.0   70.5    0.0    0.0    0.0    0.0    0.0     % fayalite (fay)
                    0.0     9.0    2.0   89.0    0.0    0.0    0.0    0.0    0.0     % magnetite (mgt)
                    0.0    38.5    2.0   34.0   25.5    0.0    0.0    0.0    0.0     % ulvospinel (ulv)
                    52.0    0.0    3.6   17.3   24.5    2.6    0.0    0.0    0.0     % enstatite (ens)
                    49.4    0.0    0.9   36.7   12.0    1.0    0.0    0.0    0.0     % hypersthene (hyp)
                    53.0    0.0    1.0   10.0   15.7   19.75   0.5    0.05   0.0     % augite (aug)
                    50.8    0.0    1.5   23.6    4.8   16.0    2.9    0.40   0.0     % pigeonite (pig)
                    44.4    0.0   35.8    0.0    0.0   19.3    0.5    0.0    0.0     % anorthite (ant)
                    67.3    0.0   20.2    0.0    0.0    0.5   11.5    0.5    0.0     % albite (alb)
                    64.8    0.0   18.3    0.0    0.0    0.0    0.5   16.4    0.0     % sanidine (san)
                    99.99   0.01   0.0    0.0    0.0    0.0    0.0    0.0    0.0     % quartz (qtz)
                     0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  0  0  0  0  0  0  0  0  0    % spinel (spn)
               0  0  0  0  1  1  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  0  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  0  1  1  1  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%               for     fay  Fe-spn Ti-spn  ens    hyp    aug    pig    ant     alb    kfs    qtz    wat
cal.cmp_mem =[   0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0  100.0    0.0    0.0    0.0    0.0    % anorthosite (ano)
                 3.5    1.5    0.1    3.7   12.7    1.5   17.0    1.5   50.0    8.5    0.0    0.0    0.0    % basalt (bas)
                 0.05   0.95   1.8    2.5    6.0   12.5    7.2   10.5   16.5   38.0    4.0    0.0    0.0    % andesite (and)
                 0.0    0.0    0.8    0.05   0.0    1.0    0.5    9.0    1.5   16.6   25.25  45.3    0.0    % rhyolite (rhy)
                 0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0  100.0];  % hydrous fluid (fld)
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
cal.T0(cal.for) =  1550;
cal.T0(cal.bas) =  1200;
cal.T0(cal.and) =  1100;
cal.T0(cal.rhy) =  850;

% set first coeff. for P-dependence of T_m^i [GPa]
cal.A(cal.for)  =   6.1;
cal.A(cal.bas)  =   4.7;
cal.A(cal.and)  =   2.85;
cal.A(cal.rhy)  =   2.7;

% set second coeff. for P-dependence of T_m^i [1]
cal.B(cal.for)  =  8.9;
cal.B(cal.bas)  =  3.3;
cal.B(cal.and)  =  2.5;
cal.B(cal.rhy)  =  2.5;

% set entropy gain of fusion DeltaS [J/K]
cal.dS          =  300;

% set coeff. for T-dependence of partition coefficients K^i [1/K]
cal.r(cal.for)  =  24.0;
cal.r(cal.bas)  =  22.0;
cal.r(cal.and)  =  20.0;
cal.r(cal.rhy)  =  16.0;

% specify melting model phase diagram parameters
% cal.cphs0   = cal.cmp_oxd(1,1)/100; % phase diagram lower bound composition [wt SiO2]
% cal.cphs1   = cal.cmp_oxd(4,1)/100; % phase diagram upper bound composition [wt SiO2]
% cal.Tphs0   = 850;                  % phase diagram lower bound temperature [degC]
% cal.Tphs1   = 1890;                 % phase diagram upper bound temperature [degC]
% cal.PhDg    = [8.5,2.2,1.3,1.3];    % phase diagram curvature factor (> 1)
% cal.perCm   = cal.cmp_oxd(3,1)/100; % peritectic liquidus composition [wt SiO2]
% cal.perCx   = cal.cmp_oxd(2,1)/100; % peritectic solidus  composition [wt SiO2]
% cal.perTm   = 1120;                 % peritectic temperature [degC]
% cal.perTx   = 1150;                 % peritectic temperature [degC]
% cal.clap    = 1e-7;                 % Clapeyron slope for P-dependence of melting T [degC/Pa]

% specify melting point dependence on H2O
cal.dTH2O   = 1400;                 % solidus shift from water content [degC/wt^pH2O]
cal.pH2O    = 0.75;                 % solidus shift from water content [degC/wt^pH2O]

% specify geochemical model parameters
cal.nte     = 4;                    % number of trace elements
cal.nir     = 2;                    % number of isotope ratios
cal.Kte_mem = [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
cal.rhox0   = [3270,4390,7700,4300,3200,3500,3250,3350,2730,2620,2520,2650]; % mineral end-member reference densities [kg/m3]
cal.rhof0   = 1000;                 % fluid reference density [kg/m3]
cal.aT      = 4e-5;                 % thermal expansivity [1/K]
cal.gH      = 0.75;                 % hydrous melt expansivity [1/(wt H2O)]
cal.bP      = 1e-8;                 % fluid compressibility [1/Pa]

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
