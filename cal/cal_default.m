% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 8;
cal.nmem   = 8;
cal.nmsy   = 5;
cal.ncmp   = 4;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O'};
     elStr = {'Si','Ti','Al','Fe','Mg','Ca','Na','K'};
cal.memStr = {'for','fay','ilm','opx','cpx','ant','alb','qtz'};
cal.msyStr = {'olv','oxd','pxn','plg','qtz'};
cal.cmpStr = {'dun','gbr','bas','eut'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
cal.mem_oxd = [  42.7    0.0    0.0    0.0   57.3    0.0    0.0   0.0      % forsterite (for)
                 29.5    0.0    0.0   70.5    0.0    0.0    0.0   0.0      % fayalite (fay)
                 0.0    18.5    2.0   79.0    0.5    0.0    0.0   0.0      % ilmenite (ilm)
                 54.9    0.2    7.5    3.0   33.9    0.4    0.1   0.0      % orthopyroxene (opx)
                 49.5    1.0    1.0   30.0    1.7   15.3    1.5   0.0      % clinopyroxene (cpx)
                 44.4    0.0   35.8    0.0    0.0   19.2    0.6   0.0      % anorthite (ant)
                 67.4    0.0   20.2    0.0    0.0    0.8   11.6   0.0      % albite (alb)
                 99.99   0.01   0.0    0.0    0.0    0.0    0.0   0.0];    % quartz (qtz)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0    % olivine (olv)
               0  0  1  0  0  0  0  0    % oxide (oxd)
               0  0  0  1  1  0  0  0    % pyroxene (pxn)
               0  0  0  0  0  1  1  0    % plagioclase (plg)
               0  0  0  0  0  0  0  1];  % quartz (qtz)

% mineral end-member composition of melting model components
%               for     fay   ilm    opx    cpx    ant    alb    qtz
cal.cmp_mem = [ 99.48   0.5   0.01   0.01    0.0    0.0   0.0   0.0      % cphs0 => dunite (dun)
                 8.0    7.0   1.0   43.0    11.0   29.0   1.0   0.0      % perCx => gabbro (gbr)
                 0.5    4.0   5.0   17.0    37.0   20.0  16.5   0.0      % perCm => basalt (bas)
                 0.0    0.0   9.0    0.0    42.0   12.0  22.0  15.0];    % cphs1 => eutectic (eut)
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

% specify melting model phase diagram parameters
cal.cphs0    =  cal.cmp_oxd(1,1)/100; % phase diagram lower bound composition [wt SiO2]
cal.cphs1    =  cal.cmp_oxd(4,1)/100; % phase diagram upper bound composition [wt SiO2]
cal.Tphs0    =   950;                 % phase diagram lower bound temperature [degC]
cal.Tphs1    =  1890;                 % phase diagram upper bound temperature [degC]
cal.PhDg     =  [8.0,3.5,1.3,1.3];    % phase diagram curvature factor (> 1)
cal.perCm    =  cal.cmp_oxd(3,1)/100; % peritectic liquidus composition [wt SiO2]
cal.perCx    =  cal.cmp_oxd(2,1)/100; % peritectic solidus  composition [wt SiO2]
cal.perT     =  1125;                 % peritectic temperature [degC]
cal.clap     =  1e-7;                 % Clapeyron slope for P-dependence of melting T [degC/Pa]
cal.dTH2O    =  [1400,1300,1200];     % solidus shift from water content [degC/wt^0.75]

% specify geochemical model parameters
cal.nte      =  4;           % number of trace elements
cal.nir      =  2;           % number of isotope ratios
cal.Kte_mem  =  [0.01;0.10;3.00;10.0].*ones(cal.nte,cal.nmem);

% specify density parameters
cal.rhox0 = [3270,4390,7500,3200,3400,2730,2620,2650];                     % mineral end-member reference densities [kg/m3]
cal.rhof0 = 1000;                                                          % fluid reference density [kg/m3]
cal.aT    =  4e-5;                                                         % thermal expansivity [1/K]
cal.gH    =  0.75;                                                         % hydrous melt expansivity [1/(wt H2O)]
cal.bP    =  1e-8;                                                         % fluid compressibility [1/Pa]

% specify mixture viscosity parameters (Costa et al., 2009)
cal.Bphi    = 2.0;              % Einstein-Roscoe powerlaw coefficient bubbles
cal.Bchi    = 2.0;              % Einstein-Roscoe powerlaw coefficient crystals
cal.chi_pck = 0.60;             % rheologically critical crystal fraction
cal.gamma   = 2.50;             % step-function steepness coefficient
cal.delta   = 27;               % solid viscosity melt-weakening slope
cal.xi      = 4.5e-4;           % solid viscosity level
cal.etaf0   = 0.1;              % fluid viscosity constant

% specify segregation coefficient parameters
cal.bm      = 100;              % melt permeability geometric factor (k0 = dx^2/bm)
cal.cm      = 0.001;            % melt percolation threshold
cal.nm      = 3;                % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.mm      = 2;                % melt permeability powerlaw (k0*(mu-cm)^nm*(1-mu)^mm)
cal.bf      = 100;              % fluid permeability geometric factor (k0 = dx^2/bm)
cal.cf      = 0.1;              % fluid percolation threshold
cal.nf      = 4;                % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
cal.mf      = 2;                % fluid permeability powerlaw (k0*(phi-cf)^nf*(1-phi)^mf)
