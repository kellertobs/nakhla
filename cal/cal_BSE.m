% specify petrological model composition parameters
% compositions specified as oxides (oxd) -> mineral endmembers (mem) 
% -> mineral systems (msy) -> melting model components (cmp)

% number of oxides, mineral end-members, mineral systems, model components
cal.noxd   = 7;
cal.nmem   = 12;
cal.nmsy   = 5;
cal.ncmp   = 5;

% label strings for all compositional representations
cal.oxdStr = {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','H$_2$O'};
     elStr = {'Si','Al','Fe','Mg','Ca','Na','H'};
cal.memStr = {'for','fay','ens','hyp','fsl','aug','pig','ant','alb','ilm','qtz','wat'};
cal.msyStr = {'olv','opx','cpx','fsp','oxs','qtz'};
cal.cmpStr = {'dun','pxn','bas','eut','fld'};

for i = 1:cal.ncmp; cal.(cal.cmpStr{i}) = i; end
for i = 1:cal.nmsy; cal.(cal.msyStr{i}) = i; end
for i = 1:cal.nmem; cal.(cal.memStr{i}) = i; end
for i = 1:cal.noxd; cal.(elStr{i}) = i; end


% oxide composition of mineral end-members
%                   Si      Al      Fe     Mg       Ca      Na     H2O
cal.mem_oxd    = [  42.7    0.0     0.0    57.3     0.0     0.0    0.0     % forsterite (for)
                    29.5    0.0    70.5     0.0     0.0     0.0    0.0     % fayalite (fay)
                    56.22   2.50    7.88   33.07    0.33    0.0    0.0     % enstatite (ens)
                    53.60   3.26   15.62   26.44    1.08    0.0    0.0     % hypersthene (hyp)
                    47.91   1.36   39.74    7.03    3.96    0.0    0.0     % ferrosillite (fsl)
                    52.40   1.00   15.55   13.16   17.66    0.23   0.0     % augite (aug)
                    49.64   1.05   30.66    4.11   13.65    0.89   0.0     % pigeonite (pig)
                    44.4    35.8    0.0     0.0    19.3     0.5    0.0     % anorthite (ant)
                    67.3    20.2    0.0     0.0     0.5    11.5    0.0     % albite (alb)
                     0.0     0.0  100.0     0.0     0.5    11.5    0.0     % ilmenite (ilm)
                   100.0     0.0    0.0     0.0     0.0     0.0    0.0     % quartz (qtz)
                     0.0     0.0    0.0     0.0     0.0     0.0  100.0];   % water (wat)
cal.mem_oxd = cal.mem_oxd./sum(cal.mem_oxd,2)*100;

% mineral end-members in mineral systems
cal.msy_mem = [1  1  0  0  0  0  0  0  0  0  0  0    % olivine (olv)
               0  0  1  1  1  0  0  0  0  0  0  0    % orthopyroxene (opx)
               0  0  0  0  0  1  1  0  0  0  0  0    % clinopyroxene (cpx)
               0  0  0  0  0  0  0  1  1  0  0  0    % feldspar (fsp)
               0  0  0  0  0  0  0  0  0  1  0  0    % oxides (oxs)
               0  0  0  0  0  0  0  0  0  0  1  0];  % quartz (qtz)

% mineral end-member composition of melting model components
%                   for           fay        ens       hyp      fsl       aug       pig       ant       alb       ilm       qtz       wat
cal.cmp_mem =   [  100.0000         0         0         0         0         0         0         0         0         0         0         0
                    33.5232   23.9360   41.6890    0.8517    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000    0.0000         0
                     0.0000    0.0000    0.0538   65.0236    0.8601    0.4994    0.2419   33.1366    0.1847    0.0000    0.0000         0
                     0.0000    0.0000    0.0000    0.1203    1.0667    0.0000   27.1951   11.8682    5.0786    0.0000   54.6712         0
                          0         0         0         0         0         0         0         0         0         0         0  100.0000];
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
cal.T0(cal.pxn) =  1450;
cal.T0(cal.bas) =  1170;
cal.T0(cal.eut) =   895;

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
cal.r(cal.dun)  =  25.0;
cal.r(cal.pxn)  =  12.8;
cal.r(cal.bas)  =  19.4;
cal.r(cal.eut)  =   9.6;

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
