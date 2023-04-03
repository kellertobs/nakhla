function  [var,cal,flag]  =  meltmodel(var,cal,type)

%*****  compute partition coefficients at P,T *****************************

if strcmp(type,'K')
    
    var.H2Om  =  max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));
    cal       =  K(var.T,var.P,var.H2Om,cal);
    
%*****  compute Tsol, Tliq at P,C^i  **************************************

elseif strcmp(type,'T')
    
    [cal,flag.Tsol]  =  Tsolidus (var,cal);
    [cal,flag.Tliq]  =  Tliquidus(var,cal);

%*****  compute equilibrium f, C_s^i, C_l^i fraction at P,T,C^i  **********

elseif strcmp(type,'E')
    
    [cal,var,flag]  =  equilibrium(var,cal);
    
end

end


function  [cal, flag]  =  Tsolidus(var,cal)

%*****  subroutine to compute solidus temperature at given bulk composition

%***  exclude invalid compositions
ii  =  sum(var.c(:,sum(var.c,1)>1),2)<=1;

%***  get T,P-,H2O-dependent partition coefficients Kci
cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  set starting guess for Tsol
if ~isfield(cal,'Tsol')
    Tsol =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c.*cal.Tm,2)));
else
    Tsol =  cal.Tsol;
end

%***  get T,P-,H2O-dependent partition coefficients Kci
var.T   = Tsol;
cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');


%***  get residual for sum(ci_b/Kci) = sum(ci_v)
r    =  sum(var.c./cal.Kc,2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-10; % tolerance for Newton residual
its_tol    =  1e3;   % maximum number of iterations
eps_T      =  1e-3;  % temperature perturbation for finite differencing, degrees
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Kci at T+eps_T/2
    var.T   = Tsol+eps_T/2;
    cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  get residual at T+eps_T/2
    rp   =  sum(var.c./cal.Kc,2)-1;
    
    %***  compute partition coefficients Kci at T-eps_T/2
    var.T   = Tsol-eps_T/2;
    cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  get residual at T-eps_T/2
    rm   =  sum(var.c./cal.Kc,2)-1;
    
    %***  finite difference drdT = (r(T+eps_T/2)-r(T-eps_T/2))/eps_T
    drdT  =  (rp - rm)/eps_T;
    
    %***  apply Newton correction to current guess of Tsol
    Tsol(ii)  =  Tsol(ii) - r(ii)./2./drdT(ii);
    
    %***  compute partition coefficients Kci at Tsol
    var.T   = Tsol;
    cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  compute residual at Tsol
    r  =  sum(var.c./cal.Kc,2)-1;
    
    %***  compute Newton residual norm
    rnorm  =  norm(r(ii),2)./sqrt(length(r(ii)));
    
    n  =  n+1;   %  update iteration count
    
    if (n==its_tol)
        %error(['!!! Newton solver for solidus T has not converged after ',num2str(its_tol),' iterations !!!']);
        flag = 0; break;
    end
    
end

cal.Tsol = Tsol;

end


function  [cal, flag]  =  Tliquidus(var,cal)

%*****  subroutine to compute liquidus temperature at given bulk composition

%***  exclude invalid compositions
ii  =  sum(var.c(:,sum(var.c,1)>1),2)<=1;

%***  get T,P-,H2O-dependent partition coefficients Kci
cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  set starting guess for Tliq
if ~isfield(cal,'Tliq')
    Tliq  =  max(min(min(cal.Tm)),min(max(max(cal.Tm)),sum(var.c.*cal.Tm,2)));
else
    Tliq =  cal.Tliq;
end

%***  get T,P-,H2O-dependent partition coefficients Kci
var.T   = Tliq;
cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  get residual for sum(ci_b*Kci) = sum(ci_v)
r  =  sum(var.c.*cal.Kc,2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-10; % tolerance for Newton residual
its_tol    =  1e3;   % maximum number of iterations
eps_T      =  1e-3;  % temperature perturbation for finite differencing, degrees
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Kci at T+eps_T/2
    var.T   = Tliq+eps_T/2;
    cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');

    %***  get residual at T+eps_T/2
    rp  =  sum(var.c.*cal.Kc,2)-1;
    
    %***  compute partition coefficients Kci at T-eps_T/2
    var.T   = Tliq-eps_T/2;
    cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  get residual at T-eps_T/2
    rm  =  sum(var.c.*cal.Kc,2)-1;
    
    %***  compute difference drdT = (r(T+eps_T/2)-r(T-eps_T/2))/eps_T
    drdT  =  (rp - rm)./eps_T;
    
    %***  apply Newton correction to Tliq
    Tliq(ii)  =  Tliq(ii) - r(ii)./2./drdT(ii);
    
    %***  compute partition coefficients Kci at Tliq
    var.T   = Tliq;
    cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');

    %***  compute residual at Tliq
    r  =  sum(var.c.*cal.Kc,2)-1;
    
    %***  compute Newton residual norm
    rnorm  =  norm(r(ii),2)./sqrt(length(r(ii)));
    
    n  =  n+1;   %  update iteration count
    
    if (n==its_tol)
        %error(['!!! Newton solver for liquidus T has not converged after ',num2str(its_tol),' iterations !!!']);
        flag = 0; break;
    end
    
end

cal.Tliq  =  Tliq;

end


function  [cal,var,flag]  =  equilibrium(var,cal)

%*****  subroutine to compute equilibrium melt fraction and phase 
%       compositions at given bulk composition, pressure and temperature

%***  get P-, H2O-dependent partition coefficients
cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  get Tsol, Tliq at c
[var,cal,flag]  =  meltmodel(var,cal,'T');

var.T = max(cal.Tsol,min(cal.Tliq,var.T));

%***  set reasonable initial guess for melt fraction
if ~isfield(var,'m') 
    var.m = max(2e-6,min(1-2e-6,(var.T-cal.Tsol)./(cal.Tliq-cal.Tsol)));
elseif ~any(var.m(:)>1e-9 & var.m(:)<1-1e-9)
    var.m = max(2e-6,min(1-2e-6,(var.T-cal.Tsol)./(cal.Tliq-cal.Tsol)));
end

%***  get P-, H2O-dependent partition coefficients
cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  set reasonable initial guess for bubble fraction
if ~isfield(var,'f') 
    var.f = max(0,min(1-var.m,var.H2O - var.m.*var.H2Om));
elseif ~any(var.f(:)>1e-9)
    var.f = max(0,min(1-var.m,var.H2O - var.m.*var.H2Om));
end

var.x = 1-var.m-var.f;
    
%***  compute residual of unity sum of components
r   =  sum(var.c./(var.m + var.x.*cal.Kc),2) - sum(var.c./(var.m./cal.Kc + var.x),2);

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-10;  % tolerance for Newton residual
its_tol    =  1e3;   % maximum number of iterations
eps_m      =  1e-6;  % temperature perturbation for finite differencing, degrees
flag.eql   =  1;     % tells us whether the Newton solver converged

sol = var.T <= cal.Tsol;
liq = var.T >= cal.Tliq;

while rnorm > rnorm_tol     % Newton iteration

    %***  get T,P-,H2O-dependent partition coefficients Kci
    varp = var;  varm = var;
    varp.m = var.m+eps_m/2;
    varm.m = var.m-eps_m/2;

    [varp,calp] = meltmodel(varp,cal,'K');
    [varm,calm] = meltmodel(varm,cal,'K');

    varp.f = (var.H2O - varp.m.*varp.H2Om)./(1 - varp.m.*varp.H2Om);
    varp.x = 1-varp.m-varp.f;
    varm.f = (var.H2O - varm.m.*varm.H2Om)./(1 - varm.m.*varm.H2Om);
    varm.x = 1-varm.m-varm.f;

    rp  =  sum(var.c./(varp.m + varp.x.*calp.Kc),2) - sum(var.c./(varp.m./calp.Kc + varp.x),2);
    rm  =  sum(var.c./(varm.m + varm.x.*calm.Kc),2) - sum(var.c./(varm.m./calm.Kc + varm.x),2);

    %***  compute analytic derivative of residual dr/dm
%     dr_dm  =  sum(-(var.c.*(cal.Kc - 1).^2)./(var.x.*cal.Kc + var.m).^2,2);
    dr_dm = (rp-rm)./(varp.m-varm.m);

    %***  apply Newton correction to melt fraction m
    var.m = max(0,min(1,var.m - r./dr_dm/4));
    var.m(liq) = 1;
    var.m(sol) = 0;

    %***  get crystal fraction x
    var.f        = max(0,min(var.H2O, (var.H2O - var.m.*var.H2Om)./(1 - var.m.*var.H2Om) ));
    unsat        = var.H2Om<cal.H2Osat;
    var.f(unsat) = 0;
    var.f(liq)   = max(0,(var.H2O(liq) - var.H2Om(liq))./(1 - var.H2Om(liq)));
    var.f(sol)   = max(0, var.H2O(sol));

    %***  get crystal fraction x
    var.x = 1-var.m;

    %***  compute residual of unity sum of components
    [var,cal] = meltmodel(var,cal,'K');
    r   =  sum(var.c./(var.m + var.x.*cal.Kc),2) - sum(var.c./(var.m./cal.Kc + var.x),2);
    r(liq | sol) = 0;

    %***  get non-linear residual norm
    rnorm  = norm(r(var.m>0 & var.x>0)./dr_dm(var.m>0 & var.x>0),2)./sqrt(length(r(var.m>0 & var.x>0))+1);
    
    n  =  n+1;  % update iteration count
    if (n==its_tol)
        disp(['!!! Newton equilibrium solver did not converge after ',num2str(its_tol),' iterations !!!']);
        flag.eql = 0; break;
    end
end

%***  get cxi, cmi as functions of Kci, ci and x, m, safeguard bounds
var.cm  =  max(0,min(1, var.c./(var.m + var.x.*cal.Kc) ));
var.cx  =  max(0,min(1, var.c./(var.m./cal.Kc + var.x) ));

var.m   =  var.m.*(1-var.f);
var.x   =  var.x.*(1-var.f);

end


function  [cal]  =  K(T,P,H2Om,cal)

%***  compute P-dependence of component melting points     
%     Parameterization as in Rudge etal (2011)
cal.Tm  =  (cal.T0 - cal.dTH2O.*H2Om.^cal.pH2O) .* (1 + P./cal.A) .^ (1./cal.B);

%***  compute T,P-dependence of equilibrium partition coefficients
%     Parameterization after Rudge, Bercovici, & Spiegelman (2010)

cal.L  = (T+273.15).*cal.dS;
cal.Kc  =  exp(cal.L./cal.r.*(1./(T+273.15) - 1./(cal.Tm+273.15)));
        
end % function


function  [cal]  =  H2Osat(T,P,SiO2m,cal)

c = (SiO2m.*100-min(cal.cmp_oxd(:,1)))./(max(cal.cmp_oxd(:,1))-min(cal.cmp_oxd(:,1)));
P = P*1e9;

H2Osat_maf = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T+273.15)-1./(1473.15))); % Katz et al., 2003; Moore et al., 1998
H2Osat_fls = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015
cal.H2Osat = (1-c).*H2Osat_maf + c.*H2Osat_fls;

end