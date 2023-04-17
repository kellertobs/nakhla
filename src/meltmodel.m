function  [var,cal,flag]  =  meltmodel(var,cal,type)

%*****  compute partition coefficients at P,T *****************************

if strcmp(type,'K')
    
    [var,cal]  =  K(var,cal);
    
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

%***  get T,P-,H2O-dependent partition coefficients Kxi
cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  set starting guess for Tsol
if ~isfield(cal,'Tsol')
    Tsol =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c(:,1:cal.ncmp-1).*cal.Tm,2)));
else
    Tsol =  cal.Tsol;
end

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.T   = Tsol;
cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  get fluid fraction and composition
f  = var.H2O;
cf = [zeros(1,cal.ncmp-1),1];

%***  get residual for sum(ci_b/Kxi) = sum(ci_v)
r  =  sum((var.c-f.*cf)./((1-f).*cal.Kx + 1e-16),2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-9; % tolerance for Newton residual
its_tol    =  200;   % maximum number of iterations
eps_T      =  1e-3;  % temperature perturbation for finite differencing, degrees
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Kxi at T+eps_T/2
    var.T     = Tsol+eps_T/2;
    cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  get residual at T+eps_T/2
    rp  =  sum((var.c-f.*cf)./((1-f).*cal.Kx + 1e-16),2)-1;
    
    %***  compute partition coefficients Kxi at T-eps_T/2
    var.T     = Tsol-eps_T/2;
    cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  get residual at T-eps_T/2
    rm  =  sum((var.c-f.*cf)./((1-f).*cal.Kx + 1e-16),2)-1;
    
    %***  finite difference drdT = (r(T+eps_T/2)-r(T-eps_T/2))/eps_T
    drdT  =  (rp - rm)/eps_T;
    
    %***  apply Newton correction to current guess of Tsol
    Tsol(ii)  =  max(min(cal.Tm(:)),Tsol(ii) - r(ii)./drdT(ii)/2);
    
    %***  compute partition coefficients Kxi at Tsol
    var.T     = Tsol;
    cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    
    %***  compute residual at Tsol
    r  =  sum((var.c-f.*cf)./((1-f).*cal.Kx + 1e-16),2)-1;
    
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

%***  get T,P-,H2O-dependent partition coefficients Kxi
cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  set starting guess for Tliq
if ~isfield(cal,'Tliq')
    Tliq =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c(:,1:cal.ncmp-1).*cal.Tm,2)));
else
    Tliq =  cal.Tliq;
end

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.T   = Tliq;
cal     = H2Osat(var.T,var.P,var.SiO2m,cal);
[var,cal] = meltmodel(var,cal,'K');

%***  get fluid fraction and composition
f  = (var.H2O-var.H2Om)./(1-var.H2Om);
cf = [zeros(1,cal.ncmp-1),1];

%***  get residual for sum(ci_b*Kxi) = sum(ci_v)
r  =  sum((var.c-f.*cf).*cal.Kx./(1-f),2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-9; % tolerance for Newton residual
its_tol    =  200;   % maximum number of iterations
eps_T      =  1e-3;  % temperature perturbation for finite differencing, degrees
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Kxi at T+eps_T/2
    var.T     = Tliq+eps_T/2;
    cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    f         = (var.H2O-var.H2Om)./(1-var.H2Om);

    %***  get residual at T+eps_T/2
    rp  =  sum((var.c-f.*cf).*cal.Kx./(1-f),2)-1;
    
    %***  compute partition coefficients Kxi at T-eps_T/2
    var.T     = Tliq-eps_T/2;
    cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    f         = (var.H2O-var.H2Om)./(1-var.H2Om);

    %***  get residual at T-eps_T/2
    rm  =  sum((var.c-f.*cf).*cal.Kx./(1-f),2)-1;
    
    %***  compute difference drdT = (r(T+eps_T/2)-r(T-eps_T/2))/eps_T
    drdT  =  (rp - rm)./eps_T;
    
    %***  apply Newton correction to Tliq
    Tliq(ii)  =  Tliq(ii) - r(ii)./drdT(ii)/2;
    
    %***  compute partition coefficients Kxi at Tliq
    var.T     = Tliq;
    cal       = H2Osat(var.T,var.P,var.SiO2m,cal);
    [var,cal] = meltmodel(var,cal,'K');
    f         = (var.H2O-var.H2Om)./(1-var.H2Om);

    %***  compute residual at Tliq
    r  =  sum((var.c-f.*cf).*cal.Kx./(1-f),2)-1;
    
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
var.cf  = zeros(size(var.c)); var.cf(:,end) = 1;
var.cm  = max(0,min(1, var.c./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16) ));
var.cx  = max(0,min(1,(var.c-var.f.*var.cf).*cal.Kx./(var.m + var.x.*cal.Kx + 1e-16) ));
rm      = sum(var.cm,2) - sum(var.cx,2);

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-9;  % tolerance for Newton residual
its_tol    =  200;   % maximum number of iterations
eps_m      =  1e-6;  % temperature perturbation for finite differencing, degrees
flag.eql   =  1;     % tells us whether the Newton solver converged

sol = var.T <= cal.Tsol;
liq = var.T >= cal.Tliq;

while rnorm > rnorm_tol     % Newton iteration

    %***  get T,P-,H2O-dependent partition coefficients Kxi
    varp      = var;
    varp.m    = var.m+eps_m/2;
    [var,cal] = meltmodel(varp,cal,'K');
    var.f     = var.H2O - var.m.*var.H2Om;
    var.x     = 1-var.m-var.f;
    var.cm    = max(0,min(1, var.c./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16) ));
    var.cx    = max(0,min(1,(var.c-var.f.*var.cf).*cal.Kx./(var.m + var.x.*cal.Kx + 1e-16) ));

    rmp       = sum(var.cm,2) - sum(var.cx,2);

    varm      = var;
    varm.m    = var.m-eps_m/2;
    [var,cal] = meltmodel(varm,cal,'K');
    var.f     = var.H2O - var.m.*var.H2Om;
    var.x     = 1-var.m-var.f;
    var.cm    = max(0,min(1, var.c./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16) ));
    var.cx    = max(0,min(1,(var.c-var.f.*var.cf).*cal.Kx./(var.m + var.x.*cal.Kx + 1e-16) ));

    rmm       = sum(var.cm,2) - sum(var.cx,2);

    %***  compute analytic derivative of residual dr/dm
    drm_dm    = (rmp-rmm)./(varp.m-varm.m);

    %***  apply Newton correction to crystal fraction x
    var.m      = max(0,min(1, var.m - rm./drm_dm/2 ));
    var.m(liq) = 1-var.f(liq);
    var.m(sol) = 0;

    %***  get droplet fraction f
    [var,cal]    = meltmodel(var,cal,'K');
    var.f        = max(0,min(var.H2O, var.H2O - var.m.*var.H2Om ));
    unsat        = var.H2Om<cal.H2Osat;
    var.f(unsat) = 0;
    var.f(liq)   = max(0,(var.H2O(liq) - var.H2Om(liq))./(1 - var.H2Om(liq)));
    var.f(sol)   = max(0, var.H2O(sol));

    %***  get melt fraction m
    var.x = max(0,min(1,1-var.m-var.f));

    %***  compute residual of unity sum of components
    var.cm  = max(0,min(1, var.c./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16) ));
    var.cx  = max(0,min(1,(var.c-var.f.*var.cf).*cal.Kx./(var.m + var.x.*cal.Kx + 1e-16) ));
    rm      = sum(var.cm,2) - sum(var.cx,2);
    rm(sol | liq) = 0;
    
    %***  get non-linear residual norm
    rnorm  = norm(rm./drm_dm,2)./sqrt(length(rm)+1);

    n  =  n+1;  % update iteration count
    if (n==its_tol)
        disp(['!!! Newton equilibrium solver did not converge after ',num2str(its_tol),' iterations !!!']);
        flag.eql = 0; break;
    end
end

end


function  [var,cal]  =  K(var,cal)

%***  get water content in melt for solidus depression

var.H2Om  =  max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));

%***  compute P-dependence of component melting points     
%     Parameterization as in Rudge etal (2011)
cal.Tm  =  (cal.T0 - cal.dTH2O.*var.H2Om.^cal.pH2O) .* (1 + var.P./cal.A) .^ (1./cal.B);

%***  compute T,P-dependence of major component equilibrium partition coefficients
%     Parameterization after Rudge, Bercovici, & Spiegelman (2010)

cal.L   = (var.T+273.15).*cal.dS;
cal.Kx  = exp(cal.L./cal.r.*(1./(var.T+273.15) - 1./(cal.Tm+273.15)));

cal.Kx(:,cal.ncmp) = 0;

%***  compute volatile component equilibrium partition coefficient

cal.Kf = zeros(size(cal.Kx));
cal.Kf(:,cal.ncmp) = 1./(var.H2Om+1e-16);
        
end % function


function  [cal]  =  H2Osat(T,P,SiO2m,cal)

c = (SiO2m.*100-min(cal.cmp_oxd(:,1)))./(max(cal.cmp_oxd(:,1))-min(cal.cmp_oxd(:,1)));
P = P*1e9;

H2Osat_maf = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T+273.15)-1./(1473.15))); % Katz et al., 2003; Moore et al., 1998
H2Osat_fls = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015
cal.H2Osat = (1-c).*H2Osat_maf + c.*H2Osat_fls;

end