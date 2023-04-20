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
var.H2Om  = cal.H2Osat .* (var.H2O>0);
[var,cal] = meltmodel(var,cal,'K');

%***  set starting guess for Tsol
if ~isfield(cal,'Tsol')
    Tsol =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c(:,1:cal.ncmp-1).*cal.Tm,2)));
else
    Tsol =  cal.Tsol;
end

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.T      = Tsol;
[var,cal]  = meltmodel(var,cal,'K');

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
    var.T      = Tsol+eps_T/2;
    [var,cal]  = meltmodel(var,cal,'K');
    
    %***  get residual at T+eps_T/2
    rp  =  sum((var.c-f.*cf)./((1-f).*cal.Kx + 1e-16),2)-1;
    
    %***  compute partition coefficients Kxi at T-eps_T/2
    var.T      = Tsol-eps_T/2;
    [var,cal]  = meltmodel(var,cal,'K');
    
    %***  get residual at T-eps_T/2
    rm  =  sum((var.c-f.*cf)./((1-f).*cal.Kx + 1e-16),2)-1;
    
    %***  finite difference drdT = (r(T+eps_T/2)-r(T-eps_T/2))/eps_T
    drdT  =  (rp - rm)/eps_T;
    
    %***  apply Newton correction to current guess of Tsol
    Tsol(ii)  =  max(min(cal.Tm(:)),Tsol(ii) - r(ii)./drdT(ii)/2);
    
    %***  compute partition coefficients Kxi at Tsol
    var.T      = Tsol;
    [var,cal]  = meltmodel(var,cal,'K');
    
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
var.H2Om   = min(var.H2O,cal.H2Osat);
[var,cal]  = meltmodel(var,cal,'K');

%***  set starting guess for Tliq
if ~isfield(cal,'Tliq')
    Tliq =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c(:,1:cal.ncmp-1).*cal.Tm,2)));
else
    Tliq =  cal.Tliq;
end

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.T      = Tliq;
[var,cal]  = meltmodel(var,cal,'K');

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
    var.T      = Tliq+eps_T/2;
    [var,cal]  = meltmodel(var,cal,'K');
    f          = (var.H2O-var.H2Om)./(1-var.H2Om);

    %***  get residual at T+eps_T/2
    rp  =  sum((var.c-f.*cf).*cal.Kx./(1-f),2)-1;
    
    %***  compute partition coefficients Kxi at T-eps_T/2
    var.T      = Tliq-eps_T/2;
    [var,cal]  = meltmodel(var,cal,'K');
    f          = (var.H2O-var.H2Om)./(1-var.H2Om);

    %***  get residual at T-eps_T/2
    rm  =  sum((var.c-f.*cf).*cal.Kx./(1-f),2)-1;
    
    %***  compute difference drdT = (r(T+eps_T/2)-r(T-eps_T/2))/eps_T
    drdT  =  (rp - rm)./eps_T;
    
    %***  apply Newton correction to Tliq
    Tliq(ii)  =  Tliq(ii) - r(ii)./drdT(ii)/2;
    
    %***  compute partition coefficients Kxi at Tliq
    var.T      = Tliq;
    [var,cal]  = meltmodel(var,cal,'K');
    f          = (var.H2O-var.H2Om)./(1-var.H2Om);

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
var.H2Om   = max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));
[var,cal]  = meltmodel(var,cal,'K');

%***  get Tsol, Tliq at c
[var,cal,flag]  =  meltmodel(var,cal,'T');

var.T = max(cal.Tsol,min(cal.Tliq,var.T));

%***  set reasonable initial guess for melt fraction
if ~isfield(var,'m') 
    var.m = max(0,min(1,(var.T-cal.Tsol)./(cal.Tliq-cal.Tsol)));
elseif ~any(var.m(:)>1e-9 & var.m(:)<1-1e-9)
    var.m = max(0,min(1,(var.T-cal.Tsol)./(cal.Tliq-cal.Tsol)));
end

%***  set reasonable initial guess for bubble fraction
if ~isfield(var,'f') 
    var.f = max(0,min(1-var.m,var.H2O - var.m.*var.H2Om));
elseif ~any(var.f(:)>1e-9)
    var.f = max(0,min(1-var.m,var.H2O - var.m.*var.H2Om));
end

var.f = max(0,min(var.H2O, var.H2O - var.m.*var.H2Om ));
var.m = min(1-var.f,var.m);
var.x = 1-var.m-var.f;
    
%***  compute residual of unity sum of components
var.H2Om  = max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));
[var,cal] = meltmodel(var,cal,'K');
var.cf    = zeros(size(var.c)); var.cf(:,end) = 1;
var.cm    = var.c./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16);
var.cx    = (var.c-var.f.*var.cf).*cal.Kx./(var.m + var.x.*cal.Kx + 1e-16);
rm        = sum(var.cm,2) - sum(var.cx,2);

rnorm     = 1;     % initialize residual norm for iterations
n         = 0;     % initialize iteration count
rnorm_tol = 1e-9;  % tolerance for Newton residual
its_tol   = 200;   % maximum number of iterations
eps_m     = 1e-9;  % temperature perturbation for finite differencing, degrees
flag.eql  = 1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol     % Newton iteration

    %***  compute numerical derivative of residual dr/dm
    varp        = var;
    varp.m      = var.m+eps_m/2;
    varp.H2Om   = max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));
    [varp,calp] = meltmodel(varp,cal,'K');
    varp.f      = varp.H2O - varp.m.*varp.H2Om;
    varp.x      = 1-varp.m-varp.f;
    varp.cm     =  varp.c./(varp.m + varp.x.*calp.Kx + varp.f.*calp.Kf + 1e-16);
    varp.cx     = (varp.c-varp.f.*varp.cf).*calp.Kx./(varp.m + varp.x.*calp.Kx + 1e-16);
    rmp         = sum(varp.cm,2) - sum(varp.cx,2);

    varm        = var;
    varm.m      = var.m-eps_m/2;
    varm.H2Om   = max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));
    [varm,calm] = meltmodel(varm,cal,'K');
    varm.f      = varm.H2O - varm.m.*varm.H2Om;
    varm.x      = 1-varm.m-varm.f;
    varm.cm     =  varm.c./(varm.m + varm.x.*calm.Kx + varm.f.*calm.Kf + 1e-16);
    varm.cx     = (varm.c-varm.f.*varm.cf).*calm.Kx./(varm.m + varm.x.*calm.Kx + 1e-16);
    rmm         = sum(varm.cm,2) - sum(varm.cx,2);

    drm_dm      = (rmp-rmm)./(varp.m-varm.m);

    %***  apply Newton correction to crystal fraction x
    var.m = max(0,min(1-var.f, var.m - rm./drm_dm/3 ));

    %***  get droplet fraction f
    var.f = max(0,min(var.H2O, var.H2O - var.m.*var.H2Om ));
    
    %***  get crystal fraction x
    var.x = max(0,min(1, 1-var.m-var.f ));

    %***  compute residual of unity sum of components
    var.H2Om  = max(0,min(1,min(var.H2O./(var.m+1e-16),cal.H2Osat)));
    [var,cal] = meltmodel(var,cal,'K');

    var.cm  =  var.c./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16);
    var.cx  = (var.c-var.f.*var.cf).*cal.Kx./(var.m + var.x.*cal.Kx + 1e-16);

    rm  = sum(var.cm,2) - sum(var.cx,2);
    rm((var.m==0 & var.T<cal.Tsol) | (var.x==0 & var.T>cal.Tliq)) = 0;
    
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

%***  compute P-dependence of component melting points     
%     Parameterization as in Rudge etal (2011)
cal.Tm  =  (cal.T0 - cal.dTH2O.*var.H2Om.^cal.pH2O) .* (1 + var.P./cal.A).^(1./cal.B);

%***  compute T,P-dependence of major component equilibrium partition coefficients
%     Parameterization after Rudge, Bercovici, & Spiegelman (2010)

cal.L   = (var.T+273.15).*cal.dS;

cal.Kx = zeros(size(var.c));
cal.Kx(:,1:end-1) = exp(cal.L./cal.r.*(1./(var.T+273.15) - 1./(cal.Tm+273.15)));

%***  compute volatile component equilibrium partition coefficient

cal.Kf = zeros(size(var.c));
cal.Kf(:,cal.ncmp) = 1./(cal.H2Osat+1e-16);
        
end % function