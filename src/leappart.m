%*****  Local Equilibrium Approximation using Pseudo-component PARTitioning

function  [var,cal,flag]  =  leappart(var,cal,type)

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
ii  =  sum(var.c,2)<=1+1e-15;

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.H2Om  = cal.H2Osat .* (var.H2O>0);
[var,cal] = leappart(var,cal,'K');

%***  set starting guess for Tsol
if ~isfield(cal,'Tsol')
    Tsol =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c(:,1:cal.ncmp-1).*cal.Tm,2)));
else
    Tsol =  cal.Tsol;
end

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.T      = Tsol;
[var,cal]  = leappart(var,cal,'K');

%***  get fluid fraction and composition
f  = var.H2O;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
its_tol    =  100;   % maximum number of iterations
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > cal.tol  % iterate down to full accuracy

    %***  get partition coefficients Kxi
    var.T      = Tsol;
    [var,cal]  = leappart(var,cal,'K');
    TK   = Tsol+273.15;
    TmK  = cal.Tm+273.15;

    %***  get residual r
    r    = sum( [-(var.c(:,1:end-1).*exp((cal.Dsx.*(TK - TmK))./(TmK.*cal.r)))./(f - 1),var.H2Om] ,2) - 1;

    %***  get analytical derivative of residual dr/dT
    drdT = sum( -(cal.Dsx.*var.c(:,1:end-1).*exp((cal.Dsx.*(TK - TmK))./(TmK.*cal.r)))./(TmK.*cal.r.*(f - 1)) ,2);
    
    %***  apply Newton correction to current guess of Tsol
    Tsol(ii)  =  Tsol(ii) - r(ii)./drdT(ii)/1.1;

    %***  compute Newton residual norm
    rnorm  =  norm(r(ii),2)./sqrt(length(r(ii)));
    
    n  =  n+1;   %  update iteration count
    
    if (n==its_tol)
        disp(['!!! Newton solver for solidus T not converged after ',num2str(its_tol),' iterations !!!']);
        flag = 0; break;
    end
    
end

cal.Tsol = Tsol;

end


function  [cal, flag]  =  Tliquidus(var,cal)

%*****  subroutine to compute liquidus temperature at given bulk composition

%***  exclude invalid compositions
ii  =  sum(var.c,2)<=1+1e-15;

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.H2Om   = min(var.H2O,cal.H2Osat);
[var,cal]  = leappart(var,cal,'K');

%***  set starting guess for Tliq
if ~isfield(cal,'Tliq')
    Tliq =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c(:,1:cal.ncmp-1).*cal.Tm,2)));
else
    Tliq =  cal.Tliq;
end

%***  get T,P-,H2O-dependent partition coefficients Kxi
var.T      = Tliq;
[var,cal]  = leappart(var,cal,'K');

%***  get fluid fraction
f  = (var.H2O-var.H2Om)./(1-var.H2Om);

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
its_tol    =  100;   % maximum number of iterations
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > cal.tol  % iterate down to full accuracy

    %***  get partition coefficients Kxi
    var.T      = Tliq;
    [var,cal]  = leappart(var,cal,'K');
    TK   = Tliq+273.15;
    TmK  = cal.Tm+273.15;

    %***  get residual r
    r    = sum( -(var.c(:,1:end-1).*exp(-(cal.Dsx.*(TK - TmK))./(TmK.*cal.r)))./(f - 1) ,2) - 1;

    %***  get get analytical derivative of residual dr/dT
    drdT = sum( (cal.Dsx.*var.c(:,1:end-1).*exp(-(cal.Dsx.*(TK - TmK))./(TmK.*cal.r)))./(TmK.*cal.r.*(f - 1)) ,2);

    %***  apply Newton correction to Tliq
    Tliq(ii)  =  Tliq(ii) - r(ii)./drdT(ii)/1.1;
    
    %***  compute Newton residual norm
    rnorm  =  norm(r(ii),2)./sqrt(length(r(ii)));
    
    n  =  n+1;   %  update iteration count
    
    if (n==its_tol)
        disp(['!!! Newton solver for liquidus T not converged after ',num2str(its_tol),' iterations !!!']);
        flag = 0; break;
    end
    
end

cal.Tliq  =  Tliq;

end


function  [cal,var,flag]  =  equilibrium(var,cal)

%*****  subroutine to compute equilibrium melt fraction and phase 
%       compositions at given bulk composition, pressure and temperature

%***  get Tsol, Tliq at c
var.H2Om        = max(0,min(cal.H2Osat, var.H2O./(var.m+1e-16) ));
[var,cal,flag]  = leappart(var,cal,'T');

if isfield(cal,'Tsol')
    var.T = max(cal.Tsol,min(cal.Tliq,var.T));
end

var.H2Om = max(0,min(cal.H2Osat, var.H2O./(var.m+1e-16) ));
var.cf   = zeros(size(var.c)); var.cf(:,cal.ncmp) = 1;
var.f    = max(0,min(var.H2O, var.H2O - cal.H2Osat ));
var.m    = max(0,min(1-var.f,var.m));
var.x    = max(0,min(1,1-var.m-var.f));

rnorm     = 1;     % initialize residual norm for iterations
n         = 0;     % initialize iteration count
its_tol   = 100;   % maximum number of iterations
flag.eql  = 1;     % tells us whether the Newton solver converged
eps       = 1e-9;

while rnorm > cal.tol     % Newton iteration

    %***  get residual of unity sum of components
    [var,cal]   = leappart(var,cal,'K');
    r           = sum( (var.c.*(1-cal.Kx)   )./(cal.Kx + (1-cal.Kx).*var.m + (cal.Kf-cal.Kx).*var.f + 1e-32)   ,2);

    % %***  get analytical derivative of residual dr/dm
    % drdm_a    = sum(-(var.c.*(1-cal.Kx).^2)./(cal.Kx + (1-cal.Kx).*var.m + (cal.Kf-cal.Kx).*var.f + 1e-32).^2,2);

    %*** get numerical derivative of residual dr/dm
    varp        = var;
    varp.m      = var.m+eps;
    varp.H2Om   = max(0,min(cal.H2Osat, varp.H2O./(varp.m+1e-32) ));
    varp.f      = max(0,min(varp.H2O, varp.H2O - varp.m.*cal.H2Osat));
    [varp,calp] = leappart(varp,cal,'K');
    rp          = sum( (varp.c.*(1-calp.Kx)   )./(calp.Kx + (1-calp.Kx).*varp.m + (calp.Kf-calp.Kx).*varp.f + 1e-32)   ,2);

    varm        = var;
    varm.m      = var.m-eps;
    varm.H2Om   = max(0,min(cal.H2Osat, varm.H2O./(varm.m+1e-32) ));
    varm.f      = max(0,min(varm.H2O, varm.H2O - varm.m.*cal.H2Osat));
    [varm,calm] = leappart(varm,cal,'K');
    rm          = sum( (varm.c.*(1-calm.Kx)   )./(calm.Kx + (1-calm.Kx).*varm.m + (calm.Kf-calm.Kx).*varm.f + 1e-32)   ,2);

    drdm_n      = (rp-rm)./2/eps;

    %***  apply Newton correction to crystal fraction x
    upd_m = - cal.alpha * r./drdm_n;
    upd_m = max(-var.m/2,min((1-var.m)/2, upd_m ));
    var.m = max(0,min(1, var.m + upd_m ));

    %***  get droplet fraction f
    var.H2Om = max(0,min(cal.H2Osat, var.H2O./(var.m+1e-32) ));
    var.f    = max(0,min(var.H2O, var.H2O - var.m.*cal.H2Osat));

    %***  get crystal fraction x
    var.x     = max(0,min(1, 1-var.m-var.f ));
    
    %***  get non-linear residual norm
    rnorm  = norm(upd_m,2)/norm(ones(size(upd_m)),2);

    n  =  n+1;  % update iteration count
    if (n==its_tol)
        disp(['!!! Newton equilibrium solver not converged after ',num2str(its_tol),' iterations !!!']);
        flag.eql = 0; break;
    end

end

var.cm = max(0,min(1, var.c        ./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16) ));
var.cx = max(0,min(1, var.c.*cal.Kx./(var.m + var.x.*cal.Kx + var.f.*cal.Kf + 1e-16) ));

end


function [r,var,cal] = resm(m,var,cal)

var.H2Om  = max(0,min(cal.H2Osat, var.H2O./(m+1e-16) ));
[var,cal] = leappart(var,cal,'K');
f         = max(0,min(var.H2O, var.H2O - m.*var.H2Om ));
x         = 1-m-f;
var.cm    = var.c        ./(m + x.*cal.Kx + f.*cal.Kf + 1e-16);
var.cx    = var.c.*cal.Kx./(m + x.*cal.Kx + f.*cal.Kf + 1e-16);
r         = sum(var.cm,2) - sum(var.cx,2);

end


function  [var,cal]  =  K(var,cal)

%***  compute P-dependence of component melting points     
%     Parameterization as in Rudge etal (2011)
cal.Tm  =  (cal.T0 - cal.dTH2O.*var.H2Om.^cal.pH2O) .* (1 + var.P./cal.A).^(1./cal.B);

%***  compute T,P-dependence of major component equilibrium partition coefficients
%     Parameterization after Rudge, Bercovici, & Spiegelman (2010)

cal.L   = (var.T+273.15).*cal.Dsx;

cal.Kx  = zeros(size(var.c));
cal.Kx(:,1:end-1) = exp(cal.L./cal.r.*(1./(var.T+273.15) - 1./(cal.Tm+273.15)));

%***  compute volatile component equilibrium partition coefficient

cal.Kf  = zeros(size(var.c));
cal.Kf(:,cal.ncmp) = 1./(cal.H2Osat+1e-16);
        
end % function