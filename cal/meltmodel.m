function  [var,cal]  =  meltmodel(var,cal,type)

%*****  compute partition coefficients at P,T *****************************

if strcmp(type,'K')
    
    cal  =  H2Osat(var.T,var.P,var.SiO2l,cal);
    cal  =  K(var.T,var.P,var.H2Ol,cal);
    var  =  [];
    
%*****  compute Tsol, Tliq at P,C^i  **************************************

elseif strcmp(type,'T')
    
    cal  =  Tsolidus( var,cal);
    cal  =  Tliquidus(var,cal);
    var  =  [];

%*****  compute equilibrium f, C_s^i, C_l^i fraction at P,T,C^i  **********

elseif strcmp(type,'E')
    
    [cal,var]  =  equilibrium(var,cal);
    
end

end


function  [cal, flag]  =  Tsolidus(var,cal)

%*****  subroutine to compute solidus temperature at given bulk composition

%***  exclude invalid compositions
ii   =  sum(var.c(:,sum(var.c,1)>1),2)<=1+1e-14;

%***  get H2O saturation limit in liquid phase
cal  =  H2Osat(var.T,var.P,var.SiO2l,cal);

%***  get P-dependent pure-component melting Temperature
cal  =  K(var.T,var.P,cal);

%***  set starting guess for Tsol
if ~isfield(cal,'Tsol')
    Tsol =  max(min(cal.Tm(:)),min(max(cal.Tm(:)),sum(var.c.*cal.Tm,2)));
else
    Tsol =  cal.Tsol;
end

%***  get T-dependent partition coefficients Ki
cal  =  K(Tsol,var.P,cal);

%***  get residual for sum(ci_b/Ki) = 1
r    =  sum(var.c./cal.K,2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-12; % tolerance for Newton residual
its_tol    =  1e3;   % maximum number of iterations
eps_T      =  0.1;   % temperature perturbation for finite differencing, degrees
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Ki at T+eps_T
    cal   =  K(Tsol+eps_T,var.P,cal);
    
    %***  get residual at T+eps_T
    rp   =  sum(var.c./cal.K,2)-1;
    
    %***  compute partition coefficients Ki at T-eps_T
    cal   =  K(Tsol-eps_T,var.P,cal);
    
    %***  get residual at T-eps_T
    rm   =  sum(var.c./cal.K,2)-1;
    
    %***  finite difference drdT = (r(T+eps_T)-r(T-eps_T))/2/eps_T
    drdT  =  (rp - rm)./2./eps_T;
    
    %***  apply Newton correction to current guess of Tsol
    Tsol(ii)  =  Tsol(ii) - r(ii)./2./drdT(ii);
    
    %***  compute partition coefficients Ki at Tsol
    cal   =  K(Tsol,var.P,cal);
    
    %***  compute residual at Tsol
    r  =  sum(var.c./cal.K,2)-1;
    
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
ii  =  sum(var.c(:,sum(var.c,1)>1),2)<=1+1e-14;

%***  get P-dependent pure component melting T
cal  =  K(var.T,var.P,cal);

%***  set starting guess for Tliq
if ~isfield(cal,'Tliq')
    Tliq  =  max(min(min(cal.Tm)),min(max(max(cal.Tm)),sum(var.c.*cal.Tm,2)));
else
    Tliq =  cal.Tliq;
end

%***  get T-dependent partition coefficients Ki
cal  =  K(Tliq,var.P,cal);

%***  get residual for sum(ci_b*Ki) = 1
r  =  sum(var.c.*cal.K,2)-1;

rnorm      =  1;     % initialize residual norm for iterations
n          =  0;     % initialize iteration count
rnorm_tol  =  1e-12; % tolerance for Newton residual
its_tol    =  1e3;   % maximum number of iterations
eps_T      =  0.1;   % temperature perturbation for finite differencing, degrees
flag       =  1;     % tells us whether the Newton solver converged

while rnorm > rnorm_tol  % iterate down to full accuracy
    
    %***  compute partition coefficients Ki at T+eps
    cal  =  K(Tliq+eps_T,var.P,cal);
    
    %***  get residual at T+eps_T
    rp  =  sum(var.c.*cal.K,2)-1;
    
    %***  compute partition coefficients Ki at T-eps_T
    cal  =  K(Tliq-eps_T,var.P,cal);
    
    %***  get residual at T-eps_T
    rm  =  sum(var.c.*cal.K,2)-1;
    
    %***  compute difference drdT = (r(T+eps_T)-r(T-eps_T))/2/eps_T
    drdT  =  (rp - rm)./2./eps_T;
    
    %***  apply Newton correction to Tliq
    Tliq(ii)  =  Tliq(ii) - r(ii)./2./drdT(ii);
    
    %***  compute partition coefficients Ki at Tliq
    cal  =  K(Tliq,var.P,cal);
    
    %***  compute residual at Tliq
    r  =  sum(var.c.*cal.K,2)-1;
    
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


function  [cal,var]  =  equilibrium(var,cal)

%*****  subroutine to compute equilibrium melt fraction and phase 
%       compositions at given bulk composition, pressure and temperature

%***  get Tsol, Tliq at c
[cal,flag]  =  Tsolidus( var,cal); if flag~=1, var.flag = flag; return; end
[cal,flag]  =  Tliquidus(var,cal); if flag~=1, var.flag = flag; return; end

%***  get H2O concentration in the liquid
[cal]    = H2Osat(var.T,var.P,var.SiO2l,cal);
var.H2Ol = max(0,min(var.H2O./var.f,var.H2Osat));

%***  get T,P-dependent partition coefficients Ki
cal  =  K(max(cal.Tsol,min(cal.Tliq,var.T)),var.P,cal);

%***  set reasonable initial condition if f does not exist or is zero
if ~isfield(var,'f') 
    var.f = max(0,min(1,(var.T-cal.Tsol)./(cal.Tliq-cal.Tsol)));
elseif ~any(var.f(:)>1e-9)
    var.f = max(0,min(1,(var.T-cal.Tsol)./(cal.Tliq-cal.Tsol)));
end
    
%***  compute residual of unity sum of components
ff   =  repmat(var.f,1,cal.ncmp);
r    =  sum(var.c./(ff+(1-ff).*cal.K),2) - sum(var.c./(ff./cal.K+(1-ff)),2);

rnorm     =  1;    % initialize residual norm for iterations
n         =  0;    % initialize iteration count
rnorm_tol = 1e-12; % tolerance for Newton residual
its_tol   = 1e3;   % maximum number of iterations
flag      =  1;    % tells us whether the Newton solver converged

var.f(var.T <= cal.Tsol) = 0;
var.f(var.T >= cal.Tliq) = 1;

ii = var.T > cal.Tsol & var.T < cal.Tliq;

while rnorm > rnorm_tol     % Newton iteration
    
    %***  compute analytic derivative of residual dr/df
    dr_df  =  - sum(var.c.*(1 -cal.K  )./(ff+(1-ff).*cal.K).^2,2) ...
              + sum(var.c.*(1./cal.K-1)./(ff./cal.K+(1-ff)).^2,2);
    
    %***  apply Newton correction to f
    var.f(ii) = max(0,min(1,var.f(ii) - r(ii)./2./dr_df(ii)));
    
    
    %***  compute residual of unity sum of components
    ff     =  repmat(var.f,1,cal.ncmp);
    r      =  sum(var.c./(ff+(1-ff).*cal.K),2) - sum(var.c./(ff./cal.K+(1-ff)),2);
    
    %***  get non-linear residual norm
    rnorm  =  norm(r(var.f>0 & var.f<1),2)./sqrt(length(r(var.f>0 & var.f<1))+1);
    
    n  =  n+1;  % update iteration count
    if (n==its_tol)
%         error(['!!! Newton solver for equilibrium f did not converge after ',num2str(its_tol),' iterations !!!']);
        flag = 0; break;
    end
end

%*** flag tells us whether Newton solver converged 
var.flag = flag;

%***  get C_s^i, C_l^i as functions of K^i, C^i and f, safeguard bounds
var.cl  =  max(0,min(1, var.c./(ff + (1-ff).*cal.K) ));
var.cs  =  max(0,min(1, var.c./(ff./cal.K + (1-ff)) ));

end


function  [cal]  =  K(T,P,H2Ol,cal)

%***  compute P-dependence of component melting points     
%     Parameterization as in Rudge etal (2011)
cal.Tm  =  (cal.T0 - cal.dTH2O.*H2Ol.^cal.pH2O) .* (1 + P./cal.A) .^ (1./cal.B);

%***  compute T,P-dependence of equilibrium partition coefficients
%     Parameterization after Rudge, Bercovici, & Spiegelman (2010)

cal.L  = (T+273.15).*cal.dS;
cal.K  =  exp(cal.L./cal.r.*(1./(T+273.15) - 1./(cal.Tm+273.15)));
        
end % function


function  [cal]  =  H2Osat(T,P,SiO2l,cal)

c = (SiO2l-min(cal.cmp_oxd(:,1)))./(max(cal.cmp_oxd(:,1))-min(cal.cmp_oxd(:,1)));

H2Osat_maf = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T+273.15)-1./(1200+273.15))); % Katz et al., 2003; Moore et al., 1998
H2Osat_fls = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015
cal.H2Osat = (1-c).*H2Osat_maf + c.*H2Osat_fls;

end