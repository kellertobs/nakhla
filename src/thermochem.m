%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy) density

% heat advection
advn_S = - advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % solid advection
         - advect(F.*sf,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

diff_S = diffus(T,ks,h,[1,2],BCD);

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1)+273.15)-T)./tau_T .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2)+273.15)-T)./tau_T .* botshape; end
if ~isnan(Twall(3)); bnd_T = bnd_T + ((Twall(3)+273.15)-T)./tau_T .* sdsshape; end
bnd_S = RHO.*cP.*bnd_T./T;

% total rate of change
dSdt  = advn_S + diff_S + diss_h + bnd_S;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo)/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% semi-implicit update of bulk entropy density
S = (a2*So+a3*Soo + (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt)/a1;

% convert entropy desnity to temperature
T = (min(cal.T0)+273.15)*exp((S - X.*Dsx - F.*Dsf)./RHO./cP + Adbt./cP.*(Pt-Ptop));


%***  update major component (SiO2) density

% major component advection
advn_C = - advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % solid advection
         - advect(F.*cf,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

% boundary layers
bnd_C = zeros(size(C));
if ~isnan(cwall(1)); bnd_C = bnd_C + (RHO.*cwall(1,:)-C).*mu./tau_a .* topshape; end
if ~isnan(cwall(2)); bnd_C = bnd_C + (RHO.*cwall(2,:)-C).*mu./tau_a .* botshape; end
if ~isnan(cwall(3)); bnd_C = bnd_C + (RHO.*cwall(3,:)-C).*mu./tau_a .* sdsshape; end

% total rate of change
dCdt = advn_C + bnd_C;                                            

% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
C = (a2*Co+a3*Coo + (b1*dCdt + b2*dCdto + b3*dCdtoo)*dt)/a1;

% apply minimum/maximum bounds
C = max(0, C );

% convert component density to concentration
c = C./sum(C,3);


%*** update phase equilibrium
eqtime = tic;

var.c     = reshape(c(:,:,1:end-1)./sum(c(:,:,1:end-1),3),Nx*Nz,cal.ncmp-1);   % component fractions [wt]
var.T     = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
var.P     = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
var.m     = reshape(mq,Nx*Nz,1);         % melt fraction [wt]
var.f     = reshape(fq,Nx*Nz,1);         % bubble fraction [wt]
var.H2O   = reshape(c(:,:,end),Nx*Nz,1); % water concentration [wt]
var.SiO2m = reshape(cm_oxd(:,:,1)./sum(cm_oxd(:,:,1:end-1),3),Nx*Nz,1); % melt silica concentration [wt]

[var,cal] =  meltmodel(var,cal,'E');

mq = reshape(var.m,Nz,Nx);
fq = reshape(var.f,Nz,Nx);
xq = reshape(var.x,Nz,Nx);

cxq = reshape([var.cx,zeros(size(var.T))    ],Nz,Nx,cal.ncmp);
cmq = reshape([var.cm.*(1-var.H2Om),var.H2Om],Nz,Nx,cal.ncmp);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;


%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_F   = - advect(F,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_F+advn_M;

% phase mass transfer rates
Gx = 2/3*Gx + 1/3*(xq.*RHO-X)./max(tau_r,4*dt);
Gf = 2/3*Gf + 1/3*(fq.*RHO-F)./max(tau_r,4*dt);
Gm = 2/3*Gm + 1/3*(mq.*RHO-M)./max(tau_r,4*dt);

% total rates of change
dXdt   = advn_X + Gx;
dFdt   = advn_F + Gf;
dMdt   = advn_M + Gm;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_F = (a1*F-a2*Fo-a3*Foo)/dt - (b1*dFdt + b2*dFdto + b3*dFdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
X = (a2*Xo+a3*Xoo + (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt)/a1;
F = (a2*Fo+a3*Foo + (b1*dFdt + b2*dFdto + b3*dFdtoo)*dt)/a1;
M = (a2*Mo+a3*Moo + (b1*dMdt + b2*dMdto + b3*dMdtoo)*dt)/a1;

% apply minimum bound
X   = max(0, X );
F   = max(0, F );
M   = max(0, M );

% get dynamically evolving mixture density 
RHO = X+F+M;

%***  update phase fractions and component concentrations

% update phase fractions
x = X./RHO;
f = F./RHO;
m = M./RHO;

% convert major component density to concentration
% c = C./RHO; c = c./sum(c,3);

% update phase entropies
sm = (S - X.*Dsx - F.*Dsf)./RHO;
sx = sm + Dsx;
sf = sm + Dsf;

% update major component phase composition
Kc  = (cxq+TINY)./(cmq+TINY);
res = 1; tol  = 1e-10;
it  = 1; mxit = 25;
while res>tol && it<maxit
    Kci = Kc;
    Kc  = Kc .* sum((c-f.*cfq)./(m + x.*Kc),3)./sum((c-f.*cfq)./(m./Kc + x),3);
    res = norm(Kci-Kc,'fro')./norm(Kc,'fro');
    it = it+1;
end
cm  = (c-f.*cfq)./(m + x.*Kc);
cx  = (c-f.*cfq)./(m./Kc + x);

TCtime = TCtime + toc - eqtime;
