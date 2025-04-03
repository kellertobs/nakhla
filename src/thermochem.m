%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

if iter==1; upd_S = 0; upd_C = 0; upd_M = 0; upd_X = 0; upd_F = 0; end

%***  update heat content (entropy) density

% heat advection
advn_S = - advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % solid advection
         - advect(F.*sf,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

diff_S = diffus(T,kT./T,h,[1,2],BCD) + diffus(Tp,ks,h,[1,2],BCD);

% heat dissipation
diss_h = diss ./ T;

% boundary layers
bnd_T = zeros(size(S));
if ~isnan(Twall(1)); bnd_T = bnd_T + ((Twall(1)+273.15)-T)./tau_T .* topshape; end
if ~isnan(Twall(2)); bnd_T = bnd_T + ((Twall(2)+273.15)-T)./tau_T .* botshape; end
if ~isnan(Twall(3)); bnd_T = bnd_T + ((Twall(3)+273.15)-T)./tau_T .* sdsshape; end
bnd_S = RHO.*cP.*bnd_T ./ T;

% total rate of change
dSdt  = advn_S + diff_S + diss_h + bnd_S;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo)/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% semi-implicit update of bulk entropy density
upd_S = - alpha*res_S*dt/a1 + beta*upd_S;
S     = S + upd_S;

% convert entropy S to natural temperature T and potential temperature Tp
[Tp,~ ] = StoT(Tp,S./rho,Pref+0*Pt,cat(3,m,x,f),[cPm;cPx;cPf],[aTm;aTx;aTf],[bPm;bPx;bPf],cat(3,rhom0,rhox0,rhof0),[sref;sref+Dsx;sref+Dsf],Tref,Pref);
[T ,si] = StoT(T ,S./rho,       Pt,cat(3,m,x,f),[cPm;cPx;cPf],[aTm;aTx;aTf],[bPm;bPx;bPf],cat(3,rhom0,rhox0,rhof0),[sref;sref+Dsx;sref+Dsf],Tref,Pref);
sm = si(:,:,1); sx = si(:,:,2); sf = si(:,:,3);  % read out phase entropies


%***  update major component densities

% major component advection
advn_C = - advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % solid advection
         - advect(F.*cf,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

% major component diffusion (regularisation)
diff_C = diffus(cm,M.*kc,h,[1,2],BCD) + diffus(cx,X.*kc,h,[1,2],BCD);

% boundary layers
bnd_C = zeros(size(C));
for i = 1:cal.ncmp
    if ~isnan(cwall(1)); bnd_C(:,:,i) = bnd_C(:,:,i) + (RHO.*cwall(1,i)-C(:,:,i)).*mu./tau_a .* topshape; end
    if ~isnan(cwall(2)); bnd_C(:,:,i) = bnd_C(:,:,i) + (RHO.*cwall(2,i)-C(:,:,i)).*mu./tau_a .* botshape; end
    if ~isnan(cwall(3)); bnd_C(:,:,i) = bnd_C(:,:,i) + (RHO.*cwall(3,i)-C(:,:,i)).*mu./tau_a .* sdsshape; end
end

% total rate of change
dCdt = advn_C + diff_C + bnd_C;                                            
  
% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
upd_C = max(-C, - alpha*res_C*dt/a1 + beta*upd_C );
C     = C + upd_C;

% convert component density to concentration
c = C./sum(C,3);

%*** update phase equilibrium
if Rcouple; phseql; end

%***  update phase fraction densities

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_F   = - advect(F,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_F+advn_M;

% total rates of change
dXdt   = advn_X + Gx;
dFdt   = advn_F + Gf;
dMdt   = advn_M + Gm;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_F = (a1*F-a2*Fo-a3*Foo)/dt - (b1*dFdt + b2*dFdto + b3*dFdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
upd_X = max(-X/2, - alpha*res_X*dt/a1 + beta*upd_X );
upd_F = max(-F/2, - alpha*res_F*dt/a1 + beta*upd_F );
upd_M = max(-M/2, - alpha*res_M*dt/a1 + beta*upd_M );

X     = X + upd_X;
F     = F + upd_F;
M     = M + upd_M;

% get dynamically evolving mixture density 
RHO = X+F+M;

%***  update phase fractions and component concentrations

% update phase fractions
x = X./RHO; 
f = F./RHO; 
m = M./RHO;

% identify subsolidus and superliquidus regions
subsol  = m<=1e-9 & T<=reshape(cal.Tsol+273.15,Nz,Nx);
supliq  = x<=1e-9 & T>=reshape(cal.Tliq+273.15,Nz,Nx);
subsolc = repmat(subsol,1,1,cal.ncmp);
supliqc = repmat(supliq,1,1,cal.ncmp);

% update major component phase composition
rnorm   = 1;  tol  = sqrt(eps);
it      = 1;  mxit = 50;
cm = cmq;  cx = cxq;
while rnorm>tol && it<mxit

    Kx = cx./(cm+eps);
    Kf = cf./(cm+eps);

    cm = c    ./(m + x.*Kx + f.*Kf + eps);
    cx = c.*Kx./(m + x.*Kx + f.*Kf + eps);

    cm = cm./sum(cm,3);
    cx = cx./sum(cx,3);

    r = x.*cx + m.*cm + f.*cf - c;
    r(subsolc) = 0; r(supliqc) = 0;

    rnorm = norm(r(:))./norm(c(:));
    it  = it+1;
end

% fix subsolidus and superliquidus conditions
cx(subsolc) = cxq(subsolc); x(subsol) = xq(subsol); f(subsol) = fq(subsol); m(subsol) = 0;
cm(supliqc) = cmq(supliqc); m(supliq) = mq(supliq); f(supliq) = fq(supliq); x(supliq) = 0;

% if (it==mxit && rnorm>tol)
%     disp(['!!! Lever rule adjustment converged to ',num2str(rnorm),' after ',num2str(mxit),' iterations !!!']);
% end

% record timing
TCtime = TCtime + toc;% - eqtime;
