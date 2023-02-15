%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

tic;

%***  update heat content (entropy)

% heat advection
advn_S = - advect(M.*sm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*sx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % solid advection
         - advect(F.*sf,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

diff_S = diffus(T,ks,h,[1,2],BCD);

% heat dissipation
diss_h = diss ./ T;

% boundary layers
if ~isnan(Twall)
    bnd_T = ((Twall+273.15)-T)./tau_T .* bndshape;
    bnd_S = rho.*cP.*bnd_T./T;
end

% total rate of change
dSdt  = advn_S + diff_S + diss_h + bnd_S;

% residual of entropy evolution
res_S = (a1*S-a2*So-a3*Soo)/dt - (b1*dSdt + b2*dSdto + b3*dSdtoo);

% semi-implicit update of bulk entropy density
S = S - alpha*res_S*dt;
% S = (a2*So+a3*Soo + (b1*dSdt + b2*dSdto + b3*dSdtoo)*dt)/a1;

% convert entropy desnity to temperature
T = (cal.Tphs1+273.15)*exp((S - X.*Dsx - F.*Dsf)./rho./cP + Adbt./cP.*(Pt-Ptop));


%***  update major component (SiO2)

% major component advection
advn_C = - advect(M.*cm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X.*cx,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % solid advection

% boundary layers
if ~isnan(cwall); bnd_C = (rho.*cwall-C)./tau_a .* bndshape; end

% total rate of change
dCdt = advn_C + bnd_C;                                            

% residual of major component evolution
res_C = (a1*C-a2*Co-a3*Coo)/dt - (b1*dCdt + b2*dCdto + b3*dCdtoo);

% semi-implicit update of major component density
C = C - alpha*res_C*dt;
% C = (a2*Co+a3*Coo + (b1*dCdt + b2*dCdto + b3*dCdtoo)*dt)/a1;

% apply minimum/maximum bounds
C = max(cal.cphs0.*rho,min(cal.cphs1.*rho,C));

% convert major component density to concentration
c = C./rho;


%***  update volatile component (H2O)
if any([v0;v1;vwall;v(:)]>10*TINY)

    % volatile component advection
    advn_V = - advect(M.*vm,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
             - advect(F.*vf,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);     % fluid advection

    % boundary layers
    if ~isnan(vwall); bnd_V = (rho.*vwall-V)./tau_a .* bndshape; end 
    
    % total rate of change
    dVdt = advn_V + bnd_V;                                                 
    
    % residual of volatile component evolution
    res_V = (a1*V-a2*Vo-a3*Voo)/dt - (b1*dVdt + b2*dVdto + b3*dVdtoo);

    % semi-implicit update of volatile component density
    V = V - alpha*res_V*dt;
%     V = (a2*Vo+a3*Voo + (b1*dVdt + b2*dVdto + b3*dVdtoo)*dt)/a1;

    % apply minimum bound
    V = max(0,V );

    % convert volatile component density to concentration
    v = V./rho;
end


%*** update phase equilibrium
eqtime = tic;

[xq,cxq,cmq,fq,vfq,vmq] = equilibrium(xq,fq,T-273.15,c,v,Pt,cal,TINY);

mq = 1-xq-fq;

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;


%***  update phase fractions

% phase advection rates
advn_X   = - advect(X,Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_F   = - advect(F,Uf(2:end-1,:),Wf(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_M   = - advect(M,Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA);
advn_rho = advn_X+advn_F+advn_M;

% phase mass transfer rates
res_Gx = Gx - (xq-x).*rho./max(tau_r,4*dt);
Gx     = Gx - alpha.*res_Gx;
res_Gf = Gf - (fq-f).*rho./max(tau_r,4*dt);
Gf     = Gf - alpha.*res_Gf;
res_Gm = Gm - (mq-m).*rho./max(tau_r,4*dt);
Gm     = Gm - alpha.*res_Gm;
% Gm = -Gx-Gf;

% total rates of change
dXdt   = advn_X + Gx;
dFdt   = advn_F + Gf;
dMdt   = advn_M + Gm;

% residual of phase density evolution
res_X = (a1*X-a2*Xo-a3*Xoo)/dt - (b1*dXdt + b2*dXdto + b3*dXdtoo);
res_F = (a1*F-a2*Fo-a3*Foo)/dt - (b1*dFdt + b2*dFdto + b3*dFdtoo);
res_M = (a1*M-a2*Mo-a3*Moo)/dt - (b1*dMdt + b2*dMdto + b3*dMdtoo);

% semi-implicit update of phase fraction densities
X = X - alpha*res_X*dt;
F = F - alpha*res_F*dt;
M = M - alpha*res_M*dt;
% X = (a2*Xo+a3*Xoo + (b1*dXdt + b2*dXdto + b3*dXdtoo)*dt)/a1;
% F = (a2*Fo+a3*Foo + (b1*dFdt + b2*dFdto + b3*dFdtoo)*dt)/a1;
% M = (a2*Mo+a3*Moo + (b1*dMdt + b2*dMdto + b3*dMdtoo)*dt)/a1;

% apply minimum bound
X = max(0, X );
F = max(0, F );
M = max(0, M );


% update phase fractions
x = X./(X+M+F);
f = F./(X+M+F);
m = M./(X+M+F);

% update phase entropies
sm = (S - X.*Dsx - F.*Dsf)./rho;
sx = sm + Dsx;
sf = sm + Dsf;

% update major component phase composition
Kc = cxq./cmq;
cm = c./(m + x.*Kc);
cx = c./(m./Kc + x);

% update volatile component phase composition
Kf = vfq./max(TINY,vmq);
vm = v./max(TINY,m + f.*Kf);
vf = v./max(TINY,m./Kf + f);

TCtime = TCtime + toc - eqtime;
