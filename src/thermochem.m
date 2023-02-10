%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

for itTC = 1:inner_TC

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
S = S - lambda*res_S*dt;


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
C = C - lambda*res_C*dt;
C          = max(cal.cphs0.*rho,min(cal.cphs1.*rho,C));
    

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
    V = V - lambda*res_V*dt;
    V          = max(0,min(rho,V));
end


% convert entropy and component densities to temperature and concentrations
T = (cal.Tphs1+273.15)*exp((S - X.*Dsx - F.*Dsf)./rho./cP + Adbt./cP.*(Pt-Ptop));
c = C./rho;
v = V./rho;

eqtime = tic;

%*** update phase equilibrium
[xq,cxq,cmq,fq,vfq,vmq] = equilibrium(xq,fq,T-273.15,c,v,Pt,cal,TINY);

mq = 1-xq-fq;

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;


%***  update phase fractions

% phase mass transfer rates
Gx = (xq.*rho-X)./max(tau_r,4*dt);
Gf = (fq.*rho-F)./max(tau_r,4*dt);
Gm = -Gx-Gf;

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

% update of phase density evolution
X = X - lambda*res_X*dt;
F = F - lambda*res_F*dt;
M = M - lambda*res_M*dt;

% apply minimum bound
X = max(0, X );
F = max(0, F );
M = max(0, M );

% update phase fractions
x = X./rho;
f = F./rho;
m = M./rho;

% normalise to unit sum
sumphs = x+f+m;
x = x./sumphs;
f = f./sumphs;
m = m./sumphs;

% update phase entropies
sm = (S - X.*Dsx - F.*Dsf)./rho;
sx = sm + Dsx;
sf = sm + Dsf;

% update phase compositions

% major component
Kc = cxq./cmq;
cm = c./(m + x.*Kc);
cx = c./(m./Kc + x);

% volatile component
Kf = vfq./max(TINY,vmq);
vm = v./max(TINY,m + f.*Kf);
vf = v./max(TINY,m./Kf + f);


TCtime = TCtime + toc;

end

