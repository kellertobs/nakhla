%% *****  THERMO-CHEMICAL EVOLUTION  **************************************
tic;

% store previous iteration
Si = S; Ci = C; Vi = V; Xi = X; Fi = F;


%***  update heat content (entropy)

% heat advection
advn_S = - advect(M(inz,inx).*sm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X(inz,inx).*sx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA) ...  % solid advection
         - advect(F(inz,inx).*sf(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);     % fluid advection

qSz    = - (ks(1:end-1,:)+ks(2:end,:))./2 .* ddz(T,h);  % z-flux
qSx    = - (ks(:,1:end-1)+ks(:,2:end))./2 .* ddx(T,h);  % x-flux
diff_S =(- ddz(qSz(:,inx),h)  ...
         - ddx(qSx(inz,:),h));

% heat dissipation
diss_h = diss ./ T(inz,inx);

% boundary layers
if ~isnan(Twall)
    bnd_T = ((Twall+273.15)-T(inz,inx))./tau_T .* bndshape;
    bnd_S = rho(inz,inx).*cP.*bnd_T./T(inz,inx);
end

% total rate of change
dSdt = advn_S + diff_S + diss_h + bnd_S;

% semi-implicit update of bulk entropy density
S(inz,inx)   = So(inz,inx) + (theta.*dSdt + (1-theta).*dSdto).*dt;    

% boundary conditions
S([1 end],:) = S([2 end-1],:);                                             
S(:,[1 end]) = S(:,[2 end-1]);


%***  update major component (SiO2)

% major component advection
advn_C = - advect(M(inz,inx).*cm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
         - advect(X(inz,inx).*cx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);     % solid advection

% major component diffusion
% qCz    = - (kc(1:end-1,:)+kc(2:end,:))./2 .* ddz(c,h);  % z-flux
% qCx    = - (kc(:,1:end-1)+kc(:,2:end))./2 .* ddx(c,h);  % x-flux
% diff_C = (- ddz(qCz(:,inx),h)  ...
%           - ddx(qCx(inz,:),h));

% boundary layers
if ~isnan(cwall); bnd_C = (rho(inz,inx).*cwall-C(inz,inx))./tau_a .* bndshape; end

% total rate of change
dCdt = advn_C + bnd_C;                                            
    
% semi-implicit update of major component density
C(inz,inx)   = Co(inz,inx) + (theta.*dCdt + (1-theta).*dCdto).*dt;         
C            = max(cal.cphs0.*rho,min(cal.cphs1.*rho,C));

% boundary conditions
C([1 end],:) = C([2 end-1],:);
C(:,[1 end]) = C(:,[2 end-1]);  
    

%***  update volatile component (H2O)
if any([v0;v1;vwall;v(:)]>10*TINY)

    % volatile component advection
    advn_V = - advect(M(inz,inx).*vm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
             - advect(F(inz,inx).*vf(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);     % fluid advection
    
    % volatile component diffusion
%     qVz    = - (kv(1:end-1,:)+kv(2:end,:))./2 .* ddz(v,h);  % z-flux
%     qVx    = - (kv(:,1:end-1)+kv(:,2:end))./2 .* ddx(v,h);  % x-flux
%     diff_V = (- ddz(qVz(:,inx),h)  ...
%               - ddx(qVx(inz,:),h));

    % boundary layers
    if ~isnan(vwall); bnd_V = (rho(inz,inx).*vwall-V(inz,inx))./tau_a .* bndshape; end 
    
    % total rate of change
    dVdt = advn_V + bnd_V;                                                 
    
    % semi-implicit update of volatile component density
    V(inz,inx)   = Vo(inz,inx) + (theta.*dVdt + (1-theta).*dVdto).*dt;
    V            = max(0,min(rho,V));

    % boundary conditions
    V([1 end],:) = V([2 end-1],:);                                         
    V(:,[1 end]) = V(:,[2 end-1]);
end


% convert entropy and component densities to temperature and concentrations
T = (cal.Tphs1+273.15)*exp((S - X.*Dsx - F.*Dsf)./rho./cP + Adbt./cP.*(Pt-Ptop));
c = C./rho;
v = V./rho;


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

eqtime = tic;

%*** update phase equilibrium
[xq(inz,inx),cxq(inz,inx),cmq(inz,inx),fq(inz,inx),vfq(inz,inx),vmq(inz,inx)] = equilibrium(xq(inz,inx),fq(inz,inx),T(inz,inx)-273.15,c(inz,inx),v(inz,inx),Pt(inz,inx),cal,TINY);

% boundary conditions
xq([1 end],:) = xq([2 end-1],:);
xq(:,[1 end]) = xq(:,[2 end-1]);
fq([1 end],:) = fq([2 end-1],:);
fq(:,[1 end]) = fq(:,[2 end-1]);

cxq([1 end],:) = cxq([2 end-1],:);
cxq(:,[1 end]) = cxq(:,[2 end-1]);
cmq([1 end],:) = cmq([2 end-1],:);
cmq(:,[1 end]) = cmq(:,[2 end-1]);

vfq([1 end],:) = vfq([2 end-1],:);
vfq(:,[1 end]) = vfq(:,[2 end-1]);
vmq([1 end],:) = vmq([2 end-1],:);
vmq(:,[1 end]) = vmq(:,[2 end-1]);

mq = 1-xq-fq;

EQtime = EQtime + toc(eqtime);


%***  update crystal fraction

% crystallisation rate
Gx = lambda * Gx + (1-lambda) * (rho(inz,inx).*xq(inz,inx)-X(inz,inx))./max(tau_r,5*dt);

% crystallinity advection
advn_X = - advect(X(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

% crystallinity diffusion
% qXz    = - (kx(1:end-1,:)+kx(2:end,:))./2 .* ddz(x,h);  % z-flux
% qXx    = - (kx(:,1:end-1)+kx(:,2:end))./2 .* ddx(x,h);  % x-flux
% diff_X = (- ddz(qXz(:,inx),h)  ...
%           - ddx(qXx(inz,:),h));

% total rate of change
dXdt   = advn_X + Gx;

% semi-implicit update of crystal fraction
X(inz,inx) = Xo(inz,inx) + (theta.*dXdt + (1-theta).*dXdto).*dt;
X = max(0, X );

% boundary conditions
X([1 end],:) = X([2 end-1],:);
X(:,[1 end]) = X(:,[2 end-1]);


%***  update bubble fraction
if any([v0;v1;vwall;v(:)]>10*TINY)

    % fluid exsolution rate
    Gf = lambda * Gf + (1-lambda) * (rho(inz,inx).*fq(inz,inx)-F(inz,inx))./max(tau_r,5*dt);

    % fluid bubble advection
    advn_F = - advect(F(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);

    % fluid bubble diffusion
%     qFz    = - (kf(1:end-1,:)+kf(2:end,:))./2 .* ddz(f,h);  % z-flux
%     qFx    = - (kf(:,1:end-1)+kf(:,2:end))./2 .* ddx(f,h);  % x-flux
%     diff_F = (- ddz(qFz(:,inx),h)  ...
%               - ddx(qFx(inz,:),h));

    % total rate of change
    dFdt   = advn_F + Gf;

    % semi-implicit update of bubble fraction
    F(inz,inx) = Fo(inz,inx) + (theta.*dFdt + (1-theta).*dFdto).*dt;
    F = max(0,min(V, F ));

    % boundary conditions
    F([1 end],:) = F([2 end-1],:);                                         
    F(:,[1 end]) = F(:,[2 end-1]);
    
end

M = rho-X-F;

% update phase fractions
x = max(0,min(1,X./rho));
f = max(0,min(1,F./rho));
m = max(0,min(1,M./rho));

% update phase entropies
sm = (S - X.*Dsx - F.*Dsf)./rho;
sx = sm + Dsx;
sf = sm + Dsf;

% update phase compositions
% major component
Kc = cxq./cmq;
cm = c./(m + x.*Kc); cm(m==0) = cmq(m==0);
cx = c./(m./Kc + x); cx(m==1) = cxq(m==1);

% volatile component
Kf = vfq./max(TINY,vmq);
vm = v./max(TINY,m + f.*Kf); vm(m==0) = vmq(m==0);
vf = v./max(TINY,m./Kf + f); vf(m==1) = vfq(m==1);


%% *****  UPDATE TC RESIDUALS  ********************************************

% get residual of thermochemical equations from iterative update
resnorm_TC = norm( S(inz,inx) - Si(inz,inx),2)./(norm(S(inz,inx),2)+TINY) ...
           + norm( C(inz,inx) - Ci(inz,inx),2)./(norm(C(inz,inx),2)+TINY) ...
           + norm( V(inz,inx) - Vi(inz,inx),2)./(norm(V(inz,inx),2)+TINY) ...
           + norm((X(inz,inx) - Xi(inz,inx)).*(x(inz,inx)>10*TINY).*(m(inz,inx)>10*TINY),2)./(norm(X(inz,inx),2)+TINY) ...
           + norm((F(inz,inx) - Fi(inz,inx)).*(f(inz,inx)>10*TINY).*(m(inz,inx)>10*TINY),2)./(norm(F(inz,inx),2)+TINY);

TCtime = TCtime + toc - toc(eqtime);

