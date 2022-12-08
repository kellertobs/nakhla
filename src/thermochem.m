%% *****  THERMO-CHEMICAL EVOLUTION  **************************************
tic;

% store previous iteration
Ti = T; ci = c; vi = v; xi = x; fi = f;


% update heat content (entropy)
advn_S = - advect(M(inz,inx).*sm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % heat advection
         - advect(X(inz,inx).*sx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA) ...
         - advect(F(inz,inx).*sf(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);

qSz    = - (ks(1:end-1,:)+ks(2:end,:))./2 .* ddz(T,h);
qSx    = - (ks(:,1:end-1)+ks(:,2:end))./2 .* ddx(T,h);
diff_S = (- ddz(qSz(:,inx),h)  ...
          - ddx(qSx(inz,:),h));

diss_h = diss ./ T(inz,inx);

if ~isnan(Twall)
    bnd_T = ((Twall+273.15)-T(inz,inx))./tau_T .* bndshape;  % impose top boundary layer
    bnd_S = rho(inz,inx).*cP.*bnd_T./T(inz,inx);
end

dSdt = advn_S + diff_S + diss_h + bnd_S;                                   % total rate of change

S(inz,inx)   = So(inz,inx) + (theta.*dSdt + (1-theta).*dSdto).*dt;         % explicit update of major component density
S([1 end],:) = S([2 end-1],:);                                             % apply zero flux boundary conditions
S(:,[1 end]) = S(:,[2 end-1]);

% update major component
advn_C = - advect(M(inz,inx).*cm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % major component advection
         - advect(X(inz,inx).*cx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

qCz    = - (kc(1:end-1,:)+kc(2:end,:))./2 .* ddz(c,h);
qCx    = - (kc(:,1:end-1)+kc(:,2:end))./2 .* ddx(c,h);
diff_C = (- ddz(qCz(:,inx),h)  ...
          - ddx(qCx(inz,:),h));

if ~isnan(cwall); bnd_C = rho(inz,inx).*(cwall-c(inz,inx))./tau_a .* bndshape; end % impose boundary layer

dCdt = advn_C + 0.*diff_C + bnd_C;                                            % total rate of change
    
C(inz,inx)   = Co(inz,inx) + (theta.*dCdt + (1-theta).*dCdto).*dt;         % explicit update of major component density
C            = max(cal.cphs0.*rho,min(cal.cphs1.*rho,C));
C([1 end],:) = C([2 end-1],:);                                             % apply boundary conditions
C(:,[1 end]) = C(:,[2 end-1]);  
    
% update volatile component
if any([v0;v1;vwall;v(:)]>10*TINY)
    advn_V = - advect(M(inz,inx).*vm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % volatile component advection
             - advect(F(inz,inx).*vf(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);
    
    qVz    = - (kc(1:end-1,:)+kc(2:end,:))./2 .* ddz(v,h);
    qVx    = - (kc(:,1:end-1)+kc(:,2:end))./2 .* ddx(v,h);
    diff_V = (- ddz(qVz(:,inx),h)  ...
              - ddx(qVx(inz,:),h));

    if ~isnan(vwall); bnd_V = rho(inz,inx).*(vwall-v(inz,inx))./tau_a .* bndshape; end % impose boundary layer
    
    dVdt = advn_V + diff_V + bnd_V;                                                 % total rate of change
    
    V(inz,inx)   = Vo(inz,inx) + (theta.*dVdt + (1-theta).*dVdto).*dt;     % explicit update of volatile component density
    V            = max(0,min(rho,V));
    V([1 end],:) = V([2 end-1],:);                                         % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);
end

% convert entropy and component densities to temperature and concentrations
T = (T0+273.15)*exp((S - X.*Dsx - F.*Dsf)./rho./cP + Adbt./cP.*(Pt-Ptop));  % convert entropy to temperature
c = C./rho;
v = V./rho;


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
eqtime = tic;

% extract indices for which equilibrium needs updating
[xq(inz,inx),cxq(inz,inx),cmq(inz,inx),fq(inz,inx),vfq(inz,inx),vmq(inz,inx)] = equilibrium(xq(inz,inx),fq(inz,inx),T(inz,inx)-273.15,c(inz,inx),v(inz,inx),Pt(inz,inx),cal,TINY);

% apply boundary conditions
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

EQtime = EQtime + toc(eqtime);

% update crystal fraction
Gx = lambda * Gx + (1-lambda) * (xq(inz,inx).*rho(inz,inx)-X(inz,inx))./(5*dt);

advn_X = - advect(X(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

qXz    = - (kx(1:end-1,:)+kx(2:end,:))./2 .* ddz(x,h);
qXx    = - (kx(:,1:end-1)+kx(:,2:end))./2 .* ddx(x,h);
diff_X = (- ddz(qXz(:,inx),h)  ...
          - ddx(qXx(inz,:),h));

dXdt   = advn_X + 0.*diff_X + Gx;                                             % total rate of change

X(inz,inx) = Xo(inz,inx) + (theta.*dXdt + (1-theta).*dXdto).*dt;           % explicit update of crystal fraction
X = max(0,min(rho-F, X ));                                                 % enforce limits
X([1 end],:) = X([2 end-1],:);                                             % apply boundary conditions
X(:,[1 end]) = X(:,[2 end-1]);

% update bubble fraction
if any([v0;v1;vwall;v(:)]>10*TINY)

    Gf = lambda * Gf + (1-lambda) * (fq(inz,inx).*rho(inz,inx)-F(inz,inx))./(5*dt);

    advn_F = - advect(F(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);

    qFz    = - (kf(1:end-1,:)+kf(2:end,:))./2 .* ddz(f,h);
    qFx    = - (kf(:,1:end-1)+kf(:,2:end))./2 .* ddx(f,h);
    diff_F = (- ddz(qFz(:,inx),h)  ...
              - ddx(qFx(inz,:),h));

    dFdt   = advn_F + diff_F + Gf;                                         % total rate of change

    F(inz,inx) = Fo(inz,inx) + (theta.*dFdt + (1-theta).*dFdto).*dt;       % explicit update of bubble fraction
    F = max(0,min(V, F ));                                                 % enforce limits
    F([1 end],:) = F([2 end-1],:);                                         % apply boundary conditions
    F(:,[1 end]) = F(:,[2 end-1]);
    
else 
    diff_F = 0;
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
resnorm_TC = norm( T(inz,inx) - Ti(inz,inx)                       ,2)./(norm(T(inz,inx),2)+TINY) ...
           + norm((c(inz,inx) - ci(inz,inx))                      ,2)./(norm(c(inz,inx),2)+TINY) ...
           + norm((v(inz,inx) - vi(inz,inx))                      ,2)./(norm(v(inz,inx),2)+TINY) ...
           + norm((x(inz,inx) - xi(inz,inx)).*(x(inz,inx)>10*TINY).*(m(inz,inx)>10*TINY),2)./(norm(x(inz,inx),2)+TINY) ...
           + norm((f(inz,inx) - fi(inz,inx)).*(f(inz,inx)>10*TINY).*(m(inz,inx)>10*TINY),2)./(norm(f(inz,inx),2)+TINY);

TCtime = TCtime + toc - toc(eqtime);

