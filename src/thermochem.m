%% *****  THERMO-CHEMICAL EVOLUTION  **************************************
tic;

% store previous iteration
Ti = T; ci = c; vi = v; xi = x; fi = f;

% update heat content (entropy)
advn_S = - advect(rho(inz,inx).*m(inz,inx).*sm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % heat advection
         - advect(rho(inz,inx).*x(inz,inx).*sx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA) ...
         - advect(rho(inz,inx).*f(inz,inx).*sf(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);

qTz    = - (ks(1:end-1,:)+ks(2:end,:))./2 .* ddz(T,h);                     % heat diffusion z-flux
qTx    = - (ks(:,1:end-1)+ks(:,2:end))./2 .* ddx(T,h);                     % heat diffusion x-flux
diff_T = (- ddz(qTz(:,inx),h)  ...                                         % heat diffusion
          - ddx(qTx(inz,:),h));

diss_h = diss ./ T(inz,inx);

if ~isnan(Twall); bnd_T = ((Twall+273.15)-T(inz,inx))./tau_T .* bndshape; end % impose top boundary layer
bnd_S = rho(inz,inx).*cP.*bnd_T./T(inz,inx);

dSdt = advn_S + diff_T + diss_h + bnd_S;                                   % total rate of change

S(inz,inx)   = So(inz,inx) + (theta.*dSdt + (1-theta).*dSdto).*dt;         % explicit update of major component density
S([1 end],:) = S([2 end-1],:);                                             % apply zero flux boundary conditions
S(:,[1 end]) = S(:,[2 end-1]);

% update major component
advn_C = - advect(rho(inz,inx).*m(inz,inx).*cm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % major component advection
         - advect(rho(inz,inx).*x(inz,inx).*cx(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

if ~isnan(cwall); bnd_C = rho(inz,inx).*(cwall-c(inz,inx))./tau_a .* bndshape; end % impose boundary layer

dCdt = advn_C + bnd_C;                                                     % total rate of change
    
C(inz,inx)   = Co(inz,inx) + (theta.*dCdt + (1-theta).*dCdto).*dt;         % explicit update of major component density
C([1 end],:) = C([2 end-1],:);                                             % apply boundary conditions
C(:,[1 end]) = C(:,[2 end-1]);  
    
% update volatile component
if any([v0;v1;vwall;v(:)]>10*TINY)
    advn_V = - advect(rho(inz,inx).*m(inz,inx).*vm(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % volatile component advection
             - advect(rho(inz,inx).*f(inz,inx).*vf(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);
    
    if ~isnan(vwall); bnd_V = rho(inz,inx).*(vwall-v(inz,inx))./tau_a .* bndshape; end % impose boundary layer
    
    dVdt = advn_V + bnd_V;                                                 % total rate of change
    
    V(inz,inx)   = Vo(inz,inx) + (theta.*dVdt + (1-theta).*dVdto).*dt;     % explicit update of volatile component density
    V([1 end],:) = V([2 end-1],:);                                         % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);
    V            = max(TINY,V);
end

% convert enthalpy and component densities to temperature and concentrations
T = (T0+273.15)*exp(S./rho./cP - x.*Dsx./cP - f.*Dsf./cP + Adbt./cP.*(Pt-Ptop));  % convert entropy to temperature
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
if diseq

    Gx = lambda * Gx + (1-lambda) * (xq(inz,inx)-x(inz,inx)).*rho(inz,inx)./(4*dt);

    advn_X = - advect(rho(inz,inx).*x(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    dXdt   = advn_X + Gx;                                                  % total rate of change
    
    X(inz,inx) = Xo(inz,inx) + (theta.*dXdt + (1-theta).*dXdto).*dt;       % explicit update of crystal fraction
    X = min(rho,max(0,X));                                                 % enforce limits
    X([1 end],:) = X([2 end-1],:);                                         % apply boundary conditions
    X(:,[1 end]) = X(:,[2 end-1]);

else
    
    X  =  lambda.*X + (1-lambda).*xq.*rho;
    Gx = (X(inz,inx)-Xo(inz,inx))./dt + advect(rho(inz,inx).*x(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);  % reconstruct crystallisation rate
    
end

% update bubble fraction
if any([v0;v1;vwall;v(:)]>10*TINY)
    if diseq

        Gf = lambda * Gf + (1-lambda) * (fq(inz,inx)-f(inz,inx)).*rho(inz,inx)./(4*dt);

        advn_F = - advect(rho(inz,inx).*f(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);

        dFdt   = advn_F + Gf;                                                  % total rate of change

        F(inz,inx) = Fo(inz,inx) + (theta.*dFdt + (1-theta).*dFdto).*dt;       % explicit update of bubble fraction
        F = min(rho-X,max(0,F));                                               % enforce limits
        F([1 end],:) = F([2 end-1],:);                                         % apply boundary conditions
        F(:,[1 end]) = F(:,[2 end-1]);

    else

        F  =  lambda.*F + (1-lambda).*fq.*rho;
        Gf = (F(inz,inx)-Fo(inz,inx))./dt + advect(rho(inz,inx).*f(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);  % reconstruct exsolution rate

    end
end

% update phase fractions
x = min(1,max(0,X./rho));
f = min(1,max(0,F./rho));
m = max(0,min(1,1-f-x));

% update phase entropies
sm = S./rho - x.*Dsx - f.*Dsf;
sx = sm + Dsx;
sf = sm + Dsf;

% update phase compositions
% major component
Kc = cxq./cmq;
cm = c./(m + x.*Kc);
cx = c./(m./Kc + x);

% volatile component
Kf = vfq./vmq;
vf = vfq;
vm = max(0.9.*vmq,min(1.1.*vmq,(v - f)./m));


%% *****  UPDATE TC RESIDUALS  ********************************************

% get residual of thermochemical equations from iterative update
resnorm_TC = norm( T(inz,inx) - Ti(inz,inx)                       ,2)./(norm(T(inz,inx),2)+TINY) ...
           + norm((c(inz,inx) - ci(inz,inx))                      ,2)./(norm(c(inz,inx),2)+TINY) ...
           + norm((v(inz,inx) - vi(inz,inx))                      ,2)./(norm(v(inz,inx),2)+TINY) ...
           + norm((x(inz,inx) - xi(inz,inx)).*(x(inz,inx)>10*TINY),2)./(norm(x(inz,inx),2)+TINY) ...
           + norm((f(inz,inx) - fi(inz,inx)).*(f(inz,inx)>10*TINY),2)./(norm(f(inz,inx),2)+TINY);

TCtime = TCtime + toc - toc(eqtime);

