%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; xi = x; fi = f;

% update heat content (entropy)
advn_S = - advection(rho.*m.*sm,Um,Wm,h,ADVN,'flx') ...                    % heat advection
         - advection(rho.*x.*sx,Ux,Wx,h,ADVN,'flx') ...
         - advection(rho.*f.*sf,Uf,Wf,h,ADVN,'flx');

qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                     % heat diffusion z-flux
qTx    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddx(T,h);                     % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h)  ...                    % heat diffusion
                           - ddx(qTx(2:end-1,:),h));

diss_h(2:end-1,2:end-1) = diss ./ T(2:end-1,2:end-1);

bndT = zeros(size(T));
if ~isnan(Twall); bndT = bndT + ((Twall+273.15)-T)./tau_T .* bndshape; end % impose top boundary layer
bndS = rho.*cP.*bndT./T;

dSdt = advn_S + diff_T + diss_h + bndS;                                    % total rate of change

S = So + (THETA.*dSdt + (1-THETA).*dSdto).*dt;                             % explicit update of major component density
S([1 end],:) = S([2 end-1],:);                                             % apply zero flux boundary conditions
S(:,[1 end]) = S(:,[2 end-1]);

    
% update major component
advn_C = advection(rho.*m.*cm,Um,Wm,h,ADVN,'flx') ...
       + advection(rho.*x.*cx,Ux,Wx,h,ADVN,'flx');
                          
qcz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(c,h);                     % major component diffusion z-flux
qcx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(c,h);                     % major component diffusion x-flux
diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                      % major component diffusion
                          - ddx(qcx(2:end-1,:),h);
    
bndC = zeros(size(c));
if ~isnan(cwall); bndC = bndC + rho.*(cwall-c)./tau_a .* bndshape; end     % impose boundary layer

dCdt = - advn_C + diff_c + bndC;                                           % total rate of change
    
C = Co + (THETA.*dCdt + (1-THETA).*dCdto).*dt;                             % explicit update of major component density
C([1 end],:) = C([2 end-1],:);                                             % apply boundary conditions
C(:,[1 end]) = C(:,[2 end-1]);  
    
% update volatile component
bndV = zeros(size(v));
if any([v0;v1;vwall;v(:)]>10*TINY)
    advn_V = advection(rho.*m.*vm,Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*f.*vf,Uf,Wf,h,ADVN,'flx');
    
    qvz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(v,h);                 % volatile component diffusion z-flux
    qvx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(v,h);                 % volatile component diffusion x-flux
    diff_v(2:end-1,2:end-1) = - ddz(qvz(:,2:end-1),h) ...                  % volatile component diffusion
                              - ddx(qvx(2:end-1,:),h);
    
    if ~isnan(vwall); bndV = bndV + rho.*(vwall-v)./tau_a .* bndshape; end % impose boundary layer
    
    dVdt = - advn_V + diff_v + bndV;                                       % total rate of change
    
    V = Vo + (THETA.*dVdt + (1-THETA).*dVdto).*dt;                         % explicit update of volatile component density
    V = max(TINY,V);
    V([1 end],:) = V([2 end-1],:);                                         % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);
end

% convert enthalpy and component densities to temperature and concentrations
T = (T0+273.15)*exp(S./rho./cP - x.*Dsx./cP - f.*Dsf./cP + aT./rhoref./cP.*(Pt-Ptop));  % convert entropy to temperature
c = C./rho;
v = V./rho;


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
[xq,cxq,cmq,fq,vfq,vmq] = equilibrium(xq,fq,T-273.15,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,TINY);

% update crystal fraction
if diseq
    
    Gx = ALPHA.*Gx + (1-ALPHA) .* (xq-x).*rho./max(3.*dt,tau_r);
    
    advn_x = - advection(rho.*x,Ux,Wx,h,ADVN,'flx');                       % get advection term
    
    dXdt   = advn_x + Gx;                                                  % total rate of change
    
    X = Xo + (THETA.*dXdt + (1-THETA).*dXdto).*dt;                         % explicit update of crystal fraction
    X = min(rho-TINY,max(TINY,X));                                         % enforce [0,1] limit
    X([1 end],:) = X([2 end-1],:);                                         % apply boundary conditions
    X(:,[1 end]) = X(:,[2 end-1]);
    
    x  = X./rho;
    
else
    
    x  =  ALPHA.*x + (1-ALPHA).*xq;
    Gx = (rho.*x-rhoo.*xo)./dt + advection(rho.*x,Ux,Wx,h,ADVN,'flx');     % reconstruct crystallisation rate
    
end

% update bubble fraction
if (diseq && any([v0;v1;vwall;v(:)]>10*TINY))
    
    Gf = ALPHA.*Gf + (1-ALPHA) .* (fq-f).*rho./max(3.*dt,tau_r);
    
    advn_f = - advection(rho.*f,Uf,Wf,h,ADVN,'flx');                       % get advection term
                          
    dFdt   = advn_f + Gf;                                                  % total rate of change
    
    F = Fo + (THETA.*dFdt + (1-THETA).*dFdto).*dt;                         % explicit update of bubble fraction
    F = min(rho-TINY,max(TINY,F));                                         % enforce [0,1-x] limit
    F([1 end],:) = F([2 end-1],:);                                         % apply boundary conditions
    F(:,[1 end]) = F(:,[2 end-1]);
    
    f  = F./rho;

else
    
    f  =  ALPHA.*f + (1-ALPHA).*fq;
    Gf = (rho.*f-rhoo.*fo)./dt + advection(rho.*f,Uf,Wf,h,ADVN,'flx');     % reconstruct exsolution rate
    
end

% update melt fraction
m = min(1-x-TINY,max(TINY,1-f-x));


if step>0

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
    
end

% get residual of thermochemical equations from iterative update
resnorm_TC = norm(T - Ti,2)./(norm(T,2)+TINY) ...
           + norm((x - xi).*(x>10*TINY),2)./(norm(x,2)+TINY) ...
           + norm((f - fi).*(f>10*TINY),2)./(norm(f,2)+TINY);
