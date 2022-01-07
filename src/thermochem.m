%% update time step
dtk = min((h/2)^2./max([kT(:)./rhoCp(:);kc]));                             % diffusive time step size
dta = min(min(h/2/max(abs([UBG(:);WBG(:);Um(:);Wm(:);Uf(:);Wf(:);Ux(:);Wx(:)]+1e-16)))); % advective time step size
dt  = min(2*dto,CFL*min([dtk,dta]));                                       % physical time step size


%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% boundary conditions shape function
switch bndmode
    case 0  % none
        bndshape = zeros(size(T));
    case 1  % top only
        bndshape = exp( ( -ZZ)/dw);
    case 2  % bot only
        bndshape = exp(-(D-ZZ)/dw);
    case 3  % top/bot only
        bndshape = exp( ( -ZZ)/dw) ...
                 + exp(-(D-ZZ)/dw);
    case 4  % all walls
        bndshape = exp( ( -ZZ)/dw) ...
                 + exp(-(D-ZZ)/dw) ...
                 + exp( ( -XX)/dw) ...
                 + exp(-(L-XX)/dw);
end
bndshape = min(1,bndshape);


% update temperature
if ~isotherm

    advn_H = advection(rho.*m.*(Cpm + 0  ).*T,Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*x.*(Cpx + Dsx).*T,Ux,Wx,h,ADVN,'flx') ...
           + advection(rho.*f.*(Cpf + Dsf).*T,Uf,Wf,h,ADVN,'flx');
                           
    qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                 % heat diffusion z-flux
    qTx    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddx(T,h);                 % heat diffusion x-flux
    diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...                 % heat diffusion
                               - ddx(qTx(2:end-1,:),h));
    
    cool = zeros(size(T));
    if ~isnan(Twall); cool = cool + rhoCp .* (Twall-T)./tau_T .* bndshape; end  % impose wall cooling

    dHdt = - advn_H + diff_T + cool;                                       % total rate of change
    
    if step>0; H = Ho + (dHdt + dHdto)/2.*dt; end                                        % explicit update of enthalpy
    H([1 end],:) = H([2 end-1],:);                                         % apply boundary conditions
    H(:,[1 end]) = H(:,[2 end-1]);    
    
end


% update composition
if ~isochem
  
    % update major component
    advn_C = advection(rho.*m.*cm,Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*x.*cx,Ux,Wx,h,ADVN,'flx');
                          
    qcz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(c,h);        % major component diffusion z-flux
    qcx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(c,h);        % major component diffusion x-flux
    diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                  % major component diffusion
                              - ddx(qcx(2:end-1,:),h);
    
    assim = zeros(size(c));
    if ~isnan(cwall); assim = assim + (cwall-c).*rho./tau_a .* bndshape; end % impose wall assimilation

    dCdt = - advn_C + diff_c + assim;                                      % total rate of change
    
    if step>0; C = Co + (dCdt + dCdto)/2.*dt; end                                        % explicit update of major component density
    C = min(rho.*cphs1-TINY,max(rho.*cphs0+TINY, C ));
    C([1 end],:) = C([2 end-1],:);                                         % apply boundary conditions
    C(:,[1 end]) = C(:,[2 end-1]);  
    
    % update volatile component
    advn_V = advection(rho.*m.*vm,Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*f.*vf,Uf,Wf,h,ADVN,'flx');
        
    qvz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(v,h);        % volatile component diffusion z-flux
    qvx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(v,h);        % volatile component diffusion x-flux
    diff_v(2:end-1,2:end-1) = - ddz(qvz(:,2:end-1),h) ...                  % volatile component diffusion
                              - ddx(qvx(2:end-1,:),h);

    assim = zeros(size(v));
    if ~isnan(vwall); assim = assim + (vwall-v).*rho./tau_a .* bndshape; end % impose wall assimilation

    dVdt = - advn_V + diff_v + assim;                             % total rate of change
    
    if step>0; V = Vo + (dVdt + dVdto)/2.*dt; end                                        % explicit update of volatile component density
    V = min(rho.*1-TINY,max(rho.*0+TINY, V ));
    V([1 end],:) = V([2 end-1],:);                                         % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);  
    
end


% convert enthalpy and component densities to temperature and concentrations
% T = alpha.*T + (1-alpha).*H./(rhoCp + rhoDs);
% c = alpha.*c + (1-alpha).*C./rho;
% v = alpha.*v + (1-alpha).*V./rho;
T = H./(rhoCp + rhoDs);
c = C./rho;
v = V./rho;


%% *****  UPDATE LOCAL PHASE EQUILIBRIUM  *********************************

[xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,(T+To)./2,(c+co)./2,(v+vo)./2,(Pt+Pto)./2,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
% [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);

    
% update crystal fraction
if diseq
    
    Gx = alpha.*Gx + (1-alpha).*((xq-x).*rho./max(2.*dt,tau_r));
    
    advn_x = advection(rho.*x,Ux,Wx,h,ADVN,'flx');                         % get advection term
    
    dxdt   = - advn_x + Gx;                                                % total rate of change
    
    if step>0; x = (rhoo.*xo + (dxdt + dxdto)/2.*dt)./rho; end                           % explicit update of crystal fraction
    x = min(1-f-TINY,max(TINY,x));                                         % enforce [0,1] limit
    x([1 end],:) = x([2 end-1],:);                                         % apply boundary conditions
    x(:,[1 end]) = x(:,[2 end-1]);
    m = 1-f-x;
    
    Kc = cxq./cmq;
    cm = c./(m + x.*Kc); cm(x < 1e-3) = cmq(x < 1e-3);
    cx = c./(m./Kc + x); cx(m < 1e-3) = cxq(m < 1e-3);
    
else
    
    x  =  alpha.*x + (1-alpha).*xq;
    cx = cxq;
    cm = cmq;
    Gx = (rho.*x-rhoo.*xo)./dt + advection(rho.*x,Ux,Wx,h,ADVN,'flx');     % reconstruct crystallisation rate
    
end


% update bubble fraction
if diseq
    
    Gf = alpha.*Gf + (1-alpha).*((fq-f).*rho./max(2.*dt,tau_r));
    
    advn_f = advection(rho.*f,Uf,Wf,h,ADVN,'flx');                         % get advection term
    
    dfdt   = - advn_f + Gf;                                                % total rate of change
    
    if step>0; f = (rhoo.*fo + (dfdt + dfdto)/2.*dt)./rho; end                           % explicit update of bubble fraction
    f = min(1-x-TINY,max(TINY,f));                                         % enforce [0,1-x] limit
    f([1 end],:) = f([2 end-1],:);                                         % apply boundary conditions
    f(:,[1 end]) = f(:,[2 end-1]);
    m = 1-f-x;

    Kv = vfq./vmq;
    vf = vfq;
    vm = (v - f.*vf)./max(TINY,m); vm(m < 1e-3) = vmq(m < 1e-3);
    
else
    
    f  =  alpha.*f + (1-alpha).*fq;
    vf = vfq;
    vm = vmq;
    Gf = (rho.*f-rhoo.*fo)./dt + advection(rho.*f,Uf,Wf,h,ADVN,'flx');     % reconstruct exsolution rate
    
end


% update melt fraction
m = 1-x-f;



%% ***** TRACE & ISOTOPE GEOCHEMISTRY  ************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
itm = it./(m + x.*KIT);
itx = it./(m./KIT + x);

% update incompatible trace element composition
advn_IT = advection(rho.*m.*itm,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*itx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(it,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(it,h);
diff_it(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(itwall); assim = assim + (itwall-it).*rho./tau_a .* bndshape; end % impose wall assimilation
    
dITdt = - advn_IT + diff_it + assim;                                       % total rate of change

if step>0; IT = ITo + (dITdt + dITdto)/2.*dt; end                                          % explicit update
IT = max(0+TINY, IT );
IT([1 end],:) = IT([2 end-1],:);                                           % boundary conditions
IT(:,[1 end]) = IT(:,[2 end-1]);

it = IT./rho;


% *****  COMPATIBLE TRACE ELEMENT  ****************************************

% update compatible trace element phase compositions
ctm = ct./(m + x.*KCT);
ctx = ct./(m./KCT + x);

% update compatible trace element composition
advn_CT = advection(rho.*m.*ctm,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*ctx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(ct,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(ct,h);
diff_ct(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(ct));
if ~isnan(ctwall); assim = assim + (ctwall-ct).*rho./tau_a .* bndshape; end % impose wall assimilation

dCTdt = - advn_CT + diff_ct + assim;                                       % total rate of change

if step>0; CT = CTo + (dCTdt + dCTdto)/2.*dt; end                                        % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);

ct = CT./rho;


% *****  STABLE ISOTOPE  **************************************************

% reactive transfer of stable isotope ratio
trns_si = Gx.*(sim.*double(Gx<0) + six.*double(Gx>=0));

% update stable isotope ratio in melt
advn_si = advection(rho.*m.*sim,Um,Wm,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(sim,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(sim,h);
diff_si(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-sim).*rho.*m./tau_a .* bndshape; end % impose wall assimilation

dSImdt = - advn_si + diff_si + assim - trns_si;                            % total rate of change

if step>0; SIm = SImo + (dSImdt + dSImdto)/2.*dt; end               % explicit update
SIm = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho.*m,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho.*m,SIm));
SIm(m < TINY)  = 0;
SIm([1 end],:) = SIm([2 end-1],:);                                         % boundary conditions
SIm(:,[1 end]) = SIm(:,[2 end-1]);

sim = SIm./rho./max(TINY,m);

% update stable isotope ratio in xtals
advn_si = advection(rho.*x.*six,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(six,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(six,h);
diff_si(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion when melt present
                           - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-six).*rho.*x./tau_a .* bndshape; end % impose wall assimilation

dSIxdt = - advn_si + diff_si + assim + trns_si;                            % total rate of change

if step>0; SIx = SIxo + (dSIxdt + dSIxdto)/2.*dt; end                                    % explicit update
SIx = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho.*x,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho.*x,SIx));
SIx(x < TINY)  = 0;
SIx([1 end],:) = SIx([2 end-1],:);                                         % boundary conditions
SIx(:,[1 end]) = SIx(:,[2 end-1]);

six = SIx./rho./max(TINY,x);

SI = SIm + SIx;
SI = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho,SI));
si = SI./rho./max(TINY,1-f);


