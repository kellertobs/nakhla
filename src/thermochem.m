%% update time step
dtk = min((h/2)^2./max([kT(:)./rhoCp(:);kc]));                             % diffusive time step size
divx = advection(x,Ux,Wx,h,ADVN,'flx');  divf = advection(f,Uf,Wf,h,ADVN,'flx');
dtf = min(0.01./max(abs(divx(:))),0.01./max(abs(divf(:))));              % phase update time step size
dta = min(min(h/2/max(abs([UBG(:);WBG(:);Um(:);Wm(:);Uf(:);Wf(:);Ux(:);Wx(:)]+1e-16)))); % advective time step size
dt  = min(2*dto,CFL*min([dtk,dtf,dta]));                                   % physical time step size


%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% get boundary shape function
switch bndmode
    case 0
        bndshape = zeros(size(T));
    case 1
        bndshape = exp( ( -ZZ)/dw);
    case 2
        bndshape = exp( ( -ZZ)/dw) ...
                 + exp(-(D-ZZ)/dw);
    case 3
        bndshape = exp( ( -ZZ)/dw) ...
                 + exp(-(D-ZZ)/dw) ...
                 + exp( ( -XX)/dw) ...
                 + exp(-(L-XX)/dw);
end
bndshape = min(1,bndshape);


% update temperature
if ~isotherm

    advn_H = advection(rhom.*mu .*(Cpm.*T + 0  ),Um,Wm,h,ADVN,'flx') ...
           + advection(rhox.*chi.*(Cpx.*T + DLx),Ux,Wx,h,ADVN,'flx') ...
           + advection(rhof.*phi.*(Cpf.*T + DLf),Uf,Wf,h,ADVN,'flx');
                           
    qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                 % heat diffusion z-flux
    qTx    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddx(T,h);                 % heat diffusion x-flux
    diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...                 % heat diffusion
                               - ddx(qTx(2:end-1,:),h));
    
    cool = zeros(size(T));
    if ~isnan(Twall); cool = cool + rhoCp .* (Twall-T)./tau_T .* bndshape; end  % impose wall cooling
        
    dHdt = - advn_H + diff_T + cool;                                       % total rate of change
    
    if step>0; H = Ho + (theta.*dHdt + (1-theta).*dHdto).*dt; end          % explicit update of enthalpy
    H([1 end],:) = H([2 end-1],:);                                         % apply boundary conditions
    H(:,[1 end]) = H(:,[2 end-1]);    
    
end


% update composition
if ~isochem
  
    % update major component
    advn_C = advection(rhom.*mu .*cm,Um,Wm,h,ADVN,'flx') ...
           + advection(rhox.*chi.*cx,Ux,Wx,h,ADVN,'flx');
                          
    qcz   = - kc.*(rhom(1:end-1,:)+rhom(2:end,:))/2.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(cm,h);        % major component diffusion z-flux
    qcx   = - kc.*(rhom(:,1:end-1)+rhom(:,2:end))/2.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(cm,h);        % major component diffusion x-flux
    diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                  % major component diffusion
                              - ddx(qcx(2:end-1,:),h);
    
    assim = zeros(size(c));
    if ~isnan(cwall); assim = assim + (cwall-c).*rho./tau_c .* bndshape; end % impose wall assimilation
    
    dCdt = - advn_C + diff_c + assim;                                      % total rate of change
    
    if step>0; C = Co + (theta.*dCdt + (1-theta).*dCdto).*dt; end          % explicit update of major component density
    C = min(rho.*cphs1-TINY,max(rho.*cphs0+TINY, C ));
    C([1 end],:) = C([2 end-1],:);                                         % apply boundary conditions
    C(:,[1 end]) = C(:,[2 end-1]);  
    
    % update volatile component
    advn_V = advection(rhom.*mu .*vm,Um,Wm,h,ADVN,'flx') ...
           + advection(rhof.*phi.*vf,Uf,Wf,h,ADVN,'flx');
        
    qvz   = - kc.*(rhom(1:end-1,:)+rhom(2:end,:))/2.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(vm,h);        % volatile component diffusion z-flux
    qvx   = - kc.*(rhom(:,1:end-1)+rhom(:,2:end))/2.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(vm,h);        % volatile component diffusion x-flux
    diff_v(2:end-1,2:end-1) = - ddz(qvz(:,2:end-1),h) ...                  % volatile component diffusion
                              - ddx(qvx(2:end-1,:),h);
    
    outgas = zeros(size(v));
    if ~isnan(fwall); outgas = outgas + min(0,fwall-f).*rho./tau_f.* bndshape; end % impose wall outgassing
    
    assim = zeros(size(v));
    if ~isnan(vwall); assim = assim + (vwall-v).*rho./tau_c .* bndshape; end % impose wall assimilation
    
    dVdt = - advn_V + diff_v + outgas + assim;                             % total rate of change
    
    if step>0; V = Vo + (theta.*dVdt + (1-theta).*dVdto).*dt; end          % explicit update of volatile component density
    V = min(rho.*1-TINY,max(rho.*0+TINY, V ));
    V([1 end],:) = V([2 end-1],:);                                         % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);  
    
end


% convert enthalpy and component densities to temperature and concentrations
T = (H - chi.*rhox.*DLx - phi.*rhof.*DLf)./rhoCp;
c = C./rho;
v = V./rho;



%% *****  UPDATE LOCAL PHASE EQUILIBRIUM  *********************************

[xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,(T+To)./2,(c+co)./2,(v+vo)./2,(Pt+Pto)./2,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);

    
% update crystal fraction
if diseq
    
    Gxi = (xq-x)./max(4.*dt,tau_r);
    Gx  = alpha.*Gx + (1-alpha).*Gxi;                                        % crystallisation rate
    
    advn_x = advection(x,Ux,Wx,h,ADVN,'flx');                              % get advection term
    
    dxdt   = - advn_x + Gx;                                                % total rate of change
    
    if step>0; x = x + (theta.*dxdt + (1-theta).*dxdto).*dt; end           % explicit update of crystal fraction
    x = min(1-f-TINY,max(TINY,x));                                         % enforce [0,1] limit
    x([1 end],:) = x([2 end-1],:);                                         % apply boundary conditions
    x(:,[1 end]) = x(:,[2 end-1]);
    
    cx  = cxq;
    cm  = max(c+TINY,min(cphs1-TINY,(c - x.*cx)./(1-x-f)));
    cm(m < 1e-3) = cmq(m < 1e-3);
    
else
    
    x  = alpha.*x + (1-alpha).*xq;
    cx = cxq;
    cm = cmq;
    Gx = (x-xo)./dt + advection(x,Ux,Wx,h,ADVN,'flx');                     % reconstruct crystallisation rate
    
end


% update bubble fraction
if diseq
    
    Gfi = (fq-f)./max(4.*dt,tau_r);
    Gf  = alpha.*Gf + (1-alpha).*Gfi;
    
    advn_f = advection(f,Uf,Wf,h,ADVN,'flx');                              % get advection term
    
    dfdt   = - advn_f + Gf;                                                % total rate of change
    
    if step>0; f = fo + (theta.*dfdt + (1-theta).*dfdto).*dt; end          % explicit update of bubble fraction
    f = min(1-x-TINY,max(TINY,f));                                         % enforce [0,1-x] limit
    f([1 end],:) = f([2 end-1],:);                                         % apply boundary conditions
    f(:,[1 end]) = f(:,[2 end-1]);
    
    vf  = vfq;
    vm  = max(v+TINY,min(1 -TINY,(v - f.*vf)./(1-f-x)));
    vm(m < 1e-3) = vmq(m < 1e-3);
    
else
    
    f  = alpha.*f + (1-alpha).*fq;
    vf = vfq;
    vm = vmq;
    Gf = (f-fo)./dt + advection(f,Uf,Wf,h,ADVN,'flx');                     % reconstruct exsolution rate
    
end


% update melt fraction
m = 1-x-f;



%% ***** TRACE & ISOTOPE GEOCHEMISTRY  ************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
itm = it./(m + x.*KIT);
itx = it./(m./KIT + x);

% update incompatible trace element composition
advn_IT = advection(rhom.*mu .*itm,Um,Wm,h,ADVN,'flx') ...
        + advection(rhox.*chi.*itx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rhom(1:end-1,:)+rhom(2:end,:))/2.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(itm,h);
qx   = - kc.*(rhom(:,1:end-1)+rhom(:,2:end))/2.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(itm,h);
diff_it(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(itwall); assim = assim + (itwall-it).*rho./tau_c .* bndshape; end % impose wall assimilation
    
dITdt = - advn_IT + diff_it + assim;                                       % total rate of change

if step>0; IT = ITo + (theta.*dITdt + (1-theta).*dITdto).*dt; end          % explicit update
IT = max(0+TINY, IT );
IT([1 end],:) = IT([2 end-1],:);                                           % boundary conditions
IT(:,[1 end]) = IT(:,[2 end-1]);

it = IT./rho;


% *****  COMPATIBLE TRACE ELEMENT  ****************************************

% update compatible trace element phase compositions
ctm = ct./(m + x.*KCT);
ctx = ct./(m./KCT + x);

% update compatible trace element composition
advn_CT = advection(rhom.*mu .*ctm,Um,Wm,h,ADVN,'flx') ...
        + advection(rhox.*chi.*ctx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rhom(1:end-1,:)+rhom(2:end,:))/2.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(ctm,h);
qx   = - kc.*(rhom(:,1:end-1)+rhom(:,2:end))/2.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(ctm,h);
diff_ct(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(ct));
if ~isnan(ctwall); assim = assim + (ctwall-ct).*rho./tau_c .* bndshape; end % impose wall assimilation

dCTdt = - advn_CT + diff_ct + assim;                                       % total rate of change

if step>0; CT = CTo + (theta.*dCTdt + (1-theta).*dCTdto).*dt; end          % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);

ct = CT./rho;


% *****  STABLE ISOTOPE  **************************************************

% reactive transfer of stable isotope ratio
trns_si = Gx.*rho.*(sim.*double(Gx<0) + six.*double(Gx>=0));

% update stable isotope ratio in melt
advn_si = advection(SIm,Um,Wm,h,ADVN,'flx');

qz   = - kc.*(rhom(1:end-1,:)+rhom(2:end,:))/2.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(sim,h);
qx   = - kc.*(rhom(:,1:end-1)+rhom(:,2:end))/2.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(sim,h);
diff_si(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-sim).*rhom.*mu./tau_c .* bndshape; end % impose wall assimilation

dSImdt = - advn_si + diff_si + assim - trns_si;                            % total rate of change

if step>0; SIm = SImo + (theta.*dSImdt + (1-theta).*dSImdto).*dt; end      % explicit update
SIm = min(max(si0-dsi,si1+dsi).*rhom.*mu,max(min(si0-dsi,si1+dsi).*rhom.*mu,SIm));
SIm(mu<1e-3)   = SI(mu<1e-3);
SIm([1 end],:) = SIm([2 end-1],:);                                         % boundary conditions
SIm(:,[1 end]) = SIm(:,[2 end-1]);

sim = SIm./rhom./max(TINY,mu);

% update stable isotope ratio in xtals
advn_si = advection(SIx,Ux,Wx,h,ADVN,'flx');
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-six).*rhox.*chi./tau_c .* bndshape; end % impose wall assimilation

dSIxdt = - advn_si + assim + trns_si;                             % total rate of change

if step>0; SIx = SIxo + (theta.*dSIxdt + (1-theta).*dSIxdto).*dt; end      % explicit update
SIx = min(max(si0-dsi,si1+dsi).*rhox.*chi,max(min(si0-dsi,si1+dsi).*rhox.*chi,SIx));
SIx(chi<1e-3)  = SI(chi<1e-3);
SIx([1 end],:) = SIx([2 end-1],:);                                         % boundary conditions
SIx(:,[1 end]) = SIx(:,[2 end-1]);

six = SIx./rhox./max(1e-3,chi);

SI = SIm + SIx;
SI = min(max(si0-dsi,si1+dsi).*rho.*(1-f),max(min(si0-dsi,si1+dsi).*rho.*(1-f),SI));
si = SI./rho./max(TINY,1-f);


% % update stable isotope composition
% advn_si = advection(mu .*si,Um,Wm,h,ADVN,'flx') ...
%         + advection(chi.*si,Ux,Wx,h,ADVN,'flx');
% 
% qz   = - kc.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(si,h);
% qx   = - kc.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(si,h);
% diff_si(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
%                            - ddx(qx(2:end-1,:),h);
%                        
% assim = zeros(size(si));
% if ~isnan(siwall); assim = assim + (siwall-si).*rho./tau_c .* bndshape; end % impose wall assimilation
% 
% dsidt = - advn_si + diff_si + assim;                                       % total rate of change
% 
% if step>0; si = sio + (theta.*dsidt + (1-theta).*dsidto).*dt; end          % explicit update
% si([1 end],:) = si([2 end-1],:);                                           % boundary conditions
% si(:,[1 end]) = si(:,[2 end-1]);


