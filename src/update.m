%% get auxiliary variables

% convert weight to volume fraction, update bulk density
res = 1;    tol = 1e-15;
while res > tol
    chii = chi;  phii = phi;
    chi  = max(TINY,min(1-TINY, x.*rho./rhox ));
    phi  = max(TINY,min(1-TINY, f.*rho./rhof ));
    mu   = 1-chi-phi;
    rho  = mu.*rhom + chi.*rhox + phi.*rhof;
    res  = (norm(chi(:)-chii(:),2) + norm(phi(:)-phii(:),2))./sqrt(2*length(chi(:)));
end
rhoref = mean(mean(rho(2:end-1,2:end-1)));                                 % reference density for magmastatic pressure
rhoBF  = (rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))./2 - rhoref;          % relative density for bouancy force term
Pt     = Ptop + P + rhoref .* g0 .* ZZ;                                        % total pressure

% update effective viscosity
eta   = eta0 .* max(1e-6,1-phi/0.5).^-A .* max(1e-6,1-chi/0.5).^-B;        % bubble-crystal-dep. magma viscosity
eta   = 1./(1./(eta + eta0./sqrt(etactr)) + 1./(eta0.*sqrt(etactr)));
% for d = 1:ceil(delta)
%     dd  = delta/ceil(delta);
%     eta = eta + dd.*(diff(eta([1,1:end,end],:),2,1)./8 + diff(eta(:,[1,1:end,end]),2,2)./8);
% end
etac  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                       % viscosity in cell corners
      +  eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;
  
% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                           % get velocity divergence
                       + ddx(U(2:end-1,:),h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update strain rates
exx(:,2:end-1)   = diff(U,1,2)./h - Div_V(:,2:end-1)./3;                   % x-normal strain rate
exx([1 end],:)   = exx([2 end-1],:);                                       % apply boundary conditions
exx(:,[1 end])   = exx(:,[2 end-1]);
ezz(2:end-1,:)   = diff(W,1,1)./h - Div_V(2:end-1,:)./3;                   % z-normal strain rate
ezz([1 end],:)   = ezz([2 end-1],:);                                       % apply boundary conditions
ezz(:,[1 end])   = ezz(:,[2 end-1]);
exz              = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                   % shear strain rate

% update stresses
txx = eta .* exx;                                                          % x-normal stress
tzz = eta .* ezz;                                                          % z-normal stress
txz = etac.* exz;                                                          % xz-shear stress

% update solution-dependent parameter and auxiliary variable fields
eII(2:end-1,2:end-1) = 1e-16 + (0.5.*(exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...  % get strain rate magnitude
    + 2.*(exz(1:end-1,1:end-1).^2.*exz(2:end,1:end-1).^2.*exz(1:end-1,2:end).^2.*exz(2:end,2:end).^2).^0.25)).^0.5;
eII(:,[1 end]) = eII(:,[2 end-1]);
eII([1 end],:) = eII([2 end-1],:);

tII(2:end-1,2:end-1) = 1e-16 + (0.5.*(txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...  % get stress magnitude
    + 2.*(txz(1:end-1,1:end-1).^2.*txz(2:end,1:end-1).^2.*txz(1:end-1,2:end).^2.*txz(2:end,2:end).^2).^0.25)).^0.5;  
tII(:,[1 end]) = tII(:,[2 end-1]);
tII([1 end],:) = tII([2 end-1],:);

% update phase segregation speeds
if coolmode==3; sds = -1;      % no slip
else;           sds = +1; end  % free slip
wf = 2/9 .* (rhof-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*df^2./((eta(1:end-1,:)+eta(2:end,:))/2); % bubble flotation speed
wf([1 end],:) = 0;
wf(:,[1 end]) = sds*wf(:,[2 end-1]);
for d = 1:delta
    wf(2:end-1,2:end-1) = wf(2:end-1,2:end-1) + diff(wf(:,2:end-1),2,1)./8 + diff(wf(2:end-1,:),2,2)./8;
    wf([1 end],:) = 0;
    wf(:,[1 end]) = sds*wf(:,[2 end-1]);
end

wx = 2/9 .* (rhox-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2./((eta(1:end-1,:)+eta(2:end,:))/2); % crystal settling speed
wx([1 end],:) = 0;
wx(:,[1 end]) = sds*wx(:,[2 end-1]);
for d = 1:delta
    wx(2:end-1,2:end-1) = wx(2:end-1,2:end-1) + diff(wx(:,2:end-1),2,1)./8 + diff(wx(2:end-1,:),2,2)./8;
    wx([1 end],:) = 0;
    wx(:,[1 end]) = sds*wx(:,[2 end-1]);
end

% update phase velocities
Wf   = W + wf;                                                             % mvp z-velocity
Uf   = U + 0.;                                                             % mvp x-velocity
Wx   = W + wx;                                                             % xtl z-velocity
Ux   = U + 0.;                                                             % xtl x-velocity
Wm   = W;                                                                  % mlt z-velocity
Um   = U;                                                                  % mlt x-velocity

fWf  = (f(1:end-1,:)+f(2:end,:))/2.*(W + wf);                              % mvp z-velocity
fUf  = (f(:,1:end-1)+f(:,2:end))/2.*(U + 0.);                              % mvp x-velocity
xWx  = (x(1:end-1,:)+x(2:end,:))/2.*(W + wx);                              % xtl z-velocity
xUx  = (x(:,1:end-1)+x(:,2:end))/2.*(U + 0.);                              % xtl x-velocity
mUm  = (m(:,1:end-1)+m(:,2:end))/2.*(U + 0.);                              % mlt x-velocity
mWm  = (m(1:end-1,:)+m(2:end,:))/2.*(W + 0.);                              % mlt z-velocity

% WT  = W + ((f(1:end-1,:)+f(2:end,:))/2.*rhof.*Cf.*wf ...                   % heat transport z-velocity
%         +  (x(1:end-1,:)+x(2:end,:))/2.*rhox.*Cx.*wx) ...
%         ./((rho(1:end-1,:)+rho(2:end,:))./2.*(Cp(1:end-1,:)+Cp(2:end,:))./2);
    
dwfdz(2:end-1,2:end-1) = ddz((f(1:end-1,2:end-1)+f(2:end,2:end-1))/2 .* wf(:,2:end-1),h); % segregation velocity divergence
dwfdz([1 end],:) = dwfdz([2 end-1],:);
dwfdz(:,[1 end]) = dwfdz(:,[2 end-1]);

dwxdz(2:end-1,2:end-1) = ddz((x(1:end-1,2:end-1)+x(2:end,2:end-1))/2 .* wx(:,2:end-1),h); % segregation velocity divergence
dwxdz([1 end],:) = dwxdz([2 end-1],:);
dwxdz(:,[1 end]) = dwxdz(:,[2 end-1]);

dWmdz(2:end-1,2:end-1) = ddz((m(1:end-1,2:end-1)+m(2:end,2:end-1))/2 .* W(:,2:end-1),h); % segregation velocity divergence
dWmdz([1 end],:) = dWmdz([2 end-1],:);
dWmdz(:,[1 end]) = dWmdz(:,[2 end-1]);

Div_fVf(2:end-1,2:end-1) = ddz(fWf(:,2:end-1),h) + ddx(fUf(2:end-1,:),h); % segregation velocity divergence
Div_fVf([1 end],:) = Div_fVf([2 end-1],:);
Div_fVf(:,[1 end]) = Div_fVf(:,[2 end-1]);

Div_xVx(2:end-1,2:end-1) = ddz(xWx(:,2:end-1),h) + ddx(xUx(2:end-1,:),h); % segregation velocity divergence
Div_xVx([1 end],:) = Div_xVx([2 end-1],:);
Div_xVx(:,[1 end]) = Div_xVx(:,[2 end-1]);

Div_mVm(2:end-1,2:end-1) = ddz(mWm(:,2:end-1),h) + ddx(mUm(2:end-1,:),h); % segregation velocity divergence
Div_mVm([1 end],:) = Div_mVm([2 end-1],:);
Div_mVm(:,[1 end]) = Div_mVm(:,[2 end-1]);

VolSrc = ((rho-rhoo)./dt + rhof.*Div_fVf+rhox.*Div_xVx+rhom.*Div_mVm)./rho - Div_V - dwxdz - dwfdz;
VolErr = mean(mean(VolSrc(2:end-1,2:end-1)));
VolSrc = VolSrc - VolErr;

rhoCp = mu*rhom*Cpm + chi*rhox*Cpx + phi*rhof*Cpf;                        % bulk heat capacity density

kT    = mu.*kTm + chi.*kTx + phi.*kTf;                                     % bubble-crystal-dep. magma thermal conductivity