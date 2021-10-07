%% get auxiliary variables

% convert weight to volume fraction, update bulk density
res = 1;    tol = 1e-15;
while res > tol
    chii = chi;  phii = phi;
    chi  = max(TINY,min(1-TINY, x.*rho./rhox ))./(chi+phi+mu);
    phi  = max(TINY,min(1-TINY, f.*rho./rhof ))./(chi+phi+mu);
    mu   = max(TINY,min(1-TINY, m.*rho./rhom ))./(chi+phi+mu);
    rho  = mu.*rhom + chi.*rhox + phi.*rhof;
    res  = (norm(chi(:)-chii(:),2) + norm(phi(:)-phii(:),2))./sqrt(2*length(chi(:)));
end
rhoBF  = (rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))./2 - rhoref;          % relative density for bouancy force term
Pt     = Ptop + rhoref.*g0.*ZZ + 0.*P;                                     % total pressure

% update effective viscosity
eta   = etam .* max(1e-6,1-min(1-TINY,phi/0.5)).^-A .* max(1e-6,1-min(1-TINY,chi/0.5)).^-B;        % bubble-crystal-dep. magma viscosity
eta   = 1./(1./(eta + etaf) + 1./etax./exp(-30.*min(0.2,1-chi)));
eta   = 1./(1./(eta + geomean(eta(:))./sqrt(etactr)) + 1./(geomean(eta(:)).*sqrt(etactr)));
% for d = 1:ceil(1)
%     dd  = delta/ceil(1);
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

wf = 2/9   .* (rhof-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*df^2./((eta(1:end-1,:)+eta(2:end,:))/2) ...  % bubble flotation speed
   + 1/250 .* (rhof-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2.*((phi(1:end-1,:)+phi(2:end,:))/2).^2./etaf; % fluid percolation speed
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

wm = 1/50 .* (rhom-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2.*(mu(1:end-1,:)+mu(2:end,:))/2.^2./etam; % melt percolation speed
wm([1 end],:) = 0;
wm(:,[1 end]) = sds*wm(:,[2 end-1]);
for d = 1:delta
    wm(2:end-1,2:end-1) = wm(2:end-1,2:end-1) + diff(wm(:,2:end-1),2,1)./8 + diff(wm(2:end-1,:),2,2)./8;
    wm([1 end],:) = 0;
    wm(:,[1 end]) = sds*wm(:,[2 end-1]);
end

% update phase velocities
Wf   = W + wf;                                                             % mvp z-velocity
Uf   = U + 0.;                                                             % mvp x-velocity
Wx   = W + wx;                                                             % xtl z-velocity
Ux   = U + 0.;                                                             % xtl x-velocity
Wm   = W + wm;                                                             % mlt z-velocity
Um   = U + 0.;                                                             % mlt x-velocity

% rhofWf = (rho(1:end-1,:)+rho(2:end,:))/2.*(f(1:end-1,:)+f(2:end,:))/2.*Wf; % mvp z-velocity
% rhofUf = (rho(:,1:end-1)+rho(:,2:end))/2.*(f(:,1:end-1)+f(:,2:end))/2.*Uf; % mvp x-velocity
% rhoxWx = (rho(1:end-1,:)+rho(2:end,:))/2.*(x(1:end-1,:)+x(2:end,:))/2.*Wx; % xtl z-velocity
% rhoxUx = (rho(:,1:end-1)+rho(:,2:end))/2.*(x(:,1:end-1)+x(:,2:end))/2.*Ux; % xtl x-velocity
% rhomWm = (rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2.*Wm; % mlt z-velocity
% rhomUm = (rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2.*Um; % mlt x-velocity
%     
% Div_rhofVf(2:end-1,2:end-1) = ddz(rhofWf(:,2:end-1),h) + ddx(rhofUf(2:end-1,:),h); % segregation velocity divergence
% Div_rhofVf([1 end],:) = Div_rhofVf([2 end-1],:);
% Div_rhofVf(:,[1 end]) = Div_rhofVf(:,[2 end-1]);
% 
% Div_rhoxVx(2:end-1,2:end-1) = ddz(rhoxWx(:,2:end-1),h) + ddx(rhoxUx(2:end-1,:),h); % segregation velocity divergence
% Div_rhoxVx([1 end],:) = Div_rhoxVx([2 end-1],:);
% Div_rhoxVx(:,[1 end]) = Div_rhoxVx(:,[2 end-1]);
% 
% Div_rhomVm(2:end-1,2:end-1) = ddz(rhomWm(:,2:end-1),h) + ddx(rhomUm(2:end-1,:),h); % segregation velocity divergence
% Div_rhomVm([1 end],:) = Div_rhomVm([2 end-1],:);
% Div_rhomVm(:,[1 end]) = Div_rhomVm(:,[2 end-1]);

% drhodt =  + advection(rhom.*mu ,Um,Wm,h,ADVN,'flx') ...
%           + advection(rhof.*phi,Uf,Wf,h,ADVN,'flx') ...
%           + advection(rhox.*chi,Ux,Wx,h,ADVN,'flx');
Div_rhov =  + advection(rhom.*mu ,0.*Um,wm,h,ADVN,'flx') ...
            + advection(rhof.*phi,0.*Uf,wf,h,ADVN,'flx') ...
            + advection(rhox.*chi,0.*Ux,wx,h,ADVN,'flx') ...
            + advection(rho      ,U    ,W ,h,ADVN,'adv');
if step>0
VolSrci = - ((rho-rhoo)./dt + Div_rhov)./rho;
% VolSrci = - ((rho-rhoo)./dt + theta.*drhodt + (1-theta).*drhodto)./rho + Div_V;
% VolSrci = - ((rho-rhoo)./dt + Div_rhofVf+Div_rhoxVx+Div_rhomVm)./rho + Div_V;
VolSrc = beta.*VolSrc + (1-beta).*VolSrci;
end
dVoldt = mean(mean(VolSrc(2:end-1,2:end-1)));
UBG    = - dVoldt./2 .* (L/2-XXu);
WBG    = - dVoldt./2 .* (D/2-ZZw);

rhoCp  = mu*rhom*Cpm + chi*rhox*Cpx + phi*rhof*Cpf;                        % bulk heat capacity density

kT     = mu.*kTm + chi.*kTx + phi.*kTf;                                     % bubble-crystal-dep. magma thermal conductivity