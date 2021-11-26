%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase densities
rhom = rhom0 .* (1 - aTm.*(T-Tphs0) - gCm.*(cm-cphs0));
rhox = rhox0 .* (1 - aTx.*(T-Tphs0) - gCx.*(cx-cphs0));
rhof = rhof0 .* (1 - aTf.*(T-Tphs0) + bPf.*(Pt-Ptop ));

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
rhoBF  = (rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))/2 - rhoref;           % relative density for bouancy force term
Pt     = Ptop + rhoref.*g0.*ZZ;% + P;                                        % total pressure

% update thermal properties
rhoCp  = mu.*rhom.*Cpm + chi.*rhox.*Cpx + phi.*rhof.*Cpf;                  % magma heat capacity density
kT     = mu.*kTm + chi.*kTx + phi.*kTf;                                    % magma thermal conductivity

% update effective viscosity
etam  = etam0 .* exp(Em./(8.3145.*(T+273.15))-Em./(8.3145.*((Tphs0+Tphs1)/2+273.15))) ...
              .* Fmc.^(cm-(cphs0+cphs1)/2) .* Fmv.^(vm./0.01) ...
              .* max(1e-12,1-min(1-TINY,phi/phic)).^-A .* max(1e-12,1-min(1-TINY,chi/chic)).^-B; % T-c-v-chi-phi-dep. melt viscosity
etaf  = etaf0.* ones(size(f));                                             % constant volatile fluid viscosity
eta   = 1./(1./(etam + etaf) + 1./etax0);                                   % effective magma viscosity
eta   = 1./(1./(eta + etamin) + 1./etamax);                                % limit viscosity range

etac  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                       % viscosity in cell corners
      +  eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;
  
% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                           % get velocity divergence
                       + ddx(U(2:end-1,:),h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update strain rates
exx(:,2:end-1) = diff(U,1,2)./h - Div_V(:,2:end-1)./3;                     % x-normal strain rate
exx([1 end],:) = exx([2 end-1],:);                                         % apply boundary conditions
exx(:,[1 end]) = exx(:,[2 end-1]);
ezz(2:end-1,:) = diff(W,1,1)./h - Div_V(2:end-1,:)./3;                     % z-normal strain rate
ezz([1 end],:) = ezz([2 end-1],:);                                         % apply boundary conditions
ezz(:,[1 end]) = ezz(:,[2 end-1]);
exz            = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                     % shear strain rate

% update stresses
txx = eta .* exx;                                                          % x-normal stress
tzz = eta .* ezz;                                                          % z-normal stress
txz = etac.* exz;                                                          % xz-shear stress

% update tensor magnitudes
eII(2:end-1,2:end-1) = 1e-16 + (0.5.*(exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...  % get strain rate magnitude
                             + 2.*(exz(1:end-1,1:end-1).^2.*exz(2:end,1:end-1).^2.*exz(1:end-1,2:end).^2.*exz(2:end,2:end).^2).^0.25)).^0.5;
eII(:,[1 end]) = eII(:,[2 end-1]);
eII([1 end],:) = eII([2 end-1],:);

tII(2:end-1,2:end-1) = 1e-16 + (0.5.*(txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...  % get stress magnitude
                             + 2.*(txz(1:end-1,1:end-1).^2.*txz(2:end,1:end-1).^2.*txz(1:end-1,2:end).^2.*txz(2:end,2:end).^2).^0.25)).^0.5;  
tII(:,[1 end]) = tII(:,[2 end-1]);
tII([1 end],:) = tII([2 end-1],:);

% update phase segregation speeds
if bndmode==3; sds = -1;      % no slip
else;          sds = +1; end  % free slip

wx = 2/9   .* ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2./((eta(1:end-1,:)+eta(2:end,:))/2); % crystal settling speed
wx([1 end],:) = 0;
wx(:,[1 end]) = sds*wx(:,[2 end-1]);

wf = 2/9   .* ((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*df^2./((eta(1:end-1,:)+eta(2:end,:))/2) ...  % bubble flotation speed
   + 1/500 .* ((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2.*((phi(1:end-1,:)+phi(2:end,:))/2).^2.*((chi(1:end-1,:)+chi(2:end,:))/2).^2./((etaf(1:end-1,:)+etaf(2:end,:))/2); % fluid percolation speed
wf([1 end],:) = 0;
wf(:,[1 end]) = sds*wf(:,[2 end-1]);

wm = 1/50  .* ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2.*((mu (1:end-1,:)+mu (2:end,:))/2).^2.*((chi(1:end-1,:)+chi(2:end,:))/2).^2./((etam(1:end-1,:)+etam(2:end,:))/2); % melt percolation speed
wm([1 end],:) = 0;
wm(:,[1 end]) = sds*wm(:,[2 end-1]);

% update phase velocities
Wf   = W + wf;                                                             % mvp z-velocity
Uf   = U + 0.;                                                             % mvp x-velocity
Wx   = W + wx;                                                             % xtl z-velocity
Ux   = U + 0.;                                                             % xtl x-velocity
Wm   = W + wm;                                                             % mlt z-velocity
Um   = U + 0.;                                                             % mlt x-velocity

% update volume source
Div_rhoV =  + advection(rhom.*mu ,0.*U,wm,h,ADVN,'flx') ...
            + advection(rhof.*phi,0.*U,wf,h,ADVN,'flx') ...
            + advection(rhox.*chi,0.*U,wx,h,ADVN,'flx') ...
            + advection(rho      ,U   ,W ,h,ADVN,'flx');
VolSrc = -((rho-rhoo)./dt + (Div_rhoV - rho.*Div_V + Div_rhoVo)/2)./(rho/2);

dVoldt = mean(mean(VolSrc(2:end-1,2:end-1)));
UBG    = - 0.*dVoldt./2 .* (L/2-XXu);
WBG    = - 0.*dVoldt./2 .* (D/2-ZZw);

VolSrc = VolSrc - dVoldt;
