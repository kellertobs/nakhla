%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase oxide compositions
wt0 = (cal.perCm-cx)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cx)./(cal.cphs1-cal.perCm);
cx_oxds = reshape((wt0(:) .* cal.oxds(1,:) + (1-wt0(:)) .* cal.oxds(2,:)) .* (cx(:)<=cal.perCm) ...
                + (wt1(:) .* cal.oxds(2,:) + (1-wt1(:)) .* cal.oxds(3,:)) .* (cx(:)> cal.perCm) ...
                  ,Nz,Nx,length(cal.oxds));
wt0 = (cal.perCm-cm)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm)./(cal.cphs1-cal.perCm);
cm_oxds = reshape((wt0(:) .* cal.oxds(1,:) + (1-wt0(:)) .* cal.oxds(2,:)) .* (cm(:)<=cal.perCm) ...
                + (wt1(:) .* cal.oxds(2,:) + (1-wt1(:)) .* cal.oxds(3,:)) .* (cm(:)> cal.perCm) ...
                  ,Nz,Nx,length(cal.oxds));
wt0 = (cal.perCm-c)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-c)./(cal.cphs1-cal.perCm);
c_oxds = reshape((wt0(:) .* cal.oxds(1,:) + (1-wt0(:)) .* cal.oxds(2,:)) .* (c(:)<=cal.perCm) ...
               + (wt1(:) .* cal.oxds(2,:) + (1-wt1(:)) .* cal.oxds(3,:)) .* (c(:)> cal.perCm) ...
                 ,Nz,Nx,length(cal.oxds));

% update phase densities
rhom = rhom0 .* (1 - aT.*(T-cal.perT-273.15) - gC.*(cm-(cal.perCx+cal.perCm)/2));
rhox = rhox0 .* (1 - aT.*(T-cal.perT-273.15) - gC.*(cx-(cal.perCx+cal.perCm)/2));
rhof = rhof0 .* (1 - aT.*(T-cal.perT-273.15) + bP.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho  = 1./(m./rhom + x./rhox + f./rhof);  

chi   = x.*rho./rhox;
phi   = f.*rho./rhof;
mu    = m.*rho./rhom;                                  

% update effective viscosity
wtm      = zeros(Nz*Nx,12);
wtm(:, 1) = reshape(cm_oxds(:,:,1),Nz*Nx,1); % SiO2
wtm(:, 3) = reshape(cm_oxds(:,:,2),Nz*Nx,1); % Al2O3
wtm(:, 4) = reshape(cm_oxds(:,:,3),Nz*Nx,1); % FeO
wtm(:, 6) = reshape(cm_oxds(:,:,4),Nz*Nx,1); % MgO
wtm(:, 7) = reshape(cm_oxds(:,:,5),Nz*Nx,1); % CaO
wtm(:, 8) = reshape(cm_oxds(:,:,6),Nz*Nx,1); % Na2O
wtm(:, 9) = reshape(cm_oxds(:,:,7),Nz*Nx,1); % K2O
wtm(:,11) = reshape(100.*vm(:,:  ),Nz*Nx,1); % H2O
etam      = reshape(grdmodel08(wtm,T(:)-273.15),Nz,Nx);

etaf  = etaf0.* ones(size(f));                                             % constant fluid viscosity
etax  = etax0.* ones(size(x));                                             % constant crysta viscosity

% get permission weights
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-9,min(1-1e-9,permute(cat(3,chi,mu,phi),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get effective viscosity
eta    = squeeze(sum(ff.*kv.*thtv,1));  if size(eta,1)~=size(T,1); eta = eta.'; end
if ~calibrt; etamax = 1e+6.*min(eta(:)); else; etamax = 1e3.*etax0; end
eta    = (1./etamax + 1./(eta.*etareg)).^-1;
etaco  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...
       +  eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;

% get segregation coefficients
Csgr = ((1-ff)./[dx;dm;df].^2.*kv.*thtv).^-1 + 1e-18;

Csgr_x = squeeze(Csgr(1,:,:)); if size(Csgr_x,1)~=size(T,1); Csgr_x = Csgr_x.'; end
Csgr_f = squeeze(Csgr(3,:,:)); if size(Csgr_f,1)~=size(T,1); Csgr_f = Csgr_f.'; end
Csgr_m = squeeze(Csgr(2,:,:)); if size(Csgr_m,1)~=size(T,1); Csgr_m = Csgr_m.'; end
Csgr_m = Csgr_m.*chi.^2; % dampen melt segregation at low crystallinity

if ~calibrt % skip the following if called from calibration script

% diffusion parameters
ks = kT./T;
kc = kT./cP./10.*mu.^0.5;
kx = kT./cP./10.*mu.^0.5.*x;
kf = kT./cP./10.*mu.^0.5.*f;
km = kT./cP./10.*mu.^0.5.*m;

% update phase segregation speeds
wm = (rhom-rho).*g0.*Csgr_m;
wm = sign((wm(1:end-1,:)+wm(2:end,:))/2).*min(abs(wm(1:end-1,:)),abs(wm(2:end,:)));
% wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_m(1:end-1,:)+1./Csgr_m(2:end,:)); % melt segregation speed
wm(1  ,:)     = min(1,1-top).*wm(1  ,:);
wm(end,:)     = min(1,1-bot).*wm(end,:);
wm(:,[1 end]) = -sds*wm(:,[2 end-1]);

wx = (rhox-rho).*g0.*Csgr_x;
wx = sign((wx(1:end-1,:)+wx(2:end,:))/2).*min(abs(wx(1:end-1,:)),abs(wx(2:end,:)));
% wx = (((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0).*2./(1./Csgr_x(1:end-1,:)+1./Csgr_x(2:end,:)); % crystal segregation speed
wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
wx(end,:)     = min(1,1-bot).*wx(end,:);
wx(:,[1 end]) = -sds*wx(:,[2 end-1]);

wf = (rhof-rho).*g0.*Csgr_f;
wf = sign((wf(1:end-1,:)+wf(2:end,:))/2).*min(abs(wf(1:end-1,:)),abs(wf(2:end,:)));
% wf = (((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0).*2./(1./Csgr_f(1:end-1,:)+1./Csgr_f(2:end,:)); % fluid segregation speed
wf(1  ,:)     = min(1,1-top+fout).*wf(1  ,:);
wf(end,:)     = min(1,1-bot+fin ).*wf(end,:);
wf(:,[1 end]) = -sds*wf(:,[2 end-1]);


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
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

% update tensor magnitudes
eII(2:end-1,2:end-1) = 1e-16 + (0.5.*(exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...  % get strain rate magnitude
                             + 2.*(exz(1:end-1,1:end-1).^2.*exz(2:end,1:end-1).^2.*exz(1:end-1,2:end).^2.*exz(2:end,2:end).^2).^0.25)).^0.5;
eII(:,[1 end]) = eII(:,[2 end-1]);
eII([1 end],:) = eII([2 end-1],:);

tII(2:end-1,2:end-1) = 1e-16 + (0.5.*(txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...  % get stress magnitude
                             + 2.*(txz(1:end-1,1:end-1).^2.*txz(2:end,1:end-1).^2.*txz(1:end-1,2:end).^2.*txz(2:end,2:end).^2).^0.25)).^0.5;  
tII(:,[1 end]) = tII(:,[2 end-1]);
tII([1 end],:) = tII([2 end-1],:);

% heat dissipation (entropy production) rate
[grdTx,grdTz] = gradient(T,h);
diss =  exx(2:end-1,2:end-1).*txx(2:end-1,2:end-1) ...
     +  ezz(2:end-1,2:end-1).*tzz(2:end-1,2:end-1) ...
     +  2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
         .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
     +  mu (2:end-1,2:end-1)./Csgr_m(2:end-1,2:end-1) .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
     +  chi(2:end-1,2:end-1)./Csgr_x(2:end-1,2:end-1) .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
     +  phi(2:end-1,2:end-1)./Csgr_f(2:end-1,2:end-1) .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2 ...
     +  ks(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2);

% update volume source
Div_rhoV =  + advection(rho.*f,0.*U,wf,h,ADVN,'flx') ...
            + advection(rho.*x,0.*U,wx,h,ADVN,'flx') ...
            + advection(rho.*m,0.*U,wm,h,ADVN,'flx') ...
            + advection(rho   ,U   ,W ,h,ADVN,'flx');
if step>0; VolSrc  = -((rho-rhoo)./dt + Div_rhoV - rho.*Div_V)./rho; end

UBG    = - 2*mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (L/2-XXu);
WBG    = - 0*mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (D/2-ZZw);
end