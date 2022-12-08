%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update oxide compositions
wt0 = (cal.perCm-cm)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm)./(cal.cphs1-cal.perCm);
cm_oxd = reshape((wt0(:) .* cal.cmp_oxd(1,:) + (1-wt0(:)) .* cal.cmp_oxd(3,:)) .* (cm(:)<=cal.perCm) ...
               + (wt1(:) .* cal.cmp_oxd(3,:) + (1-wt1(:)) .* cal.cmp_oxd(4,:)) .* (cm(:)> cal.perCm) ...
                 ,Nz,Nx,cal.nc);

wt0 = (cal.perCx-cx)./(cal.perCx-cal.cphs0);
wt1 = (cal.cphs1-cx)./(cal.cphs1-cal.perCx);
cx_oxd = reshape((wt0(:) .* cal.cmp_oxd(1,:) + (1-wt0(:)) .* cal.cmp_oxd(2,:)) .* (cx(:)<=cal.perCx) ...
               + (wt1(:) .* cal.cmp_oxd(2,:) + (1-wt1(:)) .* cal.cmp_oxd(4,:)) .* (cx(:)> cal.perCx) ...
                 ,Nz,Nx,cal.nc);

c_oxd = (m.*cm_oxd + x.*cx_oxd)./(1-f);

cm_cmp = reshape(reshape(cm_oxd,Nz*Nx,cal.nc)/cal.oxd*100,Nz,Nx,cal.nc);
cx_cmp = reshape(reshape(cx_oxd,Nz*Nx,cal.nc)/cal.oxd*100,Nz,Nx,cal.nc);
 c_cmp = reshape(reshape( c_oxd,Nz*Nx,cal.nc)/cal.oxd*100,Nz,Nx,cal.nc);

% update phase densities
rhom = squeeze(sum(permute(cm_cmp/100,[3,1,2])./cal.rhom0.')).^-1 .* (1 - cal.aT.*(T-cal.perT-273.15) - cal.gH.*vm); if size(rhom,2)~=size(T,2); rhom = rhom(1,:).'; end
rhox = squeeze(sum(permute(cx_cmp/100,[3,1,2])./cal.rhox0.')).^-1 .* (1 - cal.aT.*(T-cal.perT-273.15)             ); if size(rhox,2)~=size(T,2); rhox = rhox(1,:).'; end
rhof = cal.rhof0 .* (1 - cal.aT.*(T-cal.perT-273.15) + cal.bP.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox + f./rhof);

rhofz = (rho(1:end-1,:)+rho(2:end,:))/2;

chi   = max(TINY,min(1-TINY, x.*rho./rhox ));
phi   = max(TINY,min(1-TINY, f.*rho./rhof ));
mu    = max(TINY,min(1-TINY, m.*rho./rhom ));                                  

% update effective viscosity
wtm      = zeros(Nz*Nx,12);
wtm(:, 1) = reshape(cm_oxd(:,:,1),Nz*Nx,1); % SiO2
wtm(:, 2) = reshape(cm_oxd(:,:,2),Nz*Nx,1); % TiO2
wtm(:, 3) = reshape(cm_oxd(:,:,3),Nz*Nx,1); % Al2O3
wtm(:, 4) = reshape(cm_oxd(:,:,4),Nz*Nx,1); % FeO
wtm(:, 6) = reshape(cm_oxd(:,:,5),Nz*Nx,1); % MgO
wtm(:, 7) = reshape(cm_oxd(:,:,6),Nz*Nx,1); % CaO
wtm(:, 8) = reshape(cm_oxd(:,:,7),Nz*Nx,1); % Na2O
wtm(:, 9) = reshape(cm_oxd(:,:,8),Nz*Nx,1); % K2O
wtm(:,11) = reshape(100.*vm(:,: ),Nz*Nx,1); % H2O
etam      = reshape(grdmodel08(wtm,T(:)-273.15),Nz,Nx);

etaf = cal.etaf0.* ones(size(f));  % constant fluid viscosity
etax = cal.etax0.* ones(size(x));  % constant solid viscosity

% get permission weights
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(TINY,min(1-TINY,permute(cat(3,chi,mu,phi),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./cal.BB).^(1./cal.CC);  Sf = Sf./sum(Sf,2);
Xf = sum(cal.AA.*Sf,2).*FF + (1-sum(cal.AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get effective viscosity
eta    = squeeze(sum(ff.*kv.*thtv,1));  if size(eta,1)~=size(T,1); eta = eta.'; end
if ~calibrt; etamax = 1e+6.*min(eta(:)); else; etamax = 1e3.*cal.etax0; end
eta    = (1./etamax + 1./(eta*etareg)).^-1;
etaco  = (eta(1:end-1,1:end-1).*eta(2:end,1:end-1) ...
       .* eta(1:end-1,2:end  ).*eta(2:end,2:end  )).^0.25;

% get segregation coefficients
dd = permute(cat(3,d0*ones(size(mu)),d0*(1-mu),d0*ones(size(mu))),[3,1,2]);
Csgr = ((1-ff)./dd.^2.*kv.*thtv).^-1 + TINY^2;

Csgr_x = squeeze(Csgr(1,:,:)); if size(Csgr_x,1)~=size(T,1); Csgr_x = Csgr_x.'; end
Csgr_f = squeeze(Csgr(3,:,:)); if size(Csgr_f,1)~=size(T,1); Csgr_f = Csgr_f.'; end
Csgr_m = squeeze(Csgr(2,:,:)); if size(Csgr_m,1)~=size(T,1); Csgr_m = Csgr_m.'; end
% Csgr_m = Csgr_m.*(1-mu).^2 + TINY^2; % dampen melt segregation at high melt fraction (dm = d0.*(1-f))

if ~calibrt % skip the following if called from calibration script

wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*(Csgr_m(1:end-1,:).*Csgr_m(2:end,:)).^0.5; % melt segregation speed
% wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_m(1:end-1,:)+1./Csgr_m(2:end,:)); % melt segregation speed
wm(1  ,:)     = min(1,1-top).*wm(1  ,:);
wm(end,:)     = min(1,1-bot).*wm(end,:);
wm(:,[1 end]) = -sds*wm(:,[2 end-1]);

wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*(Csgr_x(1:end-1,:).*Csgr_x(2:end,:)).^0.5; % solid segregation speed
% wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_x(1:end-1,:)+1./Csgr_x(2:end,:)); % solid segregation speed
wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
wx(end,:)     = min(1,1-bot).*wx(end,:);
wx(:,[1 end]) = -sds*wx(:,[2 end-1]);

wf = ((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*(Csgr_f(1:end-1,:).*Csgr_f(2:end,:)).^0.5; % fluid segregation speed
% wf = ((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_f(1:end-1,:)+1./Csgr_f(2:end,:)); % fluid segregation speed
wf(1  ,:)     = min(1,1-top+fout).*wf(1  ,:);
wf(end,:)     = min(1,1-bot+fin ).*wf(end,:);
wf(:,[1 end]) = -sds*wf(:,[2 end-1]);

% diffusion parameters
kW  = Vel*h/100;                                                           % convection fluctuation diffusivity
kwm = abs((rhom-rho).*g0.*Csgr_m*d0*10);                                   % segregation fluctuation diffusivity
kwx = abs((rhox-rho).*g0.*Csgr_x*d0*10);                                   % segregation fluctuation diffusivity
kwf = abs((rhof-rho).*g0.*Csgr_f*d0*10);                                   % segregation fluctuation diffusivity
km  = m.*(kwm + kW).*rho;                                                  % melt  fraction diffusion 
kx  = x.*(kwx + kW).*rho;                                                  % solid fraction diffusion 
kf  = f.*(kwf + kW).*rho;                                                  % fluid fraction diffusion 
kc  = (x.*kwx + f.*kwf + m.*kwm + kW).*rho;                                % component diffusion
kT  = kT0 + kc.*cP;                                                        % heat diffusion
ks  = kT./T;                                                               % entropy diffusion

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
     +   mu(2:end-1,2:end-1)./Csgr_m(2:end-1,2:end-1) .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
     +  chi(2:end-1,2:end-1)./Csgr_x(2:end-1,2:end-1) .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
     +  phi(2:end-1,2:end-1)./Csgr_f(2:end-1,2:end-1) .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2 ...
     +  ks(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2);


% update volume source
if step>0
    Div_rhoV =  + advect(M(inz,inx),0.*U(inz,:),wm(:,inx),h,{ADVN,''   },[1,2],BCA) ...
                + advect(X(inz,inx),0.*U(inz,:),wx(:,inx),h,{ADVN,''   },[1,2],BCA) ...
                + advect(F(inz,inx),0.*U(inz,:),wf(:,inx),h,{ADVN,''   },[1,2],BCA) ...
                + advect(rho(inz,inx), U(inz,:), W(:,inx),h,{ADVN,'vdf'},[1,2],BCA);
    VolSrc = -((rho(inz,inx)-rhoo(inz,inx))./dt + Div_rhoV)./rho(inz,inx); 
%     VolSrc = -((rho(inz,inx)-rhoo(inz,inx))./dt + theta.*Div_rhoV + (1-theta).*Div_rhoVo)./rho(inz,inx);
end

UBG    = - 1*mean(mean(VolSrc))./2 .* (L/2-XXu);
WBG    = - 1*mean(mean(VolSrc))./2 .* (D/2-ZZw);

UDtime = UDtime + toc;
end