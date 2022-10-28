%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update phase oxide compositions
wt0 = (cal.perCx-cx)./(cal.perCx-cal.cphs0);
wt1 = (cal.perCm-cx)./(cal.perCm-cal.perCx);
wt2 = (cal.cphs1-cx)./(cal.cphs1-cal.perCm);
cx_cmps = reshape((wt0(:) .* cal.cmps(1,:) + (1-wt0(:)) .* cal.cmps(2,:)) .* (cx(:)<=cal.perCx) ...
                + (wt1(:) .* cal.cmps(2,:) + (1-wt1(:)) .* cal.cmps(3,:)) .* (cx(:)> cal.perCx & cx(:)<=cal.perCm) ...
                + (wt2(:) .* cal.cmps(3,:) + (1-wt2(:)) .* cal.cmps(4,:)) .* (cx(:)> cal.perCm) ...
                  ,Nz,Nx,length(cal.cmps));
cx_oxds = reshape(reshape(cx_cmps,Nz*Nx,length(cal.oxds))*cal.oxds/100,Nz,Nx,size(cal.oxds,2));

wt0 = (cal.perCx-cm)./(cal.perCx-cal.cphs0);
wt1 = (cal.perCm-cm)./(cal.perCm-cal.perCx);
wt2 = (cal.cphs1-cm)./(cal.cphs1-cal.perCm);
cm_cmps = reshape((wt0(:) .* cal.cmps(1,:) + (1-wt0(:)) .* cal.cmps(2,:)) .* (cm(:)<=cal.perCx) ...
                + (wt1(:) .* cal.cmps(2,:) + (1-wt1(:)) .* cal.cmps(3,:)) .* (cm(:)> cal.perCx & cm(:)<=cal.perCm) ...
                + (wt2(:) .* cal.cmps(3,:) + (1-wt2(:)) .* cal.cmps(4,:)) .* (cm(:)> cal.perCm) ...
                  ,Nz,Nx,length(cal.cmps));
cm_oxds = reshape(reshape(cm_cmps,Nz*Nx,length(cal.oxds))*cal.oxds/100,Nz,Nx,size(cal.oxds,2));

wt0 = (cal.perCx-c)./(cal.perCx-cal.cphs0);
wt1 = (cal.perCm-c)./(cal.perCm-cal.perCx);
wt2 = (cal.cphs1-c)./(cal.cphs1-cal.perCm);
c_cmps = reshape((wt0(:) .* cal.cmps(1,:) + (1-wt0(:)) .* cal.cmps(2,:)) .* (c(:)<=cal.perCx) ...
               + (wt1(:) .* cal.cmps(2,:) + (1-wt1(:)) .* cal.cmps(3,:)) .* (c(:)> cal.perCx & c(:)<=cal.perCm) ...
               + (wt2(:) .* cal.cmps(3,:) + (1-wt2(:)) .* cal.cmps(4,:)) .* (c(:)> cal.perCm) ...
                  ,Nz,Nx,length(cal.cmps));
c_oxds = reshape(reshape(c_cmps,Nz*Nx,length(cal.oxds))*cal.oxds/100,Nz,Nx,size(cal.oxds,2));

% update phase densities
rhom = rhom0 .* (1 - aT.*(T-cal.perT-273.15) - gC.*(cm-(cal.perCx+cal.perCm)/2));
rhox = rhox0 .* (1 - aT.*(T-cal.perT-273.15) - gC.*(cx-(cal.perCx+cal.perCm)/2));
rhof = rhof0 .* (1 - aT.*(T-cal.perT-273.15) + bP.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho  = 1./(m./rhom + x./rhox + f./rhof);  

chi   = max(TINY,min(1-TINY, x.*rho./rhox ));
phi   = max(TINY,min(1-TINY, f.*rho./rhof ));
mu    = max(TINY,min(1-TINY, m.*rho./rhom ));                                  

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

wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_m(1:end-1,:)+1./Csgr_m(2:end,:)); % melt segregation speed
wm(1  ,:)     = min(1,1-top).*wm(1  ,:);
wm(end,:)     = min(1,1-bot).*wm(end,:);
wm(:,[1 end]) = -sds*wm(:,[2 end-1]);

wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_x(1:end-1,:)+1./Csgr_x(2:end,:)); % melt segregation speed
wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
wx(end,:)     = min(1,1-bot).*wx(end,:);
wx(:,[1 end]) = -sds*wx(:,[2 end-1]);

wf = ((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2).*g0.*2./(1./Csgr_f(1:end-1,:)+1./Csgr_f(2:end,:)); % melt segregation speed
wf(1  ,:)     = min(1,1-top+fout).*wf(1  ,:);
wf(end,:)     = min(1,1-bot+fin ).*wf(end,:);
wf(:,[1 end]) = -sds*wf(:,[2 end-1]);


% update velocity divergence
Div_V(2:end-1,2:end-1) = ddz(W(:,2:end-1),h) ...                           % get velocity divergence
                       + ddx(U(2:end-1,:),h);
Div_V([1 end],:) = Div_V([2 end-1],:);                                     % apply boundary conditions
Div_V(:,[1 end]) = Div_V(:,[2 end-1]);

% update strain rates
exx(:,2:end-1) = diff(U,1,2)./h - Div_V(:,2:end-1)./2;                     % x-normal strain rate
exx([1 end],:) = exx([2 end-1],:);                                         % apply boundary conditions
exx(:,[1 end]) = exx(:,[2 end-1]);
ezz(2:end-1,:) = diff(W,1,1)./h - Div_V(2:end-1,:)./2;                     % z-normal strain rate
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
Div_rhoV =  + advect(rho(inz,inx).*m(inz,inx),0.*Um(inz,:),wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
            + advect(rho(inz,inx).*x(inz,inx),0.*Ux(inz,:),wx(:,inx),h,{ADVN,''   },[1,2],BCA) ...
            + advect(rho(inz,inx).*f(inz,inx),0.*Uf(inz,:),wf(:,inx),h,{ADVN,''   },[1,2],BCA) ...
            + advect(rho(inz,inx)            ,   U (inz,:),W (:,inx),h,{ADVN,'vdf'},[1,2],BCA);
% if step>0; VolSrc(inz,inx) = -((rho(inz,inx)-rhoo(inz,inx))./dt + Div_rhoV)./rho(inz,inx)*; end
if step>0; VolSrc(inz,inx) = -((rho(inz,inx)-rhoo(inz,inx))./dt + theta.*Div_rhoV + (1-theta).*Div_rhoVo)./rho(inz,inx); end

UBG    = - 1*mean(mean(VolSrc(inz,inx)))./2 .* (L/2-XXu);
WBG    = - 1*mean(mean(VolSrc(inz,inx)))./2 .* (D/2-ZZw);
end

UDtime = UDtime + toc;