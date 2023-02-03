%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update mineral-like component compositions
wt0 = (cal.perCm-cm)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm)./(cal.cphs1-cal.perCm);
cm_cmp1 = reshape((wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(3,:)),Nz,Nx,cal.nc);
cm_cmp2 = reshape((wt1(:) .* cal.cmp(3,:) + (1-wt1(:)) .* cal.cmp(4,:)),Nz,Nx,cal.nc);

wt0 = (cal.perCm-[cal.cphs0,cal.cphs1])./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-[cal.cphs0,cal.cphs1])./(cal.cphs1-cal.perCm);
mincmp = zeros(1,1,cal.nc);
maxcmp = zeros(1,1,cal.nc);
mincmp(1,:,:) = min(min((wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(3,:)),(wt1(:) .* cal.cmp(3,:) + (1-wt1(:)) .* cal.cmp(4,:))));
maxcmp(1,:,:) = max(max((wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(3,:)),(wt1(:) .* cal.cmp(3,:) + (1-wt1(:)) .* cal.cmp(4,:))));
cm_cmp1 = (cm_cmp1-mincmp)./(maxcmp-mincmp);
cm_cmp2 = (cm_cmp2-mincmp)./(maxcmp-mincmp);

a = 250;
b = 0.05;
ind1 = repmat(cm,1,1,cal.nc)<=cal.perCm-b;
ind2 = repmat(cm,1,1,cal.nc)> cal.perCm+b;
ind3 = repmat(cm,1,1,cal.nc)< cal.perCm+b & cm>=cal.perCm-b;

cm_cmp = zeros(size(cm_cmp1));
cm_cmp(ind1)  =  cm_cmp1(ind1);
cm_cmp(ind2)  =  cm_cmp2(ind2);
scl = ((cal.cmp(4,:)-cal.cmp(3,:))./(cal.cphs1-cal.perCm)./100 - (cal.cmp(3,:)-cal.cmp(1,:))./(cal.perCm-cal.cphs0)./100);
scl = reshape(sign(scl).*(cal.cmp(2,:)-squeeze(mincmp).')./(squeeze(maxcmp).'-squeeze(mincmp).').*ones(Nz*Nx,cal.nc),Nz,Nx,cal.nc);
cm_cmp(ind3)  = (cm_cmp1(ind3).^(scl(ind3)*a)+cm_cmp2(ind3).^(scl(ind3)*a)).^(1./(scl(ind3)*a));
cm_cmp = mincmp + cm_cmp.*(maxcmp-mincmp);

wt0 = (cal.perCx-cx)./(cal.perCx-cal.cphs0);
wt1 = (cal.cphs1-cx)./(cal.cphs1-cal.perCx);
cx_cmp1 = reshape((wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(2,:)),Nz,Nx,cal.nc);
cx_cmp2 = reshape((wt1(:) .* cal.cmp(2,:) + (1-wt1(:)) .* cal.cmp(4,:)),Nz,Nx,cal.nc);

wt0 = (cal.perCx-[cal.cphs0,cal.cphs1])./(cal.perCx-cal.cphs0);
wt1 = (cal.cphs1-[cal.cphs0,cal.cphs1])./(cal.cphs1-cal.perCx);
mincmp = zeros(1,1,cal.nc);
maxcmp = zeros(1,1,cal.nc);
mincmp(1,:,:) = min(min((wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(2,:)),(wt1(:) .* cal.cmp(2,:) + (1-wt1(:)) .* cal.cmp(4,:))));
maxcmp(1,:,:) = max(max((wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(2,:)),(wt1(:) .* cal.cmp(2,:) + (1-wt1(:)) .* cal.cmp(4,:))));
cx_cmp1 = (cx_cmp1-mincmp)./(maxcmp-mincmp);
cx_cmp2 = (cx_cmp2-mincmp)./(maxcmp-mincmp);

a = 250;
b = 0.05;
ind1 = repmat(cx,1,1,cal.nc)<=cal.perCx-b;
ind2 = repmat(cx,1,1,cal.nc)> cal.perCx+b;
ind3 = repmat(cx,1,1,cal.nc)< cal.perCx+b & cx>=cal.perCx-b;

cx_cmp = zeros(size(cx_cmp1));
cx_cmp(ind1)  =  cx_cmp1(ind1);
cx_cmp(ind2)  =  cx_cmp2(ind2);
scl = ((cal.cmp(4,:)-cal.cmp(2,:))./(cal.cphs1-cal.perCx)./100 - (cal.cmp(2,:)-cal.cmp(1,:))./(cal.perCx-cal.cphs0)./100);
scl = reshape(sign(scl).*(cal.cmp(2,:)-squeeze(mincmp).')./(squeeze(maxcmp).'-squeeze(mincmp).').*ones(Nz*Nx,cal.nc),Nz,Nx,cal.nc);
cx_cmp(ind3)  = (cx_cmp1(ind3).^(scl(ind3)*a)+cx_cmp2(ind3).^(scl(ind3)*a)).^(1./(scl(ind3)*a));
cx_cmp = mincmp + cx_cmp.*(maxcmp-mincmp);

c_cmp = (m.*cm_cmp + x.*cx_cmp)./(1-f);

% update oxide compositions
cm_oxd = reshape(reshape(cm_cmp,Nz*Nx,cal.nc)*cal.oxd/100,Nz,Nx,cal.nc);
cx_oxd = reshape(reshape(cx_cmp,Nz*Nx,cal.nc)*cal.oxd/100,Nz,Nx,cal.nc);
c_oxd = (m.*cm_oxd + x.*cx_oxd)./(1-f);

% update phase densities
rhom = squeeze(sum(permute(cm_cmp/100,[3,1,2])./cal.rhom0.')).^-1 .* (1 - cal.aT.*(T-cal.perT-273.15) - cal.gH.*vm); if size(rhom,2)~=size(T,2); rhom = rhom(1,:).'; end
rhox = squeeze(sum(permute(cx_cmp/100,[3,1,2])./cal.rhox0.')).^-1 .* (1 - cal.aT.*(T-cal.perT-273.15)             ); if size(rhox,2)~=size(T,2); rhox = rhox(1,:).'; end
rhof = cal.rhof0 .* (1 - cal.aT.*(T-cal.perT-273.15) + cal.bP.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho    = 1./(m./rhom + x./rhox + f./rhof);

rhofz = (rho(1:end-1,:)+rho(2:end,:))/2;
rhofx = (rho(:,1:end-1)+rho(:,2:end))/2;

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

% effective mixture shear viscosity (Costa et al., 2009)
hh     = (1-cal.xi).*erf(sqrt(pi)./(2.*(1-cal.xi)).*(chi./cal.chi_pck).*(1+(chi./cal.chi_pck).^cal.gamma));
eta    = etam .* (1+(chi./cal.chi_pck).^cal.delta) .* (1-hh).^-cal.Bchi .* (1-phi).^-cal.Bphi;

% phase segregation coefficients
Ksgr_x = 2/9*cal.dx^2./eta/sgrreg                                                                       + TINY.^2;
Ksgr_f = 2/9*cal.df^2./eta/sgrreg + cal.dx^2/cal.bf./cal.etaf0/sgrreg.*phi.^(cal.nf-1).*(1-phi).^cal.mf + TINY.^2;
Ksgr_m =                            cal.dx^2/cal.bm./    etam /sgrreg.*mu .^(cal.nm-1).*(1-mu ).^cal.mm + TINY.^2;

% bound and regularise viscosity
if ~calibrt; etamax = etacntr.*min(eta(:)); else; etamax = 1e+32.*min(eta(:)); end
eta    = (etamax.^-0.5 + (eta*cnvreg).^-0.5).^-2;
etaco  = (eta(1:end-1,1:end-1).*eta(2:end,1:end-1) ...
       .* eta(1:end-1,2:end  ).*eta(2:end,2:end  )).^0.25;

if ~calibrt % skip the following if called from calibration script

wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-mean(rhofz,2)).*g0.*(Ksgr_m(1:end-1,:).*Ksgr_m(2:end,:)).^0.5; % melt segregation speed
wm(1  ,:)     = min(1,1-top).*wm(1  ,:);
wm(end,:)     = min(1,1-bot).*wm(end,:);
wm(:,[1 end]) = -sds*wm(:,[2 end-1]);

wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-mean(rhofz,2)).*g0.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5; % solid segregation speed
wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
wx(end,:)     = min(1,1-bot).*wx(end,:);
wx(:,[1 end]) = -sds*wx(:,[2 end-1]);

wf = ((rhof(1:end-1,:)+rhof(2:end,:))/2-mean(rhofz,2)).*g0.*(Ksgr_f(1:end-1,:).*Ksgr_f(2:end,:)).^0.5; % fluid segregation speed
wf(1  ,:)     = min(1,1-top+fout).*wf(1  ,:);
wf(end,:)     = min(1,1-bot+fin ).*wf(end,:);
wf(:,[1 end]) = -sds*wf(:,[2 end-1]);

% diffusion parameters
kW  = Vel/10*h/10;                                                         % convection fluctuation diffusivity
kwx = abs((rhox-rho).*g0.*Ksgr_x*dx*10);                                   % segregation fluctuation diffusivity
kwf = abs((rhof-rho).*g0.*Ksgr_f*df*10);                                   % segregation fluctuation diffusivity
kx  = chi.*(kwx + kW + mink);                                              % solid fraction diffusion 
kf  = phi.*(kwf + kW + mink);                                              % fluid fraction diffusion 
kT  = kT0 + (phi.*kwf + chi.*kwx + kW + mink).*rho.*cP;                    % heat diffusion
ks  = kT./T;                                                               % entropy diffusion

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
if Nz==3 && Nx==3  
    diss = 0.*T(inz,inx);  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx,grdTz] = gradient(T,h);
    diss = ks(inz,inx).*(grdTz(inz,inx).^2 + grdTx(inz,inx).^2) ...
        + exx(inz,inx).*txx(inz,inx) + ezz(inz,inx).*tzz(inz,inx) ...
        + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
           .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
        +  mu(inz,inx)./Ksgr_m(inz,inx) .* ((wm(inz,inx)+wm(inz,inx))./2).^2 ...
        + chi(inz,inx)./Ksgr_x(inz,inx) .* ((wx(inz,inx)+wx(inz,inx))./2).^2 ...
        + phi(inz,inx)./Ksgr_f(inz,inx) .* ((wf(inz,inx)+wf(inz,inx))./2).^2;
end

% update volume source
if step>0 && ~restart
    drhodt  = - advect(M(inz,inx),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...  % melt  advection
              - advect(X(inz,inx),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA) ...  % xtal  advection
              - advect(F(inz,inx),Uf(inz,:),Wf(:,inx),h,{ADVN,''},[1,2],BCA);     % fluid advection
    res_rho = (a1*(rho(inz,inx)-rhoo(inz,inx))/dt + a2*(rho(inz,inx)-rhooo(inz,inx))/(dt+dto)) - (b1*drhodt + b2*drhodto);
    % res_DivV   = (alpha1*rho(inz,inx) - alpha2*rhoo(inz,inx) - alpha3*rhooo(inz,inx))./dt + (beta1*Div_rhoV + beta2*Div_rhoVo + beta3*Div_rhoVoo);  % get residual of mixture mass conservation
    % res_DivV   = (alpha1*rho(inz,inx) - alpha2*rhoo(inz,inx) - alpha3*rhooo(inz,inx))./dt + Div_rhoV;  % get residual of mixture mass conservation
    % res_DivV   = (rho(inz,inx) - rhoo(inz,inx))./dt + Div_rhoV;  % get residual of mixture mass conservation
    VolSrc   = Div_V(inz,inx) - lambda*res_rho./rho(inz,inx);  % correct volume source term by scaled residual

    UBG    = - mean(VolSrc,'all')./2 .* (L/2-XXu);
    WBG    = - mean(VolSrc,'all')./2 .* (D/2-ZZw);

end

UDtime = UDtime + toc;
end