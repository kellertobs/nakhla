%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update melting model component compositions
wt0 = (cal.perCm-cm(:))./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm(:))./(cal.cphs1-cal.perCm);
cm_cmp = reshape((wt0 .* [1 0 0 0] + (1-wt0) .* [0 0 1 0]) .* (cm(:)< cal.perCm) ...
       +         (wt1 .* [0 0 1 0] + (1-wt1) .* [0 0 0 1]) .* (cm(:)>=cal.perCm),Nz,Nx,cal.ncmp);

wt0 = (cal.perCx-cx(:))./(cal.perCx-cal.cphs0);
wt1 = (cal.cphs1-cx(:))./(cal.cphs1-cal.perCx);
cx_cmp = reshape((wt0 .* [1 0 0 0] + (1-wt0) .* [0 1 0 0]) .* (cx(:)< cal.perCx) ...
       +         (wt1 .* [0 1 0 0] + (1-wt1) .* [0 0 0 1]) .* (cx(:)>=cal.perCx),Nz,Nx,cal.ncmp);

c_cmp = (m.*cm_cmp + x.*cx_cmp)./(1-f);

% update mineral end-member compositions
cm_mem = reshape(reshape(cm_cmp,Nz*Nx,cal.ncmp)*cal.cmp_mem/100,Nz,Nx,cal.nmem);
cx_mem = reshape(reshape(cx_cmp,Nz*Nx,cal.ncmp)*cal.cmp_mem/100,Nz,Nx,cal.nmem);
c_mem  = (m.*cm_mem + x.*cx_mem)./(1-f);

% update phase oxide compositions
cm_oxd = reshape(reshape(cm_mem,Nz*Nx,cal.nmem)*cal.mem_oxd/100,Nz,Nx,cal.noxd);
cx_oxd = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.mem_oxd/100,Nz,Nx,cal.noxd);
c_oxd = (m.*cm_oxd + x.*cx_oxd)./(1-f);

% update mineral systems composition for solid assemblage
cx_msy = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
cx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    cx_msy_oxd(:,:,j,:) = reshape(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
end

% update phase densities
wtm      = zeros(Nz*Nx,9);
wtm(:,1) = reshape(cm_oxd(:,:,1),Nz*Nx,1).*100; % SiO2
wtm(:,2) = reshape(cm_oxd(:,:,2),Nz*Nx,1).*100; % TiO2
wtm(:,3) = reshape(cm_oxd(:,:,3),Nz*Nx,1).*100; % Al2O3
wtm(:,4) = reshape(cm_oxd(:,:,4),Nz*Nx,1).*100; % FeO
wtm(:,5) = reshape(cm_oxd(:,:,5),Nz*Nx,1).*100; % MgO
wtm(:,6) = reshape(cm_oxd(:,:,6),Nz*Nx,1).*100; % CaO
wtm(:,7) = reshape(cm_oxd(:,:,7),Nz*Nx,1).*100; % Na2O
wtm(:,8) = reshape(cm_oxd(:,:,8),Nz*Nx,1).*100; % K2O
wtm(:,9) = reshape(vm    (:,:  ),Nz*Nx,1).*100; % H2O
rhom   = reshape(DensityX(wtm,T(:)-273.15,Pt(:)./1e8),Nz,Nx);
rhox   = reshape(sum(reshape(cx_mem,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx) .* (1 - cal.aT.*(T-cal.perT-273.15));
rhof   = cal.rhof0 .* (1 - cal.aT.*(T-cal.perT-273.15) + cal.bP.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho    = 1./(m./rhom + x./rhox + f./rhof);

rhofz  = (rho(1:end-1,:)+rho(2:end,:))/2;

chi    = max(0,min(1, x.*rho./rhox ));
phi    = max(0,min(1, f.*rho./rhof ));
mu     = max(0,min(1, m.*rho./rhom ));

% update melt viscosity
etam   = reshape(giordano08(wtm,T(:)-273.15),Nz,Nx);

% effective mixture shear viscosity (Costa et al., 2009)
hh     = (1-cal.xi).*erf(sqrt(pi)./(2.*(1-cal.xi)).*(max(TINY^0.5,chi)./cal.chi_pck).*(1+(max(TINY^0.5,chi)./cal.chi_pck).^cal.gamma));
eta    = etam .* (1+(max(TINY^0.5,chi)./cal.chi_pck).^cal.delta) .* (1-hh).^-cal.Bchi .* max(TINY^0.5,1-phi).^-cal.Bphi;

% phase segregation coefficients
Ksgr_x = 2/9*dx^2./eta/sgrreg                                                                      ;
Ksgr_f = 2/9*df^2./eta/sgrreg + dx^2/cal.bf./cal.etaf0/sgrreg.*max(TINY^0.5,phi-cal.cf).^(cal.nf-1).*max(TINY^0.5,1-phi).^cal.mf;
Ksgr_m =                            dx^2/cal.bm./    etam /sgrreg.*max(TINY^0.5,mu -cal.cm).^(cal.nm-1).*max(TINY^0.5,1-mu ).^cal.mm;

% bound and regularise viscosity
if ~calibrt; etamax = etacntr.*min(eta(:)); else; etamax = 1e+32.*min(eta(:)); end
eta    = (etamax.^-0.5 + eta.^-0.5).^-2 .* cnvreg;
etaco  = (eta([1,1:end],[1  ,1:end]).*eta([1:end,end],[1  ,1:end]) ...
       .* eta([1,1:end],[1:end,end]).*eta([1:end,end],[1:end,end])).^0.25;

if ~calibrt % skip the following if called from calibration script

% diffusion parameters
kW  = Vel/10*h/10;                                                         % convection fluctuation diffusivity
kwx = abs((rhox-rho).*g0.*Ksgr_x*dx*10);                                   % segregation fluctuation diffusivity
kwf = abs((rhof-rho).*g0.*Ksgr_f*df*10);                                   % segregation fluctuation diffusivity
kx  = chi.*mu.*(kwx + kW + mink);                                          % solid fraction diffusion 
kf  = phi.*mu.*(kwf + kW + mink);                                          % fluid fraction diffusion 
kT  = kT0 + mu.*rho.*cP.*(phi.*kwf + chi.*kwx + kW + mink);                % heat diffusion
ks  = kT./T;                                                               % entropy diffusion

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence
      

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./2;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./2;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

% update tensor magnitudes
eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + TINY;

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end).^2)/4)).^0.5 + TINY;

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1  
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx,grdTz] = gradient(T([1,1:end,end],[1,1:end,end]),h);
    diss = ks.*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2) ...
        + exx.*txx + ezz.*tzz ...
        + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
           .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
        +  mu./Ksgr_m .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
        + chi./Ksgr_x .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
        + phi./Ksgr_f .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2;
end

UDtime = UDtime + toc;
end