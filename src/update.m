%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update phase oxide compositions
c_oxd  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% update phase mineral end-member compositions
c_mem  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cm_mem = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cx_mem = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);

% update mineral systems composition for solid assemblage
cx_msy = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
cx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    cx_msy_oxd(:,:,j,:) = reshape(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
end

cm_oxd_all = zeros(size(c,1),size(c,2),9);
cm_oxd_all(:,:,cal.ioxd) = cm_oxd;
cx_oxd_all = zeros(size(c,1),size(c,2),9);
cx_oxd_all(:,:,cal.ioxd) = cx_oxd;
 c_oxd_all = zeros(size(c,1),size(c,2),9);
 c_oxd_all(:,:,cal.ioxd) = c_oxd;

% update phase densities
rhom   = reshape(DensityX(reshape(cm_oxd_all,Nz*Nx,9),T0,Ptop./1e8)      ,Nz,Nx) .* (1 - aT.*(T-T0-273.15) + bPm.*(Pt-Ptop));
rhox   = reshape(sum(reshape(cx_mem/100,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx) .* (1 - aT.*(T-T0-273.15) + bPx.*(Pt-Ptop));
rhof   = cal.rhof0                                                               .* (1 - aT.*(T-T0-273.15) + bPf.*(Pt-Ptop));

% convert weight to volume fraction, update bulk density
rho    = 1./(m./rhom + x./rhox + f./rhof);

chi    = max(0,min(1, x.*rho./rhox ));
phi    = max(0,min(1, f.*rho./rhof ));
mu     = max(0,min(1, m.*rho./rhom ));

rhofx  = (rho(:,[end,1:end])+rho(:,[1:end,1]))/2;
rhofz  = (rho([1,1:end],:)+rho([1:end,end],:))/2;
% update lithostatic pressure
if Nz==1; Pt = Ptop.*ones(size(Tp)); else
    Pt(1,:)     = repmat(mean(rhofz(1,:),2).*g0.*h/2,1,Nx) + Ptop;
    Pt(2:end,:) = Pt(1,:) + repmat(cumsum(mean(rhofz(2:end-1,:),2).*g0.*h),1,Nx);
end

% update melt viscosity
etam   = reshape(Giordano08(reshape(cm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);

% effective mixture shear viscosity (Costa et al., 2009)
hh     = (1-cal.xi).*erf(sqrt(pi)./(2.*(1-cal.xi)).*(max(TINY^0.5,chi)./cal.chi_pck).*(1+(max(TINY^0.5,chi)./cal.chi_pck).^cal.gamma));
eta    = etam .* (1+(max(TINY^0.5,chi)./cal.chi_pck).^cal.delta) .* (1-hh).^-cal.Bchi .* max(TINY^0.5,1-phi).^-cal.Bphi;

% phase segregation coefficients
Ksgr_x = 2/9*dx^2./eta                                                                      ;
Ksgr_f = 2/9*df^2./eta + dx^2/cal.bf./cal.etaf0.*max(TINY^0.5,phi-cal.cf).^(cal.nf-1).*max(TINY^0.5,1-phi).^cal.mf;
Ksgr_m =                 dx^2/cal.bm./    etam .*max(TINY^0.5,mu -cal.cm).^(cal.nm-1).*max(TINY^0.5,1-mu ).^cal.mm;

if ~calibrt % skip the following if called from calibration script

% update velocity magnitude
Vel = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
         + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./2;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./2;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + TINY;

% update diffusion parameters
W0  = (Vel./mean(Vel(:)+TINY))./4.*mean(abs(rho-mean(rho,2)).*g0.*(D/10)^2./eta,'all');
wx0 = abs(wx(1:end-1,2:end-1)+wx(2:end,2:end-1))/2;
wf0 = abs(wf(1:end-1,2:end-1)+wf(2:end,2:end-1))/2;
Ra0 = W0.*D/10./(kT0./rho./cP);
Re0 = W0.*rho.*D/10./eta;

if Nx==1 && Nz==1; kW = 0;
else              
kW = (kW + 2.*eII.*(0.18*Delta).^2 .* (1-min(1,topshape+botshape+sdsshape)*0.9))/2;
end
kwx = wx0*dx*10;                                                           % segregation fluctuation diffusivity
kwf = wf0*df*10;                                                           % segregation fluctuation diffusivity
kx  = chi.*(kwx + kW/Prt);                                                 % solid fraction diffusion 
kf  = phi.*(kwf + kW/Prt);                                                 % fluid fraction diffusion 
ks  = rho.*cP./T.*(phi.*kwf + chi.*kwx + kW/Prt);                          % regularised heat diffusion
kc  = rho.*(phi.*kwf + chi.*kwx + kW/Prt);                                 % regularised component diffusion
eta = eta + rho.*(phi.*kwf + chi.*kwx + kW);

etamax = etacntr.*max(min(eta(:)),etamin);
eta    = (etamax.^-0.5 + eta.^-0.5).^-2 + etamin;

etaco  = (eta([1,1:end],[end,1:end]).*eta([1:end,end],[end,1:end]) ...
       .* eta([1,1:end],[1:end,1  ]).*eta([1:end,end],[1:end,1  ])).^0.25;

Ra  = Vel.*D/10./((kT0+ks.*T)./rho./cP);
Re  = Vel.*rho.*D/10./eta;

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end).^2)/4)).^0.5 + TINY;

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx ,grdTz ] = gradient(T ([1,1:end,end],[1,1:end,end]),h);
    [grdTpx,grdTpz] = gradient(Tp([1,1:end,end],[1,1:end,end]),h);
    diss = kT0./T.*(grdTz (2:end-1,2:end-1).^2 + grdTx (2:end-1,2:end-1).^2) ...
         + ks    .*(grdTpz(2:end-1,2:end-1).^2 + grdTpx(2:end-1,2:end-1).^2) ...
         + exx.*txx + ezz.*tzz ...
         + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
            .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
         +  mu./Ksgr_m .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
         + chi./Ksgr_x .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
         + phi./Ksgr_f .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2;
end

UDtime = UDtime + toc;
end