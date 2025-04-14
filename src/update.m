%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************
tic;

% update phase indicators
hasx = x >= eps^0.5;
hasf = f >= eps^0.5;
hasm = m >= eps^0.5;

% update phase oxide compositions
c_oxd  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cx_oxd = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% update phase mineral end-member compositions
c_mem  = reshape(reshape(c ,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cm_mem = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);
cx_mem = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_mem,Nz,Nx,cal.nmem);

% update phase mineral systems composition for solid assemblage
c_msy  = reshape(reshape( c_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cm_msy = reshape(reshape(cm_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);
cx_msy = reshape(reshape(cx_mem,Nz*Nx,cal.nmem)*cal.msy_mem.',Nz,Nx,cal.nmsy);

% update mineral systems oxide compositions for solid assemblage
cx_msy_oxd = zeros(Nz,Nx,cal.nmsy,cal.noxd);
for j = 1:cal.nmsy
    cx_msy_oxd(:,:,j,:) = reshape(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(reshape(cx_mem(:,:,cal.msy_mem(j,:)==1),Nz*Nx,sum(cal.msy_mem(j,:)==1))+1e-32,2),Nz,Nx,1,cal.noxd);
end

cm_oxd_all = zeros(size(c,1),size(c,2),9);
cx_oxd_all = zeros(size(c,1),size(c,2),9);
 c_oxd_all = zeros(size(c,1),size(c,2),9);
if cal.noxd>9
    cm_oxd_all = cm_oxd(:,:,cal.ioxd);
    cx_oxd_all = cx_oxd(:,:,cal.ioxd);
     c_oxd_all =  c_oxd(:,:,cal.ioxd);
else
    cm_oxd_all(:,:,cal.ioxd) = cm_oxd;
    cx_oxd_all(:,:,cal.ioxd) = cx_oxd;
     c_oxd_all(:,:,cal.ioxd) = c_oxd;
end

% get trace element phase compositions
Ktrc = zeros(Nz,Nx,cal.ntrc);
trcm = zeros(Nz,Nx,cal.ntrc);
trcx = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* c_mem(:,:,j)./100; end

    trcm(:,:,i)  = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcx(:,:,i)  = trc(:,:,i)./(m./Ktrc(:,:,i) + x);
end

% update phase densities
rhom0  = reshape(DensityX(reshape(cm_oxd_all,Nz*Nx,9),Tref,Pref./1e8)    ,Nz,Nx);
rhox0  = reshape(sum(reshape(cx_mem/100,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx);
rhof0  = cal.rhof0.*ones(size(T))                                               ;

rhom   = rhom0 .* (1 - aTm.*(T-Tref) + bPm.*(Pt-Pref));
rhox   = rhox0 .* (1 - aTx.*(T-Tref) + bPx.*(Pt-Pref));
rhof   = rhof0 .* (1 - aTf.*(T-Tref) + bPf.*(Pt-Pref));

rho0   = 1./(m./rhom0 + x./rhox0 + f./rhof0);
rho    = 1./(m./rhom  + x./rhox  + f./rhof );

rhomw  = (rhom(icz(1:end-1),:)+rhom(icz(2:end),:))/2;
rhoxw  = (rhox(icz(1:end-1),:)+rhox(icz(2:end),:))/2;
rhofw  = (rhof(icz(1:end-1),:)+rhof(icz(2:end),:))/2;

rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;

rhoref = mean(rhow,2);

Drhom  = rhomw - rhow;
Drhox  = rhoxw - rhow;
Drhof  = rhofw - rhow;
Drho   = rhow  - rhoref;

rhoW   = rhow.*W(:,2:end-1);
rhoU   = rhou.*U(2:end-1,:);

% convert weight to volume fraction, update bulk density
chi    = max(0,min(1, x.*rho./rhox ));
phi    = max(0,min(1, f.*rho./rhof ));
mu     = max(0,min(1, m.*rho./rhom ));

chi_mem = reshape(reshape(cx_mem/100.*rhox,Nz*Nx,cal.nmem)./cal.rhox0,Nz,Nx,cal.nmem);
chi_mem = chi_mem./sum(chi_mem,3);

% update thermal parameters
aT = mu.*aTm + chi.*aTx + phi.*aTf;
bP = mu.*bPm + chi.*bPx + phi.*bPf;
kT = mu.*kTm + chi.*kTx + phi.*kTf;
cP = mu.*cPm + chi.*cPx + phi.*cPf;
RhoCp = mu.*rhom.*cPm + chi.*rhox.*cPx + phi.*rhof.*cPf;
Adbt  = mu.*aTm./rhom./cPm + chi.*aTx./rhox./cPx + phi.*aTf./rhof./cPf;

% update lithostatic pressure
Pti = Pt;
if Nz==1; Pt    = max(Ptop/100,Ptop.*ones(size(Pt)) + Pcouple*(Pchmb + P(2:end-1,2:end-1))); else
    Pl(1,:)     = repmat(rhoref(1).*g0.*h/2,1,Nx) + Ptop;
    Pl(2:end,:) = Pl(1,:) + repmat(cumsum(rhoref(2:end-1).*g0.*h),1,Nx);
    Pt          = max(Ptop/100,Pl + Pcouple*(Pchmb + P(2:end-1,2:end-1)));
end
Pt = alpha.*Pt + (1-alpha).*Pti;
upd_Pt = Pt-Pti;
% dPtdt  = (Pt - Pto)/dt;
% dPtdt = ((a1*Pt-a2*Pto-a3*Ptoo)/dt - (b2*dPtdto + b3*dPtdtoo))/b1;

% update effective constituent sizes
dm = dm0.*(1-mu ).^0.5;
dx = dx0.*(1-chi).^0.5;
df = df0.*(1-phi).^0.5;

% update pure phase viscosities
etam   = reshape(Giordano08(reshape(cm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);
etax0  = reshape(prod(cal.etax0(1:end-1).^reshape(chi_mem(:,:,1:end-1)+eps,Nz*Nx,cal.nmem-1),2),Nz,Nx);
etax   = etax0 .* ones(size(chi)) .* exp(cal.Eax./(8.3145.*T)-cal.Eax./(8.3145.*(Tref+273.15)));
etaf   = cal.etaf0 .* ones(size(phi));

% get coefficient contrasts
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
% kf = permute(cat(3,dx.^2./etax,dm.^2./etam,df.^2./etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
% Mf = permute(repmat(kf,1,1,1,3),[4,1,2,3])./permute(repmat(kf,1,1,1,3),[1,4,2,3]);

% get permission weights
dd = max(eps^0.5,min(1-eps^0.5,permute(cat(3,dx ,dm ,df ),[3,1,2])));
ff = max(eps^0.5,min(1-eps^0.5,permute(cat(3,chi,mu ,phi),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./cal.BB).^(1./cal.CC);  Sf = Sf./sum(Sf,2);
Xf = sum(cal.AA.*Sf,2).*FF + (1-sum(cal.AA.*Sf,2)).*Sf;

% get momentum flux and transfer coefficients
thtv = squeeze(prod(Mv.^Xf,2));
Kv   = ff.*kv.*thtv;
Cv   = Kv.*(1-ff)./dd.^2;

% get volume flux and transfer coefficients
% thtf = squeeze(prod(Mf.^Xf,2));
% Kf   = ff.*kf.*thtf;
% Cf   = Kf.*(1-ff)./dd.^2;

% get effective viscosity
eta0 = squeeze(sum(Kv,1)); if Nx==1; eta0 = eta0.'; end

% get segregation cofficients
Ksgr   = ff./Cv;

Ksgr_x = squeeze(Ksgr(1,:,:)) + eps^2; if Nx==1; Ksgr_x = Ksgr_x.'; end
Ksgr_m = squeeze(Ksgr(2,:,:)) + eps^2; if Nx==1; Ksgr_m = Ksgr_m.'; end
Ksgr_f = squeeze(Ksgr(3,:,:)) + eps^2; if Nx==1; Ksgr_f = Ksgr_f.'; end

if ~calibrt % skip the following if called from calibration script

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V/3;                                 % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V/3;                                 % z-normal strain rate
exz = (diff(U,1,1)./h+diff(W,1,2)./h)/2;                                   % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% extract potential densities
rhomp = rhom0 .* (1 - aTm.*(Tp-Tref));
rhoxp = rhox0 .* (1 - aTx.*(Tp-Tref));
rhofp = rhof0 .* (1 - aTf.*(Tp-Tref));

rhop  = 1./(m./rhomp + x./rhoxp + f./rhofp);

% update velocity magnitude
if Nx==1 && Nz==1; Vel = 0;
elseif Nx==1
    idz = (1:Nz)';  % grid indices
    half_steps = max(1,floor(Delta_cnv ./ (2 * h)));  % half mixing length in grid steps
    
    ip = idz + half_steps;  % upper indices
    im = idz - half_steps;  % lower indices
    
    ip = min(ip, Nz);  % clamp indices to valid range [1, Nz]
    im = max(im, 1 );  % clamp indices to valid range [1, Nz]
    
    drhoz   = max(0, -(rhop(ip,:)-rhop(im,:)) ) + 1e-6.*rhop; % central density contrast across mixing length
    for i=1:10
        drhoz = drhoz + diffus(drhoz,1/8*ones(size(drhoz)),1,[1,2],BCD);
    end
    Vel     = 2/9*drhoz.*g0.*(Delta_cnv/2).^2./eta;
else
    Vel = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
             + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);
end

% update diffusion parameters
if Nx==1 && Nz==1; ke = 0; fRe1 = 1; fRe100 = 1;
elseif Nx==1
    ke     = (ke + Vel.*Delta_cnv)/2;                                      % convective mixing diffusivity
    fRe1   = 1;
    fRe100 = 1;
else
    eII0   = eta0./rho./(3*h)^2;
    eIIe   = eII .* (1-exp(-eII./eII0)+eps);
    ke     = (ke + eIIe.*Delta_cnv.^2)/2;                                               % turbulent eddy diffusivity
    fRe1   = (1-exp(-Re./1  )+eps);
    fRe100 = (1-exp(-Re./100)+eps);
end
ke = 1./(1./kmax + 1./ke) + kmin;
kwm = abs(rhom-rho).*g0.*Ksgr_m.*Delta_sgr.*hasm;                          % segregation diffusivity
kwx = abs(rhox-rho).*g0.*Ksgr_x.*Delta_sgr.*hasx;                          % segregation diffusivity
kwf = abs(rhof-rho).*g0.*Ksgr_f.*Delta_sgr.*hasf;                          % segregation diffusivity
km  = (kwm+ke.*fRe1).*mu ;                                                 % regularised melt  fraction diffusion 
kx  = (kwx+ke.*fRe1).*chi;                                                 % regularised solid fraction diffusion 
kf  = (kwf+ke.*fRe1).*phi;                                                 % regularised fluid fraction diffusion 
ks  = (ke./Prt.*fRe100 + kmin).*rho.*cP./T;                                % regularised heat diffusion
kc  =  ke./Sct.*fRe100 + kmin;                                             % regularised component diffusion
eta =  ke.*rho.*fRe1 + eta0;                                               % regularised momentum diffusion

etamax = etacntr.*max(min(eta(:)),etamin);
eta    = 1./(1./etamax + 1./eta) + etamin;

etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

% update dimensionless numbers
Ra     = Vel.*Delta_cnv./((kT+ks.*T)./rho./cP);
Re     = Vel.*Delta_cnv./( eta      ./rho    );
Rum    = abs(wm(1:end-1,2:end-1)+wm(2:end,2:end-1))/2./Vel;
Rux    = abs(wx(1:end-1,2:end-1)+wx(2:end,2:end-1))/2./Vel;
Ruf    = abs(wf(1:end-1,2:end-1)+wf(2:end,2:end-1))/2./Vel;
Pr     = (eta./rho)./((kT+ks.*T)./rho./cP);
Sc     = (eta./rho)./( kc                );
deltam = sqrt(mu .*Ksgr_m.*eta./(1-chi));
deltaf = sqrt(phi.*Ksgr_f.*eta./(1-chi));

% update stresses
txx = eta   .* exx;                                                        % x-normal stress
tzz = eta   .* ezz;                                                        % z-normal stress
txz = etaco .* exz;                                                        % xz-shear stress

tII = (0.5.*(txx.^2 + tzz.^2 ...
       + 2.*(txz(1:end-1,1:end-1).^2+txz(2:end,1:end-1).^2 ...
       +     txz(1:end-1,2:end  ).^2+txz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% heat dissipation (entropy production) rate
if Nz==1 && Nx==1
    diss = 0.*T;  % no dissipation in 0-D mode (no diffusion, no shear deformation, no segregation)
else
    [grdTx ,grdTz ] = gradient(T(icz,icx),h);
    diss = kT./T.*(grdTz (2:end-1,2:end-1).^2 + grdTx (2:end-1,2:end-1).^2) ...
         + exx.*txx + ezz.*tzz ...
         + 2.*(exz(1:end-1,1:end-1)+exz(2:end,1:end-1)+exz(1:end-1,2:end)+exz(2:end,2:end))./4 ...
            .*(txz(1:end-1,1:end-1)+txz(2:end,1:end-1)+txz(1:end-1,2:end)+txz(2:end,2:end))./4 ...
         +  mu./Ksgr_m .* ((wm(1:end-1,2:end-1)+wm(2:end,2:end-1))./2).^2 ...
         + chi./Ksgr_x .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2 ...
         + phi./Ksgr_f .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2;
end

UDtime = UDtime + toc;
end