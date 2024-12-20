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
rhom0  = reshape(DensityX(reshape(cm_oxd_all,Nz*Nx,9),Tref,Pref./1e8)    ,Nz,Nx);
rhox0  = reshape(sum(reshape(cx_mem/100,Nz*Nx,cal.nmem)./cal.rhox0,2).^-1,Nz,Nx);
rhof0  = cal.rhof0                                                              ;

rhom   = rhom0 .* (1 - aTm.*(T-Tref) + bPm.*(Pt-Pref));
rhox   = rhox0 .* (1 - aTx.*(T-Tref) + bPx.*(Pt-Pref));
rhof   = rhof0 .* (1 - aTf.*(T-Tref) + bPf.*(Pt-Pref));

rho    = 1./(m./rhom + x./rhox + f./rhof);

rhofz  = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhofx  = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;

rhoW = rhofz.*W(:,2:end-1);
rhoU = rhofx.*U(2:end-1,:);

% convert weight to volume fraction, update bulk density
chi    = max(0,min(1, x.*rho./rhox ));
phi    = max(0,min(1, f.*rho./rhof ));
mu     = max(0,min(1, m.*rho./rhom ));

phix_mem = reshape(reshape(cx_mem/100.*rhox,Nz*Nx,cal.nmem)./cal.rhox0,Nz,Nx,cal.nmem);
phix_mem = phix_mem./sum(phix_mem,3);

% update thermal parameters
aT = mu.*aTm + chi.*aTx + phi.*aTf;
kT = mu.*kTm + chi.*kTx + phi.*kTf;
cP = mu.*cPm + chi.*cPx + phi.*cPf;
RhoCp = mu.*rhom.*cPm + chi.*rhox.*cPx + phi.*rhof.*cPf;

% update lithostatic pressure
Pti = Pt;
if Nz==1; Pt    = max(1e7,(1-alpha).*Pt + alpha.*(Ptop.*ones(size(Tp)) + Pcouple*(Pchmb + P(2:end-1,2:end-1)))); else
    Pl(1,:)     = repmat(mean(rhofz(1,:),2).*g0.*h/2,1,Nx) + Ptop;
    Pl(2:end,:) = Pl(1,:) + repmat(cumsum(mean(rhofz(2:end-1,:),2).*g0.*h),1,Nx);
    Pt          = max(1e7,(1-1).*Pt + 1.*(Pl + Pcouple*(Pchmb + P(2:end-1,2:end-1))));
end
upd_Pt = Pt-Pti;

% update effective constituent sizes
dm = dm0.*(1-mu ).^0.5;
dx = dx0.*(1-chi).^0.5;
df = df0.*(1-phi).^0.5;

% update pure phase viscosities
etam   = reshape(Giordano08(reshape(cm_oxd_all,Nz*Nx,9),T(:)-273.15),Nz,Nx);
etax0  = reshape(prod(cal.etax0(1:end-1).^reshape(phix_mem(:,:,1:end-1)+eps,Nz*Nx,cal.nmem-1),2),Nz,Nx);
etax   = etax0 .* ones(size(chi)) .* exp(cal.Eax./(8.3145.*T)-cal.Eax./(8.3145.*(Tref+273.15)));
etaf   = cal.etaf0 .* ones(size(phi));

% get coefficient contrasts
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
% kf = permute(cat(3,dx.^2./etax,dm.^2./etam,df.^2./etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
% Mf = permute(repmat(kf,1,1,1,3),[4,1,2,3])./permute(repmat(kf,1,1,1,3),[1,4,2,3]);

% get permission weights
dd = max(1e-6,min(1-1e-6,permute(cat(3,dx ,dm ,df ),[3,1,2])));
ff = max(1e-6,min(1-1e-6,permute(cat(3,chi,mu ,phi),[3,1,2])));
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

% extract non-P-dependent density
rhom_nP = rhom0 .* (1 - aTm.*(Tp-Tref));
rhox_nP = rhox0 .* (1 - aTx.*(Tp-Tref));
rhof_nP = rhof0 .* (1 - aTf.*(Tp-Tref));

rho_nP  = 1./(m./rhom_nP + x./rhox_nP + f./rhof_nP);

% detect convection layers
drhoz    = gradient(mean(rho_nP,2));
[~,zpks] = findpeaks(drhoz,'MinPeakHeight',10,'MinPeakProminence',1);
nlay = length(zpks)+1;
zlay = zeros(1,nlay+1);

zt = max(ZZ(eta>=etamax/10 & ZZ<D/2));
if isempty(zt); zlay(1) = 0; else; zlay(1) = zt; end
for iz = 1:nlay-1
    if zpks(iz)>5 && zpks(iz)<Nz-5
        zlay(iz+1) = zpks(iz)*h+h/2;
    else
        zlay(iz+1) = [];
        nlay = nlay-1;
    end
end
zb = min(ZZ(eta>=etamax/10 & ZZ>zlay(end-1) ));
if isempty(zb); zlay(end) = D; else; zlay(end) = zb; end

% limit correlation length for convective mixing to distance from layer
% and domain boundaries
Delta_cnv =zeros(Nz,1);
for iz = 1:nlay
    Delta_cnv = Delta_cnv + max(0,min(Delta_cnv0,min(ZZ-zlay(iz),zlay(iz+1)-ZZ)));
end
ind0 = Delta_cnv==0;
for i=1:10
    Delta_cnv = Delta_cnv + diffus(Delta_cnv,1/8*ones(size(Delta_cnv)),1,[1,2],BCD);
    % Delta_cnv(ind0) = 0;
    Delta_cnv([1 end]) = [h/2,h/2];
end

% update velocity magnitude
if Nx==1 && Nz==1; Vel = 0;
elseif Nx==1
    ip      = min(Nz,round((1:Nz).'+(Delta_cnv./h-1/2)/2));
    im      = max( 1,round((1:Nz).'-(Delta_cnv./h-1/2)/2));
    drhoz   = max(0, -(rho_nP(ip,:)-rho_nP(im,:)) ) + 1e-6.*rho_nP;
    for i=1:10
        drhoz = drhoz + diffus(drhoz,1/8*ones(size(drhoz)),1,[1,2],BCD);
    end
    Vel     = 2/9*drhoz.*g0.*(Delta_cnv/2).^2./eta;
else
    Vel = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
             + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);
end

% update velocity divergence
Div_V = ddz(W(:,2:end-1),h) + ddx(U(2:end-1,:),h);                         % get velocity divergence

% update strain rates
exx = diff(U(2:end-1,:),1,2)./h - Div_V./3;                                % x-normal strain rate
ezz = diff(W(:,2:end-1),1,1)./h - Div_V./3;                                % z-normal strain rate
exz = 1/2.*(diff(U,1,1)./h+diff(W,1,2)./h);                                % shear strain rate

eII = (0.5.*(exx.^2 + ezz.^2 ...
       + 2.*(exz(1:end-1,1:end-1).^2+exz(2:end,1:end-1).^2 ...
       +     exz(1:end-1,2:end  ).^2+exz(2:end,2:end  ).^2)/4)).^0.5 + eps;

% update diffusion parameters
if Nx==1 && Nz==1; kW = 0; Pr = Prt; Sc = Sct;
elseif Nx==1
    kW  = Vel.*Delta_cnv;                                                  % convective mixing diffusivity
    Pr  = Prt;
    Sc  = Sct;
else
    kW  = eII.*Delta_cnv.^2;                                               % turbulent eddy diffusivity
    Pr  = Prt ./ (1-exp(-Re./10)+eps);
    Sc  = Sct ./ (1-exp(-Re./10)+eps);
end
kwm = abs(rhom-rho).*g0.*Ksgr_m.*Delta_sgr + kmin;                         % segregation diffusivity
kwx = abs(rhox-rho).*g0.*Ksgr_x.*Delta_sgr + kmin;                         % segregation diffusivity
kwf = abs(rhof-rho).*g0.*Ksgr_f.*Delta_sgr + kmin;                         % segregation diffusivity
km  = (kwm+kW).*mu ;                                                       % regularised melt  fraction diffusion 
kx  = (kwx+kW).*chi;                                                       % regularised solid fraction diffusion 
kf  = (kwf+kW).*phi;                                                       % regularised fluid fraction diffusion 
ks  = (kW./Pr + kmin).*rho.*cP./T;                                         % regularised heat diffusion
kc  =  kW./Sc + kmin;                                                      % regularised component diffusion
eta = (kW.*rho + eta0)/2 + eta/2;                                          % regularised momentum diffusion

etamax = etacntr.*max(min(eta(:)),etamin);
eta    = 1./(1./etamax + 1./eta) + etamin;

etaco  = (eta(icz(1:end-1),icx(1:end-1)).*eta(icz(2:end),icx(1:end-1)) ...
       .* eta(icz(1:end-1),icx(2:end  )).*eta(icz(2:end),icx(2:end  ))).^0.25;

% update dimensionless numbers
Ra     = Vel.*D/10./((kT+ks.*T)./rho./cP);
Re     = Vel.*D/10./( eta       ./rho    );
Rum    = abs(wm(1:end-1,2:end-1)+wm(2:end,2:end-1))/2./Vel;
Rux    = abs(wx(1:end-1,2:end-1)+wx(2:end,2:end-1))/2./Vel;
Ruf    = abs(wf(1:end-1,2:end-1)+wf(2:end,2:end-1))/2./Vel;
Pr     = (eta./rho)./((kT+ks.*T)./rho./cP);
Sc     = (eta./rho)./( kc                 );
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