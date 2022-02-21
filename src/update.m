%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase densities
rhom = rhom0 .* (1 - aTm.*(T-Tphs0) - gCm.*(cm-cphs0));
rhox = rhox0 .* (1 - aTx.*(T-Tphs0) - gCx.*(cx-cphs0));
rhof = rhof0 .* (1 - aTf.*(T-Tphs0) + bPf.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox + f./rhof);  rho([1 end],:) = rho([2 end-1],:);  rho(:,[1 end]) = rho(:,[2 end-1]);
chi   = x.*rho./rhox;
phi   = f.*rho./rhof;
mu    = m.*rho./rhom;
rhoBF = ((rho (2:end-2,2:end-1)+rho (3:end-1,2:end-1))/4 + (rhoo(2:end-2,2:end-1)+rhoo(3:end-1,2:end-1))/4 - rhoref);  % taken at mid-point in time
if Nx <= 10; rhoBF = repmat(mean(rhoBF,2),1,Nx-2); end


% update thermal properties
rhoCp = rho.*(m.*Cpm + x.*Cpx + f.*Cpf);
rhoDs = rho.*(m.*0   + x.*Dsx + f.*Dsf);
kT    = mu.*kTm + chi.*kTx + phi.*kTf;                                     % magma thermal conductivity

% update effective viscosity
etam  = etam0 .* exp(Em./(8.3145.*(T+273.15))-Em./(8.3145.*((Tphs0+Tphs1)/2+273.15))) ...
              .* Fmc.^(cm-(cphs0+cphs1)/2) .* Fmv.^(vm./0.01);             % variable melt viscosity
etaf  = etaf0.* ones(size(f));                                             % constant fluid viscosity
etax  = etax0.* ones(size(x));                                             % constant crysta viscosity

% get permission weights
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-16,permute(cat(3,chi,mu,phi),[3,1,2]));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum and volume flux and transfer coefficients
Kv =    ff .*kv.*thtv;
Cv = (1-ff)./[dx;dm;df].^2.*Kv;

% compose effective viscosity, segregation coefficients
eta = squeeze(sum(Kv,1));                                                  
eta = max(etamin,min(etamax,eta));                                         % limit viscosity range
if step>0; etact = (eta + etao)/2;                                         % effective viscosity in cell centres
else;      etact =  eta;  end
etaco  = (etact(1:end-1,1:end-1)+etact(2:end,1:end-1) ...                  % effective viscosity in cell corners
       +  etact(1:end-1,2:end  )+etact(2:end,2:end  ))./4;

Ksgr_x = max(1e-18,min(1e-6,chi./squeeze(Cv(1,:,:))));
Ksgr_m = max(1e-18,min(1e-6,mu ./squeeze(Cv(2,:,:))));
Ksgr_f = max(1e-18,min(1e-6,phi./squeeze(Cv(3,:,:))));

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
txx = etact .* exx;                                                        % x-normal stress
tzz = etact .* ezz;                                                        % z-normal stress
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

% update phase segregation speeds
if bndmode==4; sds = -1;      % no slip for 'all sides(4)'
else;          sds = +1; end  % free slip for other types

wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0.*((Ksgr_x(1:end-1,:)+Ksgr_x(2:end,:))/2); % crystal segregation speed
wx([1 end],:) = 0;
wx(:,[1 end]) = sds*wx(:,[2 end-1]);

wf = any(v(:)>1e-6).*((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0.*((Ksgr_f(1:end-1,:)+Ksgr_f(2:end,:))/2); % fluid segregation speed
wf([1 end],:) = [fout;fin].*wf([2 end-1],:);
wf(:,[1 end]) = sds*wf(:,[2 end-1]);

wm = 0.*((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0.*((Ksgr_m(1:end-1,:)+Ksgr_m(2:end,:))/2); % melt segregation speed
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
Div_rhoV =  + advection(rho.*f,0.*U,wf,h,ADVN,'flx') ...
            + advection(rho.*x,0.*U,wx,h,ADVN,'flx') ...
            + advection(rho   ,U   ,W ,h,ADVN,'flx');
% VolSrc = -((rho-rhoo)./dt + Div_rhoV - rho.*Div_V)./rho;
VolSrc = -((rho-rhoo)/dt + (Div_rhoV - rho.*Div_V + Div_rhoVo)/2)./rho;

UBG    = + mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (L/2-XXu);
WBG    = + mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (D/2-ZZw);
