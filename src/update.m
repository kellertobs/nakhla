%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase densities
rhom = rhom0 .* (1 - aTm.*(T-perT) - gCm.*(cm-(perCx+perCm)/2));
rhox = rhox0 .* (1 - aTx.*(T-perT) - gCx.*(cx-(perCx+perCm)/2));
rhof = rhof0 .* (1 - aTf.*(T-perT) + bPf.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho   = 1./(m./rhom + x./rhox + f./rhof);  rho([1 end],:) = rho([2 end-1],:);  rho(:,[1 end]) = rho(:,[2 end-1]);

chi   = x.*rho./rhox;
phi   = f.*rho./rhof;
mu    = m.*rho./rhom;
rhoBF = (rho (2:end-2,2:end-1)+rho (3:end-1,2:end-1))/2 - rhoref;
if Nx <= 10; rhoBF = repmat(mean(rhoBF,2),1,Nx-2); end

% update thermal properties
rhoCp = rho.*(m.*Cpm + x.*Cpx + f.*Cpf);
rhoDs = rho.*(m.*0   + x.*Dsx + f.*Dsf);
kT    = mu.*kTm + chi.*kTx + phi.*kTf;                                     % magma thermal conductivity

% update effective viscosity
etam  = etam0 .* exp(Em./(8.3145.*(T+273.15))-Em./(8.3145.*(perT+273.15))) ...
              .* Fmc.^((cm-(perCx+perCm)/2)./(cphs1-cphs0)) .* Fmv.^(vm./0.01);  % variable melt viscosity
etaf  = etaf0.* ones(size(f));                                             % constant fluid viscosity
etax  = etax0.* ones(size(x));                                             % constant crysta viscosity

% get permission weights
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-4,min(1-1e-4,permute(cat(3,chi,mu,phi),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get effective viscosity
eta    = squeeze(sum(ff.*kv.*thtv,1));                                                  
eta    = (1./etamax + 1./eta).^-1 + etamin;                                % limit viscosity range
eta    = log10(eta);
for i = 1:round(delta/5)
    eta(2:end-1,2:end-1) = eta(2:end-1,2:end-1) + diff(eta(:,2:end-1),2,1)./8 + diff(eta(2:end-1,:),2,2)./8;
    eta([1 end],:) = eta([2 end-1],:);
    eta(:,[1 end]) = eta(:,[2 end-1]);
end
etact  = 10.^eta;

etaco  = (etact(1:end-1,1:end-1)+etact(2:end,1:end-1) ...                  % effective viscosity in cell corners
       +  etact(1:end-1,2:end  )+etact(2:end,2:end  ))./4;

% get segregation coefficients
Csgr = ((1-ff)./[dx;dm;df].^2.*kv.*thtv).^-1;

Csgr_x = squeeze(Csgr(1,:,:)) + 1e-16;
Csgr_m = squeeze(Csgr(2,:,:)) + 1e-16;
Csgr_f = squeeze(Csgr(3,:,:)) + 1e-16;

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
wx = ((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0.*((Csgr_x(1:end-1,:)+Csgr_x(2:end,:))/2); % crystal segregation speed
for i = 1:round(delta)
    wx(2:end-1,2:end-1) = wx(2:end-1,2:end-1) + diff(wx(:,2:end-1),2,1)./8 + diff(wx(2:end-1,:),2,2)./8;
    wx([1 end],:) = 0;
    wx(:,[1 end]) = -sds*wx(:,[2 end-1]);
end

wf = any(v(:)>1e-6).*((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0.*((Csgr_f(1:end-1,:)+Csgr_f(2:end,:))/2); % fluid segregation speed
for i = 1:round(delta)
    wf(2:end-1,2:end-1) = wf(2:end-1,2:end-1) + diff(wf(:,2:end-1),2,1)./8 + diff(wf(2:end-1,:),2,2)./8;
    wf([1 end],:) = [fout;fin].*wf([2 end-1],:);
    wf(:,[1 end]) = -sds*wf(:,[2 end-1]);
end

wm = ((rhom(1:end-1,:)+rhom(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0.*((Csgr_m(1:end-1,:)+Csgr_m(2:end,:))/2); % melt segregation speed
for i = 1:round(delta)
    wm(2:end-1,2:end-1) = wm(2:end-1,2:end-1) + diff(wm(:,2:end-1),2,1)./8 + diff(wm(2:end-1,:),2,2)./8;
    wm([1 end],:) = 0;
    wm(:,[1 end]) = -sds*wm(:,[2 end-1]);
end

% update volume source
Div_rhoV =  + advection(rho.*f,0.*U,wf,h,ADVN,'flx') ...
            + advection(rho.*x,0.*U,wx,h,ADVN,'flx') ...
            + advection(rho.*m,0.*U,wm,h,ADVN,'flx') ...
            + advection(rho   ,U   ,W ,h,ADVN,'flx');
% VolSrc  = ALPHA.*VolSrc + (1-ALPHA).* (-((rho-rhoo)./dt + THETA*(Div_rhoV - rho.*Div_V) + (1-THETA)*Div_rhoVo)./(THETA*rho));
VolSrc  = ALPHA.*VolSrc + (1-ALPHA).* (-((rho-rhoo)./dt + Div_rhoV - rho.*Div_V)./rho);

UBG    = - mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (L/2-XXu);
WBG    = - mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (D/2-ZZw);
