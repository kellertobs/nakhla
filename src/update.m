%%*****  UPDATE PARAMETERS & AUXILIARY FIELDS  ****************************

% update phase densities
rhom = rhom0 .* (1 - aT.*(T-perT-273.15) - gC.*(cm-(perCx+perCm)/2));
rhox = rhox0 .* (1 - aT.*(T-perT-273.15) - gC.*(cx-(perCx+perCm)/2));
rhof = rhof0 .* (1 - aT.*(T-perT-273.15) + bP.*(Pt-Ptop ));

% convert weight to volume fraction, update bulk density
rho  = 1./(m./rhom + x./rhox + f./rhof);  
rho([1 end],:) = rho([2 end-1],:);  
rho(:,[1 end]) = rho(:,[2 end-1]);

chi   = x.*rho./rhox;
phi   = f.*rho./rhof;
mu    = m.*rho./rhom;

% update thermal properties
Ds    = x.*Dsx + f.*Dsf;
kT    = (mu.*kTm + chi.*kTx + phi.*kTf)./T;                                % magma thermal conductivity

% update effective viscosity
etam  = etam0 .* exp(Em./(8.3145.*T)-Em./(8.3145.*(perT+273.15))) ...
              .* Fmc.^((cm-(perCx+perCm)/2)./(cphs1-cphs0)) ...
              .* Fmv.^(vm./0.01);                                          % variable melt viscosity
etaf  = etaf0.* ones(size(f));                                             % constant fluid viscosity
etax  = etax0.* ones(size(x));                                             % constant crysta viscosity

% get permission weights
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]);
 
ff = max(1e-6,min(1-1e-6,permute(cat(3,chi,mu,phi),[3,1,2])));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get effective viscosity
eta    = squeeze(sum(ff.*kv.*thtv,1));
eta    = (1./etamax + 1./eta).^-1 + etamin;
eta([1 end],:) = eta([2 end-1],:);  
eta(:,[1 end]) = eta(:,[2 end-1]);
etaco  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...
       +  eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;

% get segregation coefficients
Csgr = ((1-ff)./[dx;1e-16;df].^2.*kv.*thtv).^-1;

Csgr_x = squeeze(Csgr(1,:,:)) + 1e-18;  Csgr_x([1 end],:) = Csgr_x([2 end-1],:);  Csgr_x(:,[1 end]) = Csgr_x(:,[2 end-1]);
Csgr_f = squeeze(Csgr(3,:,:)) + 1e-18;  Csgr_f([1 end],:) = Csgr_f([2 end-1],:);  Csgr_f(:,[1 end]) = Csgr_f(:,[2 end-1]);

% update phase segregation speeds
wx = (((rhox(1:end-1,:)+rhox(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0).*2./(1./Csgr_x(1:end-1,:)+1./Csgr_x(2:end,:)); % crystal segregation speed
wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
wx(end,:)     = min(1,1-bot).*wx(end,:);
wx(:,[1 end]) = -sds*wx(:,[2 end-1]);
for i = 1:round(delta)
    wx(2:end-1,2:end-1) = wx(2:end-1,2:end-1) + diff(wx(:,2:end-1),2,1)./8 + diff(wx(2:end-1,:),2,2)./8;
    wx(1  ,:)     = min(1,1-top).*wx(1  ,:);
    wx(end,:)     = min(1,1-bot).*wx(end,:);
    wx(:,[1 end]) = -sds*wx(:,[2 end-1]);
end

wf = (((rhof(1:end-1,:)+rhof(2:end,:))/2-(rho(1:end-1,:)+rho(2:end,:))/2)*g0).*2./(1./Csgr_f(1:end-1,:)+1./Csgr_f(2:end,:)); % fluid segregation speed
wf(1  ,:) = min(1,1-top+fout).*wf(1  ,:);
wf(end,:) = min(1,1-bot+fin ).*wf(end,:);
wf(:,[1 end]) = -sds*wf(:,[2 end-1]);
for i = 1:round(delta)
    wf(2:end-1,2:end-1) = wf(2:end-1,2:end-1) + diff(wf(:,2:end-1),2,1)./8 + diff(wf(2:end-1,:),2,2)./8;
    wf(1  ,:) = min(1,1-top+fout).*wf(1  ,:);
    wf(end,:) = min(1,1-bot+fin ).*wf(end,:);
    wf(:,[1 end]) = -sds*wf(:,[2 end-1]);
end

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
txx = eta .* exx;                                                          % x-normal stress
tzz = eta .* ezz;                                                          % z-normal stress
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
     +  x(2:end-1,2:end-1).^2./Csgr_x(2:end-1,2:end-1) .* ((wx(1:end-1,2:end-1)+wx(2:end,2:end-1))./2).^2  ...
     +  f(2:end-1,2:end-1).^2./Csgr_f(2:end-1,2:end-1) .* ((wf(1:end-1,2:end-1)+wf(2:end,2:end-1))./2).^2 ...
     +  kT(2:end-1,2:end-1).*(grdTz(2:end-1,2:end-1).^2 + grdTx(2:end-1,2:end-1).^2);

% update volume source
Div_rhoV =  + advection(rho.*f,0.*U,wf,h,ADVN,'flx') ...
            + advection(rho.*x,0.*U,wx,h,ADVN,'flx') ...
            + advection(rho   ,U   ,W ,h,ADVN,'flx');
if step>0; VolSrc  = -((rho-rhoo)./dt + Div_rhoV - rho.*Div_V)./rho; end

UBG    = - 2*mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (L/2-XXu);
WBG    = - 0*mean(mean(VolSrc(2:end-1,2:end-1)))./2 .* (D/2-ZZw);