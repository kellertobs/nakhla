% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 && save_op == 1
    parfile = [opdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN NAKHLA MODEL | %s  *************\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);

load ocean;                  % load custom colormap
run(['../cal/cal_',calID]);  % load melt model calibration
if periodic % periodic sides
    BCA     =  {'','periodic'};  % boundary condition on advection (top/bot, sides)
    BCD     =  {'','periodic'};  % boundary condition on advection (top/bot, sides)
else % closed sides
    BCA     =  {'',''};  % boundary condition on advection (top/bot, sides)
    BCD     =  {'',''};  % boundary condition on advection (top/bot, sides) 
end

Dsx = -cal.dS;
Dsf = cal.dS+100;

% normalise major components to anhydrous unit sum, rescale to hydrous
c0(1:end-1) = c0(1:end-1)./sum(c0(1:end-1)).*(1-c0(end));
c1(1:end-1) = c1(1:end-1)./sum(c1(1:end-1)).*(1-c1(end));
cwall(:,1:end-1) = cwall(:,1:end-1)./sum(cwall(:,1:end-1),2).*(1-cwall(:,end));
dcg   = dcg-mean(dcg);
dcr   = dcr-mean(dcr);

% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
Xf        = (Xc(1:end-1)+Xc(2:end))./2;
Zf        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu] = meshgrid(Xf,Zc);
[XXw,ZZw] = meshgrid(Xc,Zf);
Xc        = Xc(2:end-1);
Zc        = Zc(2:end-1);
[XX,ZZ]   = meshgrid(Xc,Zc);

Nx = length(Xc);
Nz = length(Zc);

% get smoothed initialisation field
rng(seed);
smth = smth*Nx*Nz*1e-4;
rp   = randn(Nz,Nx);
for i = 1:round(smth)
    rp = rp + diffus(rp,1/8*ones(size(rp)),1,[1,2],BCD);
    rp = rp - mean(mean(rp));
end
rp = rp./max(abs(rp(:)));

gp = exp(-(XX-D/2).^2/(D/8)^2 - (ZZ-D/2).^2/(D/8)^2);

% get mapping arrays
NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set up shape functions for initial boundary layers
topinit = zeros(size(ZZ));
botinit = zeros(size(ZZ));
sdsinit = zeros(size(XX));
if any(bnd_h)
    switch bndmode
        case 0  % none
        case 1  % top only
            topinit = (1+erf( ( -ZZ+bnd_h(1))/bnd_w))/2;
        case 2  % bot only
            botinit = (1+erf(-(D-ZZ-bnd_h(2))/bnd_w))/2;
        case 3  % top/bot only
            topinit = (1+erf( ( -ZZ+bnd_h(1))/bnd_w))/2;
            botinit = (1+erf(-(D-ZZ-bnd_h(2))/bnd_w))/2;
        case 4 % all walls
            topinit = (1+erf( ( -ZZ+bnd_h(1))/bnd_w))/2;
            botinit = (1+erf(-(D-ZZ-bnd_h(2))/bnd_w))/2;
            sdsinit = (1+erf( ( -XX+bnd_h(3))/bnd_w))/2 ...
                    + (1+erf(-(L-XX-bnd_h(3))/bnd_w))/2;
        case 5 % only walls
            sdsinit = (1+erf( ( -XX+bnd_h(3))/bnd_w))/2 ...
                    + (1+erf(-(L-XX-bnd_h(3))/bnd_w))/2;
    end
    sdsinit = max(0,sdsinit-topinit-botinit);
end

% set up shape functions for transient boundary layers
topshape = zeros(size(ZZ));
botshape = zeros(size(ZZ));
sdsshape = zeros(size(XX));
if ~any(bnd_h)
    switch bndmode
        case 0  % none
        case 1  % top only
            topshape = exp( ( -ZZ+h/2)/bnd_w);
        case 2  % bot only
            botshape = exp(-(D-ZZ-h/2)/bnd_w);
        case 3  % top/bot only
            topshape = exp( ( -ZZ+h/2)/bnd_w);
            botshape = exp(-(D-ZZ-h/2)/bnd_w);
        case 4 % all walls
            topshape = exp( ( -ZZ+h/2)/bnd_w);
            botshape = exp(-(D-ZZ-h/2)/bnd_w);
            sdsshape = exp( ( -XX+h/2)/bnd_w) ...
                     + exp(-(L-XX-h/2)/bnd_w);
        case 5 % only walls
            sdsshape = exp( ( -XX+h/2)/bnd_w) ...
                     + exp(-(L-XX-h/2)/bnd_w);
    end
    sdsshape = max(0,sdsshape - topshape - botshape);
end

bnd_S = zeros(Nz,Nx);
bnd_C = zeros(Nz,Nx,cal.ncmp);
bnd_V = zeros(Nz,Nx);

% set specified boundaries to no slip, else to free slip
if bndmode>=4;               sds = +1;      % no slip sides for 'all sides(4)'
else;                        sds = -1; end  % free slip sides for other types
if bndmode==1 || bndmode>=3; top = +1;      % no slip top for 'top only(1)', 'top/bot(3)', 'all sides(4)'
else;                        top = -1; end  % free slip for other types
if bndmode>=2;               bot = +1;      % no slip bot for 'bot only(2)', 'top/bot(3)', 'all sides(4)'
else;                        bot = -1; end  % free slip for other types
if bndmode==5;               top = -1; bot = -1; end % free slip top/bot for 'only walls(5)'

% set ghosted index arrays
if periodic
    icx = [Nx,1:Nx,1];
    icz = [Nz,1:Nz,1];
    ifx = [Nx,1:Nx+1,2];
    ifz = [Nz,1:Nz+1,2];
else
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
    ifx = [1,1:Nx+1,Nx+1];
    ifz = [1,1:Nz+1,Nz+1];
end

% initialise solution fields
switch Tinit
    case 'layer'
        Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_T))/2 + dTr.*rp + dTg.*gp;  % potential temperature [C]
    case 'linear'
        Tp  =  T0 + (T1-T0) .* (ZZ/D) + dTr.*rp + dTg.*gp;  % potential temperature [C]
end
c = zeros(Nz,Nx,cal.ncmp);
for i = 1:cal.ncmp
    c(:,:,i)  =  c0(i) + (c1(i)-c0(i)) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_c))/2 + dcr(i).*rp + dcg(i).*gp;  % trace elements
end

te = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    te(:,:,i)  =  te0(i) + (te1(i)-te0(i)) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_c))/2 + dter(i).*rp + dteg(i).*gp;  % trace elements
end
ir = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir
    ir(:,:,i)  =  ir0(i) + (ir1(i)-ir0(i)) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_c))/2 + dirr(i).*rp + dirg(i).*gp;  % isotope ratios  
end

% apply initial boundary layers
if any(topinit(:)) && ~isnan(Twall(1)); Tp = Tp + (Twall(1)-Tp).*topinit; end
if any(botinit(:)) && ~isnan(Twall(2)); Tp = Tp + (Twall(2)-Tp).*botinit; end
if any(sdsinit(:)) && ~isnan(Twall(3)); Tp = Tp + (Twall(3)-Tp).*sdsinit; end
Tin = Tp;

for i = 1:cal.ncmp
    if any(topinit(:)) && ~any(isnan(cwall(1,:))); c(:,:,i) = c(:,:,i) + (cwall(1,i)-c(:,:,i)).*topinit; end
    if any(botinit(:)) && ~any(isnan(cwall(2,:))); c(:,:,i) = c(:,:,i) + (cwall(2,i)-c(:,:,i)).*botinit; end
    if any(sdsinit(:)) && ~any(isnan(cwall(3,:))); c(:,:,i) = c(:,:,i) + (cwall(3,i)-c(:,:,i)).*sdsinit; end
end
cin = c;

for i = 1:cal.nte
    if any(topinit(:)) && ~isnan(tewall(1,i)); te(:,:,i) = te(:,:,i) + (tewall(1,i)-te(:,:,i)).*topinit; end
    if any(botinit(:)) && ~isnan(tewall(2,i)); te(:,:,i) = te(:,:,i) + (tewall(2,i)-te(:,:,i)).*botinit; end
    if any(sdsinit(:)) && ~isnan(tewall(3,i)); te(:,:,i) = te(:,:,i) + (tewall(3,i)-te(:,:,i)).*sdsinit; end
end
tein = te; 

for i = 1:cal.nir
    if any(topinit(:)) && ~isnan(irwall(1,i)); ir(:,:,i) = ir(:,:,i) + (irwall(1,i)-ir(:,:,i)).*topinit; end
    if any(botinit(:)) && ~isnan(irwall(2,i)); ir(:,:,i) = ir(:,:,i) + (irwall(2,i)-ir(:,:,i)).*botinit; end
    if any(sdsinit(:)) && ~isnan(irwall(3,i)); ir(:,:,i) = ir(:,:,i) + (irwall(3,i)-ir(:,:,i)).*sdsinit; end
end
irin = ir;

U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wf = 0.*W; wx = 0.*W; wm = 0.*W;
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*Tp; %Div_rhoV = 0.*P;  DD = sparse(length(P(:)),length([W(:);U(:)]));
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wf  = W;  Uf  = U; 
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref = 1e-6;  
Div_V  = 0.*Tp;  advn_rho = 0.*Tp;  drhodt = 0.*Tp;  drhodto = drhodt;
exx    = 0.*Tp;  ezz = 0.*Tp;  exz = zeros(Nz-1,Nx-1);  eII = 0.*Tp;  
txx    = 0.*Tp;  tzz = 0.*Tp;  txz = zeros(Nz-1,Nx-1);  tII = 0.*Tp; 
VolSrc = 0.*Tp; 
kW     = 0.*Tp;
Tref   = min(cal.T0)+273.15;
Pref   = 1e5;
rhom   = mean(cal.rhox0-500).*ones(size(Tp)); 
rhox   = mean(cal.rhox0).*ones(size(Tp));
rhof   = cal.rhof0.*ones(size(Tp));
Pchmb  = Pchmb0;  Pchmbo = Pchmb;  Pchmboo = Pchmbo;  dPchmbdt = Pchmb;  dPchmbdto = dPchmbdt; dPchmbdtoo = dPchmbdto;  upd_Pchmb = dPchmbdt;
Pt     = Ptop + Pchmb + mean(rhom,'all').*g0.*ZZ;  Pl = Pt;
rhof   = rhof.*(1+bPf.*(Pt-Pref));
rhox   = rhox.*(1+bPx.*(Pt-Pref));
rhom   = rhom.*(1+bPm.*(Pt-Pref));
rho    = rhom;
rhofz  = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhofx  = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
rhoWo  = rhofz.*W(:,2:end-1); rhoWoo = rhoWo; advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhofx.*U(2:end-1,:); rhoUoo = rhoUo; advn_mx = 0.*rhoUo;
rhoref = mean(rho,'all');
T      =  (Tp+273.15).*exp(aT./rhoref./cP.*(Pt-Pref));
fq     = zeros(size(Tp));  mq = ones(size(Tp));  xq = 1-mq-fq; 
cmq    = c; cxq = c;  cfq = 0.*c;  cfq(:,:,end) = 1;  cf = cfq;
cm_oxd = reshape(reshape(cmq,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

% get volume fractions and bulk density
step    = 0;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;
res  = 1;  tol = 1e-12;  it = 1;
while res > tol
    Pti = Pt; Ti = T; xi = xq; fi = fq;
    
    rhoref = mean(rho,'all');
    Adbt   = aT./rhoref;
    rhofz  = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
    rhofx  = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
    if Nz==1; Pt = Ptop.*ones(size(Tp)) + Pchmb; else
        Pl(1,:)     = repmat(mean(rhofz(1,:),2).*g0.*h/2,1,Nx) + Ptop;
        Pl(2:end,:) = Pl(1,:) + repmat(cumsum(mean(rhofz(2:end-1,:),2).*g0.*h),1,Nx);
        Pt          = Pl + Pchmb;
    end
    
    T  =  (Tp+273.15).*exp(Adbt./cP.*(Pt-Pref));

    eqtime = tic;

    var.c      = reshape(c,Nx*Nz,cal.ncmp);   % component fractions [wt]
    var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
    var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
    var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt]
    var.f      = reshape(fq,Nx*Nz,1);         % bubble fraction [wt]
    var.H2O    = reshape(c(:,:,end),Nx*Nz,1); % water concentration [wt]
    var.SiO2m  = reshape(cm_oxd(:,:,1)./sum(cm_oxd(:,:,1:end-1),3),Nx*Nz,1); % melt silica concentration [wt]
    cal.H2Osat = fluidsat(var.T,var.P,var.SiO2m,cal);

    [var,cal] = meltmodel(var,cal,'E');


    mq = reshape(var.m,Nz,Nx);
    fq = reshape(var.f,Nz,Nx);
    xq = reshape(var.x,Nz,Nx);
    x  = xq;  m = mq;  f = fq;

    cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
    cmq = reshape(var.cm,Nz,Nx,cal.ncmp);
    cm  = cmq; cx = cxq;

    cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

    eqtime = toc(eqtime);
    EQtime = EQtime + eqtime;

    update;

    res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2) ...
         + norm( T(:)- Ti(:),2)./norm( T(:),2) ...
         + norm((x(:)- xi(:)).*(x(:)>TINY^0.5),2)./(norm(x(:),2)+TINY) ...
         + norm((f(:)- fi(:)).*(f(:)>TINY^0.5),2)./(norm(f(:),2)+TINY);

    it = it+1;
end
rhoo = rho;
dto  = dt; 

m0  = mean(m(:)); 
x0  = mean(x(:)); 
f0  = mean(f(:)); 
cx0 = squeeze(mean(mean(cx,1),2)).'; 
cm0 = squeeze(mean(mean(cm,1),2)).'; 

cm0_mem = cm0*cal.cmp_mem;
cx0_mem = cx0*cal.cmp_mem;
 c0_mem =  c0*cal.cmp_mem;

cm0_oxd = cm0*cal.cmp_oxd;
cx0_oxd = cx0*cal.cmp_oxd;
 c0_oxd =  c0*cal.cmp_oxd;

rhof0 = cal.rhof0;
rhox0 = mean(rhox(:));

cm0_oxd_all = zeros(1,9);
cm0_oxd_all(:,cal.ioxd) = cm0_oxd;
etam0 = Giordano08(cm0_oxd_all,T0);
rhom0 = DensityX(cm0_oxd_all,T0,Ptop);

cm1_oxd = (0.99.*cm0_mem(1,:) + 0.01.*cal.cmp_mem(end-1,:))*cal.mem_oxd/100;
cm2_oxd = (0.99.*cm0_mem(1,:) - 0.01.*cal.cmp_mem(end-1,:))*cal.mem_oxd/100;

cm1_oxd_all = zeros(1,9);
cm1_oxd_all(:,cal.ioxd) = cm1_oxd;
cm2_oxd_all = zeros(1,9);
cm2_oxd_all(:,cal.ioxd) = cm2_oxd;

rho0 = (x0./rhox0 + m0./rhom0).^-1;

fprintf('    initial T   : %4.3f \n'  ,T0);
fprintf('    initial SiO2: %4.3f \n'  ,c0_oxd(1)./sum(c0_oxd(1:end-1)).*100);
fprintf('    initial H2O : %4.3f \n'  ,c0_oxd(end));
fprintf('    initial x   : %4.3f \n'  ,x0);
fprintf('    initial f   : %4.3f \n\n',f0);

% rhom1 = DensityX(cm1_oxd_all,T0,Ptop/1e8);
% rhom2 = DensityX(cm2_oxd_all,T0,Ptop/1e8);
% 
% DrhoT = rhom0.*aT*max([abs(T0-Twall)/100,abs(T0-T1),T0/100]);
% Drhoc = abs(rhom1-rhom2);
% Drhox = 0.01*abs(rhox0-rhom0);
% Drhof = 0.01*abs(cal.rhof0-rhom0) * (max([c0(end),c1(end),max(cwall(:,end))])>TINY);
% Drho0 = DrhoT + Drhoc + Drhox + Drhof;
% 
% uT    = DrhoT*g0*(D/10)^2/etam0;
% uc    = Drhoc*g0*(D/10)^2/etam0;
% ux    = Drhox*g0*(D/10)^2/etam0;
% uf    = Drhof*g0*(D/10)^2/etam0 * (max([c0(end),c1(end),max(cwall(:,end))])>TINY);
% u0    = Drho0*g0*(D/10)^2/etam0;
% 
% wx0   = abs(rhox0-rhom0)*g0*dx0^2/etam0;
% wf0   = abs(rhof0-rhom0)*g0*df0^2/etam0 * (max([c0(end),c1(end),max(cwall(:,end))])>TINY);
% 
% ud0   = kT0/rhom0/cP/(D/10);
% 
% uwT   = bnd_w/tau_T; 
% uwc   = bnd_w/tau_a; 
% 
% RaT   = uT/ud0;
% Rac   = uc/ud0;
% Rax   = ux/ud0;
% Raf   = uf/ud0;
% Ra    = u0/ud0;
% 
% Rux   = wx0/u0;
% Ruf   = wf0/u0;
% 
% RwT   = uwT/u0;
% Rwc   = uwc/u0;
% 
% Re    = u0*rhom0*(D/10)/etam0;
% Rex   = wx0*rhom0*dx0/etam0;
% Ref   = wf0*rhom0*df0/etam0;
% 
% fprintf('    crystal Re: %1.3e \n'  ,Rex);
% fprintf('     bubble Re: %1.3e \n'  ,Ref);
% fprintf('     system Re: %1.3e \n\n',Re );
% 
% fprintf('    thermal Ra: %1.3e \n'  ,RaT);
% fprintf('   chemical Ra: %1.3e \n'  ,Rac);
% fprintf('    crystal Ra: %1.3e \n'  ,Rax);
% fprintf('     bubble Ra: %1.3e \n'  ,Raf);
% fprintf('   combined Ra: %1.3e \n\n',Ra );
% 
% fprintf('    crystal Ru: %1.3e \n'  ,Rux);
% fprintf('     bubble Ru: %1.3e \n\n',Ruf);
% 
% fprintf('    thermal Rw: %1.3e \n'  ,RwT);
% fprintf('   chemical Rw: %1.3e \n\n',Rwc);

% get bulk enthalpy, silica, volatile content densities
Tp   = Tp+273.15;
S    = rho.*(cP.*log(T/Tref) + x.*Dsx + f.*Dsf - Adbt.*(Pt-Pref));  So = S;  res_S = 0.*S;
S0   = rho.*(cP.*log(  Tref) + x.*Dsx + f.*Dsf - Adbt.*    Pref );  
C    = rho.*(m.*cm + x.*cx + f.*cf); Co = C;  res_C = 0.*C;
X    = rho.*x; Xo = X;  res_X = 0.*X;
F    = rho.*f; Fo = F;  res_F = 0.*F;
M    = rho.*m; Mo = M;  res_M = 0.*M;
RHO  = X+M+F;

% get phase entropies
sm = S./rho - x.*Dsx - f.*Dsf;
sx = sm + Dsx;
sf = sm + Dsf;

% get trace element phase compositions
Kte = zeros(Nz,Nx,cal.nte);
tem = zeros(Nz,Nx,cal.nte);
tex = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    for j=1:cal.nmem; Kte(:,:,i) = Kte(:,:,i) + cal.Kte_mem(i,j) .* c_mem(:,:,j)./100; end

    tem(:,:,i)  = te(:,:,i)./(m + x.*Kte(:,:,i));
    tex(:,:,i)  = te(:,:,i)./(m./Kte(:,:,i) + x);
end

irm = zeros(Nz,Nx,cal.nir);
irx = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir
    irm(:,:,i)  = ir(:,:,i)./(1-f);
    irx(:,:,i)  = ir(:,:,i)./(1-f);
end

% get geochemical component densities
TE = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    TE(:,:,i)  = rho.*(m.*tem(:,:,i) + x.*tex(:,:,i));
end
TEo = TE;  TEi = TE;

IR = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir
    IR(:,:,i)  = rho.*(m.*irm(:,:,i) + x.*irx(:,:,i));
end
IRo = IR;  IRi = IR;

% initialise phase change rates
Gx  = 0.*x;  Gf  = 0.*f;  Gm  = 0.*m;

% initialise auxiliary variables 
dSdt   = 0.*T;  dSdto  = dSdt; diss_h = 0.*T;
dCdt   = 0.*c;  dCdto  = dCdt;
dFdt   = 0.*f;  dFdto  = dFdt;
dXdt   = 0.*x;  dXdto  = dXdt;
dMdt   = 0.*m;  dMdto  = dMdt;
bnd_TE = zeros(Nz,Nx,cal.nte);
adv_TE = zeros(Nz,Nx,cal.nte);
dff_TE = zeros(Nz,Nx,cal.nte);
Kte    = zeros(Nz,Nx,cal.nte);
dTEdt  = 0.*te; dTEdto = dTEdt;
bnd_IR = zeros(Nz,Nx,cal.nir);
adv_IR = zeros(Nz,Nx,cal.nir);
dff_IR = zeros(Nz,Nx,cal.nir);
dIRdt  = 0.*ir; dIRdto = dIRdt;
upd_S  = 0.*S;
upd_C  = 0.*C;
upd_X  = 0.*X;
upd_F  = 0.*F;
upd_M  = 0.*M;
upd_rho= 0.*rho;
upd_TE = 0.*TE;
upd_IR = 0.*IR;

% initialise timing and iterative parameters
frst    = 1;
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumMdt = 0; dsumMdto = 0;
dsumSdt = 0; dsumSdto = 0;
dsumCdt = 0; dsumCdto = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','c','cm','cx','cf','TE','IR','te','ir','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTEdt','dIRdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
        name = [opdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        SOL = [W(:);U(:);P(:)];
        RHO = X+M+F;
        Tp = Tref*exp((S - X.*Dsx - F.*Dsf)./RHO./cP);
        sm = (S - X.*Dsx - F.*Dsf)./RHO;
        sx = sm + Dsx;
        sf = sm + Dsf;
        
        So = S;
        Co = C;
        Xo = X;
        Fo = F;
        Mo = M;
        rhoo = rho;
        TEo = TE;
        IRo = IR;
        dSdto = dSdt;
        dCdto = dCdt;
        dXdto = dXdt;
        dFdto = dFdt;
        dMdto = dMdt;
        drhodto = drhodt;
        dTEdto = dTEdt;
        dIRdto = dIRdt;

        update; output;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        fluidmech;
        update;
        history;
        output;
    end
else
    % complete, plot, and save initial condition
    fluidmech;
    update;
    history;
    output;
end

restart = 0;
time    = time+dt;
step    = step+1;
