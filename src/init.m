%% *****  MODEL INITIALISATION  *******************************************

if ~postprc
% create output directory
if ~isfolder([outdir,'/',runID])
    mkdir([outdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 && save_op == 1 && ~postprc
    parfile = [outdir,'/',runID,'/',runID,'_par'];
    save(parfile);
end

fprintf('\n\n')
fprintf('*************************************************************\n');
fprintf('*****  RUN NAKHLA MODEL | %s  *************\n',datetime('now'));
fprintf('*************************************************************\n');
fprintf('\n   run ID: %s \n\n',runID);
end

run(['../cal/cal_',calID]);  % load melt model calibration

% load and process custom colormaps
switch colourmap
    case 'ocean'
        load ./colmap/ocean.mat
        colmap = ocean;   
    case 'batlow'
        load ./colmap/batlow.mat 
        colmap = batlow;
    case 'batlowW'
        load ./colmap/batlowW.mat
        colmap = batlowW;
    case 'batlowK'
        load ./colmap/batlowK.mat
        colmap = batlowK;
    case 'lipari'
        load ./colmap/lipari.mat
        colmap = lipari;  
    case 'glasgow'
        load ./colmap/glasgow.mat;
        colmap = glasgow; 
    case 'navia'
        load ./colmap/navia.mat;
        colmap = navia;  
    case 'lapaz'
        load ./colmap/lapaz.mat;
        colmap = lapaz;  
    case 'lajolla'
        load ./colmap/lajolla.mat;
        colmap = lajolla;  
end
nclmp = length(colmap);

colmapcmp = colmap(max(1,min(nclmp,1+round((1:nclmp)/(nclmp/(cal.ncmp-1)))*round(nclmp/(cal.ncmp-1)))),:);
colmapoxd = colmap(max(1,min(nclmp,1+round((1:nclmp)/(nclmp/(cal.noxd-1)))*round(nclmp/(cal.noxd-1)))),:);
colmapmem = colmap(max(1,min(nclmp,1+round((1:nclmp)/(nclmp/(cal.nmem-1)))*round(nclmp/(cal.nmem-1)))),:);
colmapmsy = colmap(max(1,min(nclmp,1+round((1:nclmp)/(nclmp/(cal.nmsy-1)))*round(nclmp/(cal.nmsy-1)))),:);

if periodic % periodic sides
    BCA     =  {'closed','periodic'};  % boundary condition on advection (top/bot, sides)
    BCD     =  {'closed','periodic'};  % boundary condition on advection (top/bot, sides)
else % closed sides
    BCA     =  {'closed','closed'};  % boundary condition on advection (top/bot, sides)
    BCD     =  {'closed','closed'};  % boundary condition on advection (top/bot, sides) 
end

Dsx = -cal.Dsx;
Dsf =  cal.Dsf;

% normalise major components to anhydrous unit sum, rescale to hydrous
c0(1:end-1) = c0(1:end-1)./sum(c0(1:end-1)).*(1-c0(end));
c1(1:end-1) = c1(1:end-1)./sum(c1(1:end-1)).*(1-c1(end));
cwall(:,1:end-1) = cwall(:,1:end-1)./sum(cwall(:,1:end-1),2).*(1-cwall(:,end));
dcg   = dcg-round(mean(dcg),16);
dcr   = dcr-round(mean(dcr),16);

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

gp = exp(-(XX-L/2  ).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2) ...
   + exp(-(XX-L/2+L).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2) ...
   + exp(-(XX-L/2-L).^2/(max(L,D)/8)^2 - (ZZ-D/2).^2/(max(L,D)/8)^2);

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
            topshape = exp( ( -ZZ)/max(h,bnd_w));
        case 2  % bot only
            botshape = exp(-(D-ZZ)/max(h,bnd_w));
        case 3  % top/bot only
            topshape = exp( ( -ZZ)/max(h,bnd_w));
            botshape = exp(-(D-ZZ)/max(h,bnd_w));
        case 4 % all walls
            topshape = exp( ( -ZZ)/max(h,bnd_w));
            botshape = exp(-(D-ZZ)/max(h,bnd_w));
            sdsshape = exp( ( -XX)/max(h,bnd_w)) ...
                     + exp(-(L-XX)/max(h,bnd_w));
        case 5 % only walls
            sdsshape = exp( ( -XX)/max(h,bnd_w)) ...
                     + exp(-(L-XX)/max(h,bnd_w));
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
if periodic  % periodic side boundaries
    icx = [Nx,1:Nx,1];
    icz = [1,1:Nz,Nz];
    ifx = [Nx,1:Nx+1,2];
    ifz = [2,1:Nz+1,Nz];
else         % closed side boundaries
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
    ifx = [2,1:Nx+1,Nx];
    ifz = [2,1:Nz+1,Nz];
end

% initialise solution fields
switch init_mode
    case {'constant' , 'liquidus'}
        Tp  =  max(cal.T0) + 0.*rp;  % potential temperature [C]
        c = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
            c(:,:,i)  =  c0(i) + dcr(i).*rp + dcg(i).*gp;  % trace elements
        end
        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trc(:,:,i)  =  trc0(i) + dr_trc(i).*rp + dg_trc(i).*gp;  % trace elements
        end
    case 'layer'
        Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay+rp*dlay)/wlay_T))/2 + dTr.*rp + dTg.*gp;  % potential temperature [C]
        c = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
            c(:,:,i)  =  c0(i) + (c1(i)-c0(i)) .* (1+erf((ZZ/D-zlay+rp*dlay)/wlay_c))/2 + dcr(i).*rp + dcg(i).*gp;  % trace elements
        end
        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trc(:,:,i)  =  trc0(i) + (trc1(i)-trc0(i)) .* (1+erf((ZZ/D-zlay+rp*dlay)/wlay_c))/2 + dr_trc(i).*rp + dg_trc(i).*gp;  % trace elements
        end
    case 'linear'
        Tp  =  T0 + (T1-T0) .* (ZZ/D) + dTr.*rp + dTg.*gp;  % potential temperature [C]
        c = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
            c(:,:,i)  =  c0(i) + (c1(i)-c0(i)) .* (ZZ/D) + dcr(i).*rp + dcg(i).*gp;  % trace elements
        end
        trc = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trc(:,:,i)  =  trc0(i) + (trc1(i)-trc0(i)) .* (ZZ/D) + dr_trc(i).*rp + dg_trc(i).*gp;  % trace elements
        end
    case 'read_1D'
        initname = [outdir,'/',runID,'/',runID,'_init.mat'];
        load(initname,'m','Tp','c','trc');

        hi  = D./(size(Tp,1));
        Xci = +hi/2:hi:L-hi/2;
        Zci = +hi/2:hi:D-hi/2;
        m   = repmat(interp1(Zci,m,Zc).',1,Nx);
        indmix = m>0.75;

        Tp  = repmat(interp1(Zci,Tp-273.15,Zc).',1,Nx);
        Tp(indmix) = mean(Tp(indmix));
        Tp = Tp + dTr.*rp + dTg.*gp;

        ci = zeros(Nz,Nx,cal.ncmp);
        for i = 1:cal.ncmp
            cii  = repmat(interp1(Zci,c(:,:,i),Zc).',1,Nx);
            cii(indmix) = mean(cii(indmix));
            ci(:,:,i)   = cii + dcr(i).*rp + dcg(i).*gp;
        end
        c = ci;

        trci = zeros(Nz,Nx,cal.ntrc);
        for i = 1:cal.ntrc
            trcii   = repmat(interp1(Zci,trc(:,:,i),Zc).',1,Nx);
            trcii(indmix) = mean(trcii(indmix));
            trci(:,:,i)   = trcii + dr_trc(i).*rp + dg_trc(i).*gp;
        end
        trc = trci;
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

for i = 1:cal.ntrc
    if any(topinit(:)) && ~isnan(trcwall(1,i)); trc(:,:,i) = trc(:,:,i) + (trcwall(1,i)-trc(:,:,i)).*topinit; end
    if any(botinit(:)) && ~isnan(trcwall(2,i)); trc(:,:,i) = trc(:,:,i) + (trcwall(2,i)-trc(:,:,i)).*botinit; end
    if any(sdsinit(:)) && ~isnan(trcwall(3,i)); trc(:,:,i) = trc(:,:,i) + (trcwall(3,i)-trc(:,:,i)).*sdsinit; end
end
tein = trc;

U   =  zeros(Nz+2,Nx+1);  UBG = U; Ui = U; upd_U = 0*U;
W   =  zeros(Nz+1,Nx+2);  WBG = W; Wi = W; wf = 0.*W; wx = 0.*W; wm = 0.*W; upd_W = 0*W;
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*Tp; upd_P = 0*P; %Div_rhoV = 0.*P;  DD = sparse(length(P(:)),length([W(:);U(:)]));
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wf  = W;  Uf  = U; 
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

Re     = eps;  
Div_V  = 0.*Tp;  advn_rho = 0.*Tp;  advn_X = 0.*Tp; advn_M = 0.*Tp; advn_F = 0.*Tp; drhodt = 0.*Tp;  drhodto = drhodt;
exx    = 0.*Tp;  ezz = 0.*Tp;  exz = zeros(Nz-1,Nx-1);  eII = 0.*Tp;  
txx    = 0.*Tp;  tzz = 0.*Tp;  txz = zeros(Nz-1,Nx-1);  tII = 0.*Tp; 
eta    = ones(Nz,Nx);
dV = 0.*Tp; 
ke     = 0.*Tp;
Tref   = min(cal.T0)+273.15;
Pref   = 1e5;
sref   = 0e3;
c0_oxd = c0*cal.cmp_oxd;
c0_oxd_all = zeros(size(c0,1),9);
if cal.noxd>9
    c0_oxd_all = c0_oxd(:,cal.ioxd);
else
    c0_oxd_all(:,cal.ioxd) = c0_oxd;
end
rhom0  = mean(cal.rhox0-500).*ones(size(Tp)); 
rhox0  = mean(cal.rhox0).*ones(size(Tp));
rhof0  = cal.rhof0.*ones(size(Tp));
Pchmb  = Pchmb0;  Pchmbo = Pchmb;  Pchmboo = Pchmbo;  dPchmbdt = Pchmb;  dPchmbdto = dPchmbdt; dPchmbdtoo = dPchmbdto;  upd_Pchmb = dPchmbdt;
Pt     = Ptop + Pchmb + mean(rhom0,'all').*g0.*ZZ;  Pl = Pt;  Pto = Pt; Ptoo = Pt; dPtdt = 0*Pt; dPtdto = dPtdt; dPtdtoo = dPtdto;
rhof   = rhof0.*(1+bPf.*(Pt-Pref));
rhox   = rhox0.*(1+bPx.*(Pt-Pref));
rhom   = rhom0.*(1+bPm.*(Pt-Pref));
rho    = rhom;
rhow   = (rho(icz(1:end-1),:)+rho(icz(2:end),:))/2;
rhou   = (rho(:,icx(1:end-1))+rho(:,icx(2:end)))/2;
rhoWo  = rhow.*W(:,2:end-1); rhoWoo = rhoWo; advn_mz = 0.*rhoWo(2:end-1,:);
rhoUo  = rhou.*U(2:end-1,:); rhoUoo = rhoUo; advn_mx = 0.*rhoUo;
fq     = zeros(size(Tp));  mq = ones(size(Tp))/2;  xq = 1-mq-fq; 
cmq    = c; cxq = c;  cfq = 0.*c;  cfq(:,:,end) = 1;  cf = cfq;
cm_oxd = reshape(reshape(c,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
cm_oxd_all = zeros(Nz,Nx,9);
if cal.noxd>9
    cm_oxd_all = cm_oxd(:,:,cal.ioxd);
else
    cm_oxd_all(:,:,cal.ioxd) = cm_oxd;
end
aT   = aTm;
kT   = kTm;
cP   = cPm; RhoCp = rho.*cP;
Adbt = aT./RhoCp;
Tp   = (Tp+273.15); %T = Tp;
T    = Tp.*exp(Adbt.*(Pt-Pref));
sm   = cPm.*log(Tp./Tref);  sx = cPx.*log(Tp./Tref) + Dsx;  sf = cPf.*log(Tp./Tref) + Dsf;
x    = xq;  m = mq;  f = fq; mu = m; chi = x; phi = f;

% get volume fractions and bulk density
step    = 0;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
a1      = 1; a2 = 0; a3 = 0; b1 = 1; b2 = 0; b3 = 0;

if ~postprc
res  = 1;  tol = 1e-12;  it = 1;
while res > tol
    Ptii = Pt; Ti = T; xi = xq; fi = fq;

    eqtime = tic;

    var.c      = reshape(c,Nx*Nz,cal.ncmp);   % component fractions [wt]
    var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
    var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
    var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt]
    var.f      = reshape(fq,Nx*Nz,1);         % bubble fraction [wt]
    var.H2O    = var.c(:,end);                % water concentration [wt]
    var.X      = reshape(cm_oxd_all,Nz*Nx,9); % melt oxide fractions [wt %]
    cal.H2Osat = fluidsat(var); % water saturation [wt]

    [var,cal] = leappart(var,cal,'E');

    Tsol   = reshape(cal.Tsol,Nz,Nx);
    Tliq   = reshape(cal.Tliq,Nz,Nx);
    H2Osat = reshape(cal.H2Osat,Nz,Nx);

    mq = reshape(var.m.*(var.m>eps^0.5),Nz,Nx);
    fq = reshape(var.f.*(var.f>eps^0.5),Nz,Nx);
    xq = reshape(var.x.*(var.x>eps^0.5),Nz,Nx);
    mq = mq./(mq+xq+fq);
    xq = xq./(mq+xq+fq);
    fq = fq./(mq+xq+fq);
    x  = xq;  m = mq;  f = fq;

    cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
    cfq = reshape(var.cf,Nz,Nx,cal.ncmp);
    cmq = reshape(var.cm,Nz,Nx,cal.ncmp);
    cm  = cmq; cx = cxq;  cf = cfq;

    eqtime = toc(eqtime);
    EQtime = EQtime + eqtime;

    update;

    % T  = Tp.*exp(Adbt.*(Pt-Pref));
    % sm = cPm.*log(Tp./T0);  sx = sm+Dsx;  sf = sx+Dsf;

    X    = rho.*x; Xo = X;  res_X = 0.*X;
    F    = rho.*f; Fo = F;  res_F = 0.*F;
    M    = rho.*m; Mo = M;  res_M = 0.*M;
    RHO  = X+M+F;
    C    = M.*cm + X.*cx + F.*cf;

    switch init_mode
        case 'liquidus'
            T  = Tliq+273.15+T0;
            a  = aTm.*(T -Tref);
            b  = bPm.*(Pt-Pref);
            sm = sref     + cPm.*log(T/Tref) - (aTm./(bPm.*rhom0)) .* log((1-a+b)./(1-a));
            a  = aTx.*(T -Tref);
            b  = bPx.*(Pt-Pref);
            sx = sref+Dsx + cPx.*log(T/Tref) - (aTx./(bPx.*rhox0)) .* log((1-a+b)./(1-a));
            a  = aTf.*(T -Tref);
            b  = bPf.*(Pt-Pref);
            sf = sref+Dsf + cPf.*log(T/Tref) - (aTf./(bPf.*rhof0)) .* log((1-a+b)./(1-a));
            S  = M.*sm + X.*sx + F.*sf;
            [Tp,~ ] = StoT(Tp,S./rho,Pref+0*Pt,cat(3,m,x,f),[cPm;cPx;cPf],[aTm;aTx;aTf],[bPm;bPx;bPf],cat(3,rhom0,rhox0,rhof0),[sref;sref+Dsx;sref+Dsf],Tref,Pref);
            [T ,si] = StoT(T ,S./rho,       Pt,cat(3,m,x,f),[cPm;cPx;cPf],[aTm;aTx;aTf],[bPm;bPx;bPf],cat(3,rhom0,rhox0,rhof0),[sref;sref+Dsx;sref+Dsf],Tref,Pref);
            sm = si(:,:,1); sx = si(:,:,2); sf = si(:,:,3);
        otherwise
            sm = cPm.*log(Tp./Tref);  sx = cPx.*log(Tp./Tref) + Dsx;  sf = cPf.*log(Tp./Tref) + Dsf;
            S  = M.*sm + X.*sx + F.*sf;

            [Tp,~ ] = StoT(Tp,S./rho,Pref+0*Pt,cat(3,m,x,f),[cPm;cPx;cPf],[aTm;aTx;aTf],[bPm;bPx;bPf],cat(3,rhom0,rhox0,rhof0),[sref;sref+Dsx;sref+Dsf],Tref,Pref);
            [T ,si] = StoT(T ,S./rho,       Pt,cat(3,m,x,f),[cPm;cPx;cPf],[aTm;aTx;aTf],[bPm;bPx;bPf],cat(3,rhom0,rhox0,rhof0),[sref;sref+Dsx;sref+Dsf],Tref,Pref);
            sm = si(:,:,1); sx = si(:,:,2); sf = si(:,:,3);
    end

    res  = norm(Pt(:)-Ptii(:),2)./norm(Pt(:),2) ...
         + norm( T(:)-Ti  (:),2)./norm( T(:),2);

    it = it+1;
end
% end
To   = T;  
Tpo  = Tp;
So   = S;  
Mo   = M;
Co   = C;
Xo   = X;
Fo   = F;
rhoo = rho;
dto  = dt; 

% initialise correlation length for convective/turbulent regularisation
Delta_cnv0 = Delta_cnv;
corrl;

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
if cal.noxd>9
    cm0_oxd_all = cm0_oxd(:,cal.ioxd);
else
    cm0_oxd_all(:,cal.ioxd) = cm0_oxd;
end
etam0 = Giordano08(cm0_oxd_all,T0);
rhom0 = DensityX(cm0_oxd_all,T0,Ptop);

cm1_oxd = (0.99.*cm0_mem(1,:) + 0.01.*cal.cmp_mem(end-1,:))*cal.mem_oxd/100;
cm2_oxd = (0.99.*cm0_mem(1,:) - 0.01.*cal.cmp_mem(end-1,:))*cal.mem_oxd/100;

cm1_oxd_all = zeros(1,9);
cm2_oxd_all = zeros(1,9);
if cal.noxd>9
    cm1_oxd_all = cm1_oxd(:,cal.ioxd);
    cm2_oxd_all = cm2_oxd(:,cal.ioxd);
else
    cm1_oxd_all(:,cal.ioxd) = cm1_oxd;
    cm2_oxd_all(:,cal.ioxd) = cm2_oxd;
end

rho0 = (x0./rhox0 + m0./rhom0).^-1;

fprintf('    initial T   : %4.3f \n'  ,T0);
fprintf('    initial SiO2: %4.3f \n'  ,c0_oxd(1)./sum(c0_oxd(1:end-1)).*100);
fprintf('    initial H2O : %4.3f \n'  ,c0_oxd(end));
fprintf('    initial x   : %4.3f \n'  ,x0);
fprintf('    initial f   : %4.3f \n\n',f0);

% get trace element phase compositions
Ktrc = zeros(Nz,Nx,cal.ntrc);
trcm = zeros(Nz,Nx,cal.ntrc);
trcx = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* c_mem(:,:,j)./100; end

    trcm(:,:,i)  = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcx(:,:,i)  = trc(:,:,i)./(m./Ktrc(:,:,i) + x);
end

% get geochemical component densities
TRC = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    TRC(:,:,i)  = rho.*(m.*trcm(:,:,i) + x.*trcx(:,:,i));
end
TRCo = TRC;

% initialise phase change rates
Gx  = 0.*x;  Gf  = 0.*f;  Gm  = 0.*m;

% initialise auxiliary variables 
dSdt   = 0.*T;  dSdto  = dSdt; diss_h = 0.*T; 
dTpdt  = 0.*T;  dTpdto = dTpdt;
dTdt   = 0.*T;  dTdto  = dTdt;
dCdt   = 0.*c;  dCdto  = dCdt;
dFdt   = 0.*f;  dFdto  = dFdt;
dXdt   = 0.*x;  dXdto  = dXdt;
dMdt   = 0.*m;  dMdto  = dMdt;
bnd_TRC = zeros(Nz,Nx,cal.ntrc);
adv_TRC = zeros(Nz,Nx,cal.ntrc);
dff_TRC = zeros(Nz,Nx,cal.ntrc);
K_trc     = zeros(Nz,Nx,cal.ntrc);
dTRCdt  = 0.*trc; dTRCdto = dTRCdt;
upd_S  = 0.*S;
upd_Tp = 0.*Tp;
upd_T  = 0.*T;
upd_C  = 0.*C;
upd_X  = 0.*X;
upd_F  = 0.*F;
upd_M  = 0.*M;
upd_rho= 0.*rho;
upd_TRC = 0.*TRC;
upd_eta = 0.*eta;

% initialise timing and iterative parameters
frst    = 1;
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumSdt = 0; dsumSdto = 0;
dsumBdt = 0; dsumBdto = 0;
dsumMdt = 0; dsumMdto = 0;
dsumXdt = 0; dsumXdto = 0;
dsumFdt = 0; dsumFdto = 0;
dsumCdt = 0; dsumCdto = 0;
dsumTdt = 0; dsumTdto = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [outdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [outdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','Pchmb','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','Tp','c','cm','cx','cf','sm','sx','sf','TRC','trc','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTRCdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','dV','wf','wx','wm','cal');
        name = [outdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        SOL = [W(:);U(:);P(:)];
        RHO = X+M+F;

        update;
        phseql;
        corrl;
        update;
        store;
        fluidmech;
        update;
        output;

        time    = time+dt;
        step    = step+1;
        restart = 0;

    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        restart = 0;
        store;
        fluidmech;
        update;
        history;
        output;
        step = step+1;
    end
else
    % complete, plot, and save initial condition
    store;
    fluidmech;
    update;
    history;
    output;
    step = step+1;
end
else
    % M   = m.*rho;
    % X   = x.*rho;
    % F   = f.*rho;
    Delta_cnv0 = Delta_cnv;
    RHO = rho;
end

