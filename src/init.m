load ocean;  % load custom colormap
run(['../cal/cal_',calID]);

% minimum cutoff phase, component fractions
TINY     =  1e-16;               

% get coordinate arrays
Xc        = -h/2:h:L+h/2;
Zc        = -h/2:h:D+h/2;
[XX,ZZ]   = meshgrid(Xc,Zc);
Xf        = (Xc(1:end-1)+Xc(2:end))./2;
Zf        = (Zc(1:end-1)+Zc(2:end))./2;
[XXu,ZZu] = meshgrid(Xf,Zc);
[XXw,ZZw] = meshgrid(Xc,Zf);

Nx = length(Xc);
Nz = length(Zc);

% get smoothed initialisation field
rng(seed);
rp = randn(Nz,Nx);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 + diff(rp(2:end-1,:),2,2)./8;
    rp = rp - mean(mean(rp(2:end-1,2:end-1)));
    rp([1 end],:) = 0;
    rp(:,[1 end]) = 0;
end
rp = rp./max(abs(rp(:)));

% get mapping arrays
NP =  Nz   * Nx   ;
NW = (Nz-1)* Nx   ;
NU =  Nz   *(Nx-1);
MapP = reshape(1:NP,Nz  ,Nx  );
MapW = reshape(1:NW,Nz-1,Nx  );
MapU = reshape(1:NU,Nz  ,Nx-1) + NW;

if bndinit
    switch bndmode
        case 0  % none
            bndinit = zeros(size(ZZ));
        case 1  % top only
            bndinit = (1+erf( ( -ZZ+dw)/(dw/5)))/2;
        case 2  % bot only
            bndinit = (1+erf(-(D-ZZ-dw)/(dw/5)))/2;
        case 3  % top/bot only
            bndinit = (1+erf( ( -ZZ+dw)/(dw/5)))/2 ...
                    + (1+erf(-(D-ZZ-dw)/(dw/5)))/2;
        case 4 % all walls
            bndinit = (1+erf( ( -ZZ+dw)/(dw/5)))/2 ...
                    + (1+erf(-(D-ZZ-dw)/(dw/5)))/2 ...
                    + (1+erf( ( -XX+dw)/(dw/5)))/2 ...
                    + (1+erf(-(L-XX-dw)/(dw/5)))/2;
    end
end
switch bndmode
    case 0  % none
        bndshape = zeros(size(ZZ));
    case 1  % top only
        bndshape = exp( ( -ZZ+h/2)/dw);
    case 2  % bot only
        bndshape = exp(-(D-ZZ-h/2)/dw);
    case 3  % top/bot only
        bndshape = exp( ( -ZZ+h/2)/dw) ...
            + exp(-(D-ZZ-h/2)/dw);
    case 4 % all walls
        bndshape = exp( ( -ZZ+h/2)/dw) ...
            + exp(-(D-ZZ-h/2)/dw) ...
            + exp( ( -XX+h/2)/dw) ...
            + exp(-(L-XX-h/2)/dw);
end
bndshape = max(0,min(1,bndshape));
bndshape([1 end],:) = bndshape([2 end-1],:);
bndshape(:,[1 end]) = bndshape(:,[2 end-1]);

bndS = 0.*bndshape;  bndC =  0.*bndshape;  bndV =  0.*bndshape;

% set specified boundaries to no slip, else to free slip
if bndmode==4;               sds = +1;      % no slip sides for 'all sides(4)'
else;                        sds = -1; end  % free slip sides for other types
if bndmode==1 || bndmode>=3; top = +1;      % no slip top for 'top only(1)', 'top/bot(3)', 'all sides(4)'
else;                        top = -1; end  % free slip for other types
if bndmode>=2;               bot = +1;      % no slip bot for 'bot only(2)', 'top/bot(3)', 'all sides(4)'
else;                        bot = -1; end  % free slip for other types

% initialise solution fields
Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay)/wlay_T))/2 + dT.*rp;  if any(bndinit(:)) && ~isnan(Twall); Tp = Tp + (Twall-Tp).*bndinit; end % potential temperature [C]
c   =  c0 + (c1-c0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dc.*rp;  if any(bndinit(:)) && ~isnan(cwall); c  = c  + (cwall-c ).*bndinit; end % major component
v   =  v0 + (v1-v0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dv.*rp;  if any(bndinit(:)) && ~isnan(vwall); v  = v  + (vwall-v ).*bndinit; end % volatile component

it  =  it0 + (it1-it0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dit.*rp;  if any(bndinit(:)) && ~isnan(itwall); it  = it  + (itwall-it ).*bndinit; end % incompatible trace element
ct  =  ct0 + (ct1-ct0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dct.*rp;  if any(bndinit(:)) && ~isnan(ctwall); ct  = ct  + (ctwall-ct ).*bndinit; end % compatible trace element
si  =  si0 + (si1-si0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dsi.*rp;  if any(bndinit(:)) && ~isnan(siwall); si  = si  + (siwall-si ).*bndinit; end % stable isotope ratio
rip =  ri0 + (ri1-ri0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dri.*rp;  if any(bndinit(:)) && ~isnan(riwall); rip = rip + (riwall-rip).*bndinit; end % radiogenic isotope parent
rid =  rip.*HLRID./HLRIP;                                           % radiogenic isotope daughter

U   =  zeros(size((XX(:,1:end-1)+XX(:,2:end))));  Ui = U;  res_U = 0.*U;
W   =  zeros(size((XX(1:end-1,:)+XX(2:end,:))));  Wi = W;  res_W = 0.*W; wf = 0.*W; wc = 0.*W;
P   =  0.*Tp;  Pi = P;  res_P = 0.*P;  meanQ = 0;  
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wf  = W;  Uf  = U; 
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_rhoV = 0.*P;  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(Nz-1,Nx-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(Nz-1,Nx-1);  tII = 0.*P; 
VolSrc =  0.*P;  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;

rhoo =  rhom0.*ones(size(Tp)); rhoref = rhom0;  %#ok<NASGU>
dto  =  dt;
Pt   =  rhoref.*g0.*ZZ + Ptop;  
if Nx<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))).*ones(size(Pt)); end
T    =  (Tp+273.15).*exp(aT./rhoref./cP.*Pt); % real temperature [K]

% get volume fractions and bulk density
step = 0;
res  = 1;  tol = 1e-15;  x = ones(size(T))./10;  f = v/2;
while res > tol
    xi = x;  fi = f;
    
    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T-273.15,c,v,Pt,cal,TINY);
    
    x  = xq;  f = fq;  m = 1-x-f;
    cm = cmq; cx = cxq;
    vm = vmq; vf = vfq;
    Kc = cxq./cmq;
    Kf = vfq./vmq;

    update;
    
    rhoref  = mean(mean(rho(2:end-1,2:end-1)));
    Pt      = Ptop + rhoref.*g0.*ZZ;
    if Nx<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))); end
    
    T    =  (Tp+273.15).*exp(aT./rhoref./cP.*Pt);

    res  = (norm(x(:)-xi(:),2) + norm(f(:)-fi(:),2))./sqrt(2*length(x(:)));
end
rhoBF   = (rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))/2 - rhoref;
rhoo    = rho;
Pto     = Pt;

% get geochemical phase compositions
itm  = it ./(m + x.*KIT ); itx  = it ./(m./KIT  + x);
ctm  = ct ./(m + x.*KCT ); ctx  = ct ./(m./KCT  + x);
ripm = rip./(m + x.*KRIP); ripx = rip./(m./KRIP + x);
ridm = rid./(m + x.*KRID); ridx = rid./(m./KRID + x);
  
% get bulk enthalpy, silica, volatile content densities
S = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx + f.*Dsf - aT./rhoref.*(Pt-Ptop));  
C = rho.*(m.*cm + x.*cx);
V = rho.*(m.*vm + f.*vf);
X = rho.*x;
F = rho.*f;

% get phase entropies
sm = S./rho - x.*Dsx - f.*Dsf;
sx = sm + Dsx;
sf = sm + Dsf;

% get geochemical content densities
IT  = rho.*(m.*itm + x.*itx);
CT  = rho.*(m.*ctm + x.*ctx);
SI  = rho.*si;
RIP = rho.*(m.*ripm + x.*ripx);
RID = rho.*(m.*ridm + x.*ridx);

% initialise reaction/decay rates
Gx = 0.*x;  Gf = 0.*f;
dcy_rip = rho.*rip./HLRIP.*log(2);
dcy_rid = rho.*rid./HLRID.*log(2);

% initialise auxiliary variables 
dSdt   = 0.*T;  diff_T = 0.*T;  diss_h = 0.*T;
dCdt   = 0.*c;  diff_c = 0.*c;
dVdt   = 0.*v;  diff_v = 0.*v;  
dFdt   = 0.*f;  diff_f = 0.*f;  
dXdt   = 0.*x;  diff_x = 0.*x;  
dSIdt  = 0.*si;  diff_si  = 0.*SI;
dITdt  = 0.*IT;  diff_it  = 0.*IT;
dCTdt  = 0.*CT;  diff_ct  = 0.*CT;
dRIPdt = 0.*RIP; diff_rip  = 0.*RIP;
dRIDdt = 0.*RID; diff_rid  = 0.*RID;

% initialise timing and iterative parameters
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumMdt = 0;
dsumSdt = 0;
dsumCdt = 0;
dsumVdt = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','S','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dSdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dFdt','dXdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
        
        xq = x; fq = f;
        SOL = [W(:);U(:);P(:)];
        dcy_rip = rho.*rip./HLRIP.*log(2);
        dcy_rid = rho.*rid./HLRID.*log(2);
        Pto = Pt; etao = eta; rhoo = rho; Div_rhoVo = Div_rhoV;
        update; output;
        time  = time+dt;
        step  = step+1;
    else % continuation file does not exist, start from scratch
        % update coefficients and run initial fluidmech solve
        if ~bnchm
            update;
            fluidmech;
            history;
            output;
        end
        time  = time+dt;
        step  = step+1;
    end
else
    % update coefficients and run initial fluidmech solve
    if ~bnchm
        update;
        fluidmech;
        history;
        output;
    end
    time  = time+dt;
    step  = step+1;
end