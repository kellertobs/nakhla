load ocean;  % load custom colormap

% get coordinate arrays
X         = -h/2:h:L+h/2;
Z         = -h/2:h:D+h/2;
[XX,ZZ]   = meshgrid(X,Z);
Xfc       = (X(1:end-1)+X(2:end))./2;
Zfc       = (Z(1:end-1)+Z(2:end))./2;
[XXu,ZZu] = meshgrid(Xfc,Z);
[XXw,ZZw] = meshgrid(X,Zfc);

Nx = length(X);
Nz = length(Z);

% get smoothed initialisation field
rng(seed);
rp = randn(Nz,Nx);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 + diff(rp(2:end-1,:),2,2)./8;
    rp([1 2 end-1 end],:) = 0;
    rp(:,[1 2 end-1 end]) = 0;
end
rp = rp./max(abs(rp(:)));
rp = rp - mean(mean(rp(2:end-1,2:end-1)));

% get mapping arrays
NP =  Nz   * Nx   ;
NW = (Nz-1)* Nx   ;
NU =  Nz   *(Nx-1);
MapP = reshape(1:NP,Nz  ,Nx  );
MapW = reshape(1:NW,Nz-1,Nx  );
MapU = reshape(1:NU,Nz  ,Nx-1) + NW;

% boundary conditions shape function
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

% set initial solution fields
T   =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay)/wlay_T))/2 + dT.*rp;  if bndinit && ~isnan(Twall); T = T + (Twall-T).*bndshape; end % temperature
c   =  c0 + (c1-c0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dc.*rp;  if bndinit && ~isnan(cwall); c = c + (cwall-c).*bndshape; end % major component
v   =  v0 + (v1-v0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dv.*rp;  if bndinit && ~isnan(vwall); v = v + (vwall-v).*bndshape; end % volatile component

it  =  it0 + (it1-it0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dit.*rp;  if bndinit && ~isnan(itwall); it  = it  + (itwall-it ).*bndshape; end % incompatible trace element
ct  =  ct0 + (ct1-ct0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dct.*rp;  if bndinit && ~isnan(ctwall); ct  = ct  + (ctwall-ct ).*bndshape; end % compactible trace element
si  =  si0 + (si1-si0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dsi.*rp;  if bndinit && ~isnan(siwall); si  = si  + (siwall-si ).*bndshape; end % stable isotope ratio
rip =  ri0 + (ri1-ri0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dri.*rp;  if bndinit && ~isnan(riwall); rip = rip + (riwall-rip).*bndshape; end % radiogenic isotope parent
rid =  rip.*HLRID./HLRIP;                                           % radiogenic isotope daughter

U   =  zeros(size((XX(:,1:end-1)+XX(:,2:end))));  Ui = U;  res_U = 0.*U;
W   =  zeros(size((XX(1:end-1,:)+XX(2:end,:))));  Wi = W;  res_W = 0.*W; wf = 0.*W; wc = 0.*W;
P   =  0.*c;  Pi = P;  res_P = 0.*P;  meanQ = 0;  Pt = rhom0.*g0.*ZZ + Ptop;
S   = [W(:);U(:);P(:)];

% get volume fractions and bulk density
res = 1;  tol = 1e-15;  iter = 0;  x = 0;  f = 0;
while res > tol
    xi = x;  fi = f;
    
    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
    mq = 1-xq-fq;
    
    x  = xq;  f = fq;  m = mq;
    cm = cmq; cx = cxq;
    vm = vmq; vf = vfq;
        
    rhom = rhom0 .* (1 - aTm.*(T-Tphs0) - gCm.*(cm-cphs0));
    rhox = rhox0 .* (1 - aTx.*(T-Tphs0) - gCx.*(cx-cphs0));
    rhof = rhof0 .* (1 - aTf.*(T-Tphs0) + bPf.*(Pt-Ptop));
    
    rho  = 1./(m./rhom + x./rhox + f./rhof);

    chi  = x.*rho./rhox;
    phi  = f.*rho./rhof;
    mu   = m.*rho./rhom;
    
    rhoref  = mean(mean(rho(2:end-1,2:end-1)));
    Pt      = Ptop + rhoref.*g0.*ZZ;
    if Nx<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))); end
    
    res  = (norm(x(:)-xi(:),2) + norm(f(:)-fi(:),2))./sqrt(2*length(x(:)));
    iter = iter+1;
end
rhoBF   = (rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))/2 - rhoref;
rhoo    = rho;
Pto     = Pt;
  
% get bulk enthalpy, silica, volatile content densities
if ~react;  Dsx = 0;  Dsf = 0;  end
rhoCp = rho.*(m.*Cpm + x.*Cpx + f.*Cpf);
rhoDs = rho.*(m.*0   + x.*Dsx + f.*Dsf);
H = (rhoCp + rhoDs).*T;    sumH0 = sum(sum(H(2:end-1,2:end-1)));
C = rho.*(m.*cm + x.*cx);  sumC0 = sum(sum(C(2:end-1,2:end-1)));
V = rho.*(m.*vm + f.*vf);  sumV0 = sum(sum(V(2:end-1,2:end-1)));

itm  = it./(m + x.*KIT); itx = it./(m./KIT + x);
ctm  = ct./(m + x.*KCT); ctx = ct./(m./KCT + x);
sim  = si;               six = si;
ripm = rip./(m + x.*KRIP); ripx = rip./(m./KRIP + x);
ridm = rid./(m + x.*KRID); ridx = rid./(m./KRID + x);

IT  = rho.*(m.*itm + x.*itx);
CT  = rho.*(m.*ctm + x.*ctx);
SIm = rho.* m.*sim;  
SIx = rho.* x.*six;
SI  = SIm + SIx;
RIP = rho.*(m.*ripm + x.*ripx);
RID = rho.*(m.*ridm + x.*ridx);

% initialise reaction/decay rates
Gx = 0.*x;  Gf = 0.*f;
dcy_rip = rho.*rip./HLRIP.*log(2);
dcy_rid = rho.*rid./HLRID.*log(2);

% initialise auxiliary variables 
dHdt   = 0.*T;  diff_T = 0.*T;
dCdt   = 0.*c;  diff_c = 0.*c;
dVdt   = 0.*v;  diff_v = 0.*v;  
dfdt   = 0.*f;
dxdt   = 0.*x;
dSImdt = 0.*si;  diff_sim = 0.*sim;
dSIxdt = 0.*si;  diff_six = 0.*SIx;
dITdt  = 0.*IT;  diff_it  = 0.*IT;
dCTdt  = 0.*CT;  diff_ct  = 0.*CT;
dRIPdt = 0.*RIP; diff_rip  = 0.*RIP;
dRIDdt = 0.*RID; diff_rid  = 0.*RID;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_rhoV = 0.*P;  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(size(f)-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(size(f)-1);  tII = 0.*P; 
VolSrc =  0.*P;  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;

% initialise timing and iterative parameters
step   =  0;
time   =  0;
iter   =  0;
hist   = [];

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart)];
    end
    load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SIm','SIx','SI','RIP','RID','it','ct','sim','six','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSImdt','dSIxdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist');
%     name = [opdir,runID,'/',runID,'_par'];
%     load(name);
    
    Pscale = geomean(eta(:))/h;
    S = [W(:);U(:);P(:)/Pscale];
    Pto = Pt; etao = eta; rhoo = rho; Div_rhoVo = Div_rhoV;
    update; output;
    time = time+dt;
    step = step+1;
else
    update; output;
end