% create output directory
if ~isfolder([opdir,'/',runID])
    mkdir([opdir,'/',runID]);
end

% save input parameters and runtime options (unless restarting)
if restart == 0 
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
calibrt =  0;                % not in calibrate mode
TINY    =  1e-16;            % minimum cutoff phase, component fractions

% calculate dimensionless numbers characterising the system dynamics
[x0,cx0,cm0,f0,vf0,vm0] = equilibrium(0,0,T0,c0,v0,Ptop,cal,TINY);

wt0 = (cal.perCm-cm0)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm0)./(cal.cphs1-cal.perCm);
cm0_oxds = (wt0(:) .* cal.oxds(1,:) + (1-wt0(:)) .* cal.oxds(2,:)) .* (cm0(:)<=cal.perCm) ...
         + (wt1(:) .* cal.oxds(2,:) + (1-wt1(:)) .* cal.oxds(3,:)) .* (cm0(:)> cal.perCm);

wtm([1 3 4 6 7 8 9 11 12]) = [cm0_oxds,100.*vm0,0]; % SiO2
etam0     = grdmodel08(wtm,T0);

DrhoT = aT*max([abs(T0-Twall),abs(T0-T1),T0/100]);
Drhoc = gC*max([abs(c0-cwall),abs(c0-c1),c0/100]); 
Drhox = x0/100*abs(rhox0-rhom0);
Drhof = f0/100*abs(rhof0-rhom0);
Drho0 = DrhoT + Drhoc + Drhox + Drhof;

uT    = DrhoT*g0*(D/10)^2/etam0/etareg;
uc    = Drhoc*g0*(D/10)^2/etam0/etareg;
ux    = Drhox*g0*(D/10)^2/etam0/etareg;
uf    = Drhof*g0*(D/10)^2/etam0/etareg .* (max([v0,v1,vwall])>TINY);
u0    = Drho0*g0*(D/10)^2/etam0/etareg;

wx0   = abs(rhox0-rhom0)*g0*dx^2/etam0;
wf0   = abs(rhof0-rhom0)*g0*df^2/etam0 .* (max([v0,v1,vwall])>TINY);

ud0   = kT/rhom0/cP/(D/10);

uwT   = dw/tau_T; 
uwc   = dw/tau_a; 

RaT   = uT/ud0;
Rac   = uc/ud0;
Rax   = ux/ud0;
Raf   = uf/ud0;
Ra    = u0/ud0;

Rux   = wx0/u0;
Ruf   = wf0/u0;

RwT   = uwT/u0;
Rwc   = uwc/u0;

Re    = u0*rhom0*(D/10)/etam0/etareg;
Rex   = wx0*rhom0*dx/etam0;
Ref   = wf0*rhom0*df/etam0;

fprintf('    crystal Re: %1.3e \n'  ,Rex);
fprintf('     bubble Re: %1.3e \n'  ,Ref);
fprintf('     system Re: %1.3e \n\n',Re );

fprintf('    thermal Ra: %1.3e \n'  ,RaT);
fprintf('   chemical Ra: %1.3e \n'  ,Rac);
fprintf('    crystal Ra: %1.3e \n'  ,Rax);
fprintf('     bubble Ra: %1.3e \n'  ,Raf);
fprintf('   combined Ra: %1.3e \n\n',Ra );

fprintf('    crystal Ru: %1.3e \n'  ,Rux);
fprintf('     bubble Ru: %1.3e \n\n',Ruf);

fprintf('    thermal Rw: %1.3e \n'  ,RwT);
fprintf('   chemical Rw: %1.3e \n\n',Rwc);


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

inz = 2:Nz-1;
inx = 2:Nx-1;

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

bndS = 0.*bndshape;  bndC =  0.*bndshape;  bndV =  0.*bndshape;  bndSI =  0.*bndshape;

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
Div_V  =  0.*P;  Div_rhoV = 0.*P(inz,inx);  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(Nz-1,Nx-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(Nz-1,Nx-1);  tII = 0.*P; 
VolSrc =  0.*P;  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;

rhoo =  rhom0.*ones(size(Tp)); rhoref = rhom0;  %#ok<NASGU>
dto  =  dt;
Pt   =  rhoref.*g0.*ZZ + Ptop;  
if Nx<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))).*ones(size(Pt)); end
T    =  (Tp+273.15).*exp(aT./rhoref./cP.*Pt); % real temperature [K]

% get volume fractions and bulk density
step    = 0;
theta   = 1/2;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
res  = 1;  tol = 1e-12;  x = ones(size(T))./10;  f = v/2;
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
    if Nz<=10; Pt = mean(mean(Pt(2:end-1,2:end-1))); end
    
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
dSdt   = 0.*T(inz,inx);  dTdt   = 0.*T(inz,inx);  diss_h = 0.*T(inz,inx);
dCdt   = 0.*c(inz,inx);
dVdt   = 0.*v(inz,inx);
dFdt   = 0.*f(inz,inx);
dXdt   = 0.*x(inz,inx);
dSIdt  = 0.*si(inz,inx);
dITdt  = 0.*IT(inz,inx);
dCTdt  = 0.*CT(inz,inx);
dRIPdt = 0.*RIP(inz,inx);
dRIDdt = 0.*RID(inz,inx);

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
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','S','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dSdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dFdt','dXdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
        
        xq = x; fq = f;
        SOL = [W(:);U(:);P(:)];
        dcy_rip = rho.*rip./HLRIP.*log(2);
        dcy_rid = rho.*rid./HLRID.*log(2);
        Pto = Pt; etao = eta; rhoo = rho; Div_rhoVo = Div_rhoV;
        update; output;
    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',name);
    end
else
    % complete, plot, and save initial condition
    update;
    fluidmech;
    history;
    output;
end

time  = time+dt;
step  = step+1;
