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
calibrt =  0;                % not in calibrate mode
TINY    =  1e-16;            % minimum cutoff phase, component fractions
BCA     =  {'',''};          % boundary condition on advection (top/bot, sides)
bnchm   =  0;                % not a benchmark run

% calculate dimensionless numbers characterising the system dynamics
res = 1; x0 = 0; f0 = 0;
while res>1e-16
    ci = c0*(1-f0);
    [x0,cx0,cm0,f0,vf0,vm0] = equilibrium(x0,f0,T0,c0*(1-f0),v0,Ptop,cal,TINY);
    m0 = 1-x0-f0;
    res = abs(c0*(1-f0)-ci)/ci;
end

fprintf('    initial T: %4.3f \n'  ,T0);
fprintf('    initial c: %4.3f \n'  ,c0);
fprintf('    initial v: %4.3f \n'  ,v0);
fprintf('    initial x: %4.3f \n'  ,x0);
fprintf('    initial f: %4.3f \n\n',f0);

% update oxide compositions
wt0 = (cal.perCm-cm0)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm0)./(cal.cphs1-cal.perCm);
cm0_oxd = (wt0(:) .* cal.cmp_oxd(1,:) + (1-wt0(:)) .* cal.cmp_oxd(3,:)) .* (cm0< cal.perCm) ...
        + (wt1(:) .* cal.cmp_oxd(3,:) + (1-wt1(:)) .* cal.cmp_oxd(4,:)) .* (cm0>=cal.perCm);

wt0 = (cal.perCx-cx0)./(cal.perCx-cal.cphs0);
wt1 = (cal.cphs1-cx0)./(cal.cphs1-cal.perCx);
cx0_oxd = (wt0(:) .* cal.cmp_oxd(1,:) + (1-wt0(:)) .* cal.cmp_oxd(2,:)) .* (cx0< cal.perCx) ...
        + (wt1(:) .* cal.cmp_oxd(2,:) + (1-wt1(:)) .* cal.cmp_oxd(4,:)) .* (cx0>=cal.perCx);

c0_oxd = (m0.*cm0_oxd + x0.*cx0_oxd)./(1-f0);

cm0_cmp = cm0_oxd/cal.oxd*100;
cx0_cmp = cx0_oxd/cal.oxd*100;
 c0_cmp =  c0_oxd/cal.oxd*100;

rhof0 = cal.rhof0;
rhom0 = sum(cm0_cmp/100./cal.rhom0).^-1;
rhox0 = sum(cx0_cmp/100./cal.rhox0).^-1;

cm1_oxd = (0.999.*cm0_cmp + 0.001.*cal.cmp(1))*cal.oxd./100;
cm2_oxd = (0.999.*cm0_cmp + 0.001.*cal.cmp(4))*cal.oxd./100;

wtm = [];
wtm([1 2 3 4 6 7 8 9 11 12]) = [cm0_oxd,100.*vm0,0];
etam0 = grdmodel08(wtm,T0);

DrhoT = cal.aT*max([abs(T0-Twall),abs(T0-T1),T0/100]);
Drhoc = abs(sum(cm1_oxd/100./cal.rhom0).^-1-sum(cm2_oxd/100./cal.rhom0).^-1);
Drhox = 0.01*abs(rhox0-rhom0);
Drhof = 0.01*abs(cal.rhof0-rhom0) * (max([v0,v1,vwall])>TINY);
Drho0 = DrhoT + Drhoc + Drhox + Drhof;

uT    = DrhoT*g0*(D/10)^2/etam0/etareg;
uc    = Drhoc*g0*(D/10)^2/etam0/etareg;
ux    = Drhox*g0*(D/10)^2/etam0/etareg;
uf    = Drhof*g0*(D/10)^2/etam0/etareg * (max([v0,v1,vwall])>TINY);
u0    = Drho0*g0*(D/10)^2/etam0/etareg;

wx0   = abs(rhox0-rhom0)*g0*dx^2/etam0;
wf0   = abs(rhof0-rhom0)*g0*df^2/etam0 * (max([v0,v1,vwall])>TINY);

ud0   = kT0/rhom0/cP/(D/10);

uwT   = bnd_w/tau_T; 
uwc   = bnd_w/tau_a; 

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
smth = smth*Nx*Nz*1e-4;
rp   = randn(Nz,Nx);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 + diff(rp(2:end-1,:),2,2)./8;
    rp = rp - mean(mean(rp(2:end-1,2:end-1)));
    rp([1 end],:) = 0;
    rp(:,[1 end]) = 0;
end
rp = rp./max(abs(rp(:)));

gp = exp(-(XX-D/2).^2/(D/8)^2 - (ZZ-D/2).^2/(D/8)^2);

% get mapping arrays
NP =  Nz   * Nx   ;
NW = (Nz-1)* Nx   ;
NU =  Nz   *(Nx-1);
MapP = reshape(1:NP,Nz  ,Nx  );
MapW = reshape(1:NW,Nz-1,Nx  );
MapU = reshape(1:NU,Nz  ,Nx-1) + NW;

if bnd_h>0
    switch bndmode
        case 0  % none
            bndinit = zeros(size(ZZ));
        case 1  % top only
            bndinit = (1+erf( ( -ZZ+bnd_h)/bnd_w))/2;
        case 2  % bot only
            bndinit = (1+erf(-(D-ZZ-bnd_h)/bnd_w))/2;
        case 3  % top/bot only
            bndinit = (1+erf( ( -ZZ+bnd_h)/bnd_w))/2 ...
                + (1+erf(-(D-ZZ-bnd_h)/bnd_w))/2;
        case 4 % all walls
            bndinit = (1+erf( ( -ZZ+bnd_h)/bnd_w))/2 ...
                + (1+erf(-(D-ZZ-bnd_h)/bnd_w))/2 ...
                + (1+erf( ( -XX+bnd_h)/bnd_w))/2 ...
                + (1+erf(-(L-XX-bnd_h)/bnd_w))/2;
                + (1+erf(-(D-ZZ-bnd_h)/bnd_w))/2;
        case 5 % only walls
            bndinit = (1+erf( ( -XX+bnd_h)/bnd_w))/2 ...
                    + (1+erf(-(L-XX-bnd_h)/bnd_w))/2;
    end
    bndinit = max(0,min(1,bndinit));
else
    bndinit = zeros(size(ZZ));
end

switch bndmode
    case 0  % none
        bndshape = zeros(size(ZZ(inz,inx)));
    case 1  % top only
        bndshape = exp( ( -ZZ(inz,inx))/bnd_w);
    case 2  % bot only
        bndshape = exp(-(D-ZZ(inz,inx))/bnd_w);
    case 3  % top/bot only
        bndshape = exp( ( -ZZ(inz,inx))/bnd_w) ...
                 + exp(-(D-ZZ(inz,inx))/bnd_w);
    case 4 % all walls
        bndshape = exp( ( -ZZ(inz,inx))/bnd_w) ...
                 + exp(-(D-ZZ(inz,inx))/bnd_w) ...
                 + exp( ( -XX(inz,inx))/bnd_w) ...
                 + exp(-(L-XX(inz,inx))/bnd_w);
    case 5 % only walls
        bndshape = exp( ( -XX(inz,inx))/bnd_w) ...
                 + exp(-(L-XX(inz,inx))/bnd_w);
end
bndshape = max(0,min(1,bndshape));

bnd_S = zeros(size(bndshape));
bnd_C = zeros(size(bndshape));
bnd_V = zeros(size(bndshape));

% set specified boundaries to no slip, else to free slip
if bndmode>=4;               sds = +1;      % no slip sides for 'all sides(4)'
else;                        sds = -1; end  % free slip sides for other types
if bndmode==1 || bndmode>=3; top = +1;      % no slip top for 'top only(1)', 'top/bot(3)', 'all sides(4)'
else;                        top = -1; end  % free slip for other types
if bndmode>=2;               bot = +1;      % no slip bot for 'bot only(2)', 'top/bot(3)', 'all sides(4)'
else;                        bot = -1; end  % free slip for other types
if bndmode==5;               top = -1; bot = -1; end % free slip top/bot for 'only walls(5)'

% initialise solution fields
Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay)/wlay_T))/2 + dTr.*rp + dTg.*gp;  if any(bndinit(:)) && ~isnan(Twall); Tp = Tp + (Twall-Tp).*bndinit; end % potential temperature [C]
c   =  c0 + (c1-c0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dcr.*rp + dcg.*gp;  if any(bndinit(:)) && ~isnan(cwall); c  = c  + (cwall-c ).*bndinit; end; cin = c; % major component
v   =  v0 + (v1-v0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dvr.*rp + dvg.*gp;  if any(bndinit(:)) && ~isnan(vwall); v  = v  + (vwall-v ).*bndinit; end; vin = v; % volatile component

te = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    te(:,:,i)  =  te0(i) + (te1(i)-te0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dter(i).*rp + dteg(i).*gp;  % trace elements
    if any(bndinit(:)) && ~isnan(tewall(i)); te(:,:,i)  = te(:,:,i) + (tewall(i)-te(:,:,i)).*bndinit; end; tein = te; 
end
ir = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir
    ir(:,:,i)  =  ir0(i) + (ir1(i)-ir0(i)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dirr(i).*rp + dirg(i).*gp;  % isotope ratios  
    if any(bndinit(:)) && ~isnan(irwall(i)); ir(:,:,i)  = ir(:,:,i) + (irwall(i)-ir(:,:,i)).*bndinit; end; irin = ir;
end

U   =  zeros(size((XX(:,1:end-1)+XX(:,2:end))));  Ui = U;  res_U = 0.*U;
W   =  zeros(size((XX(1:end-1,:)+XX(2:end,:))));  Wi = W;  res_W = 0.*W; wf = 0.*W; wc = 0.*W;
P   =  0.*Tp;  Pi = P;  res_P = 0.*P;  meanQ = 0;  Vel = 0.*P;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wf  = W;  Uf  = U; 
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_Vo = 0;  Div_rhoV = 0.*P(inz,inx);  Div_rhoVo = Div_rhoV;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(Nz-1,Nx-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(Nz-1,Nx-1);  tII = 0.*P; 
VolSrc =  0.*P(inz,inx);  MassErr = 0;  drhodt = 0.*P;  drhodto = 0.*P;
rho    =  rhom0.*ones(size(Tp));
rhoref =  mean(rho(inz,inx),'all');
Pt     =  Ptop + rhoref.*g0.*ZZ .* 1.01;

% get volume fractions and bulk density
step    = 0;
theta   = 1/2;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
res  = 1;  tol = 1e-13;  x = zeros(size(Tp));  f = v/2;
while res > tol
    Pti = Pt; xi = x; fi = f;
    
    rhoref =  mean(rho(inz,inx),'all');
    rhofz  = (rho(1:end-1,:)+rho(2:end,:))/2;
    Adbt   =  cal.aT./rhoref;
    Pt( 2:end, :) = repmat(cumsum(mean(rhofz,2).*g0.*h),1,Nx);
    Pt            = Pt - Pt(2,:)/2 + Ptop;
    Pt([1 end],:) = Pt([2 end-1],:);
    Pt(:,[1 end]) = Pt(:,[2 end-1]);
    if Nz<=10; Pt = Ptop.*ones(size(Tp)); end

    T    =  (Tp+273.15).*exp(Adbt./cP.*Pt);

    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T-273.15,c,v,Pt,cal,TINY);
    
    v  = vin .*(1-x);
    c  = cin .*(1-f);
    te = tein.*(1-f);
    ir = irin.*(1-f);

    x  = xq;  f = fq;  m = 1-x-f;
    cm = cmq; cx = cxq;
    vm = vmq; vf = vfq;
    Kc = cxq./cmq;
    Kf = vfq./vmq;

    update;

    res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2) ...
         + norm(x(:)-xi(:),2)./(norm(x(:),2)+TINY) ...
         + norm(f(:)-fi(:),2)./(norm(f(:),2)+TINY);
end
rhoo = rho;
dto  = dt; 

% get bulk enthalpy, silica, volatile content densities
S  = rho.*(cP.*log(T/(cal.Tphs1+273.15)) + x.*Dsx + f.*Dsf - Adbt.*(Pt-Ptop));  
S0 = rho.*(cP.*log(cal.Tphs1+273.15) + x.*Dsx + f.*Dsf - Adbt.*Ptop);  
C  = rho.*(m.*cm + x.*cx);
V  = rho.*(m.*vm + f.*vf);
X  = rho.*x;
F  = rho.*f;
M  = rho.*m;

% get phase entropies
sm = S./rho - x.*Dsx - f.*Dsf;
sx = sm + Dsx;
sf = sm + Dsf;

% get trace element phase compositions
Kte = zeros(Nz,Nx,cal.nte);
tem = zeros(Nz,Nx,cal.nte);
tex = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    for j=1:cal.nc; Kte(:,:,i) = Kte(:,:,i) + cal.Kte_cmp(i,j) .* c_cmp(:,:,j)./100; end

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

IR = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir
    IR(:,:,i)  = rho.*(m.*irm(:,:,i) + x.*irx(:,:,i));
end

% initialise phase change rates
Gx = 0.*x(inz,inx);  Gf = 0.*f(inz,inx);  Gm = 0.*m(inz,inx);

% initialise auxiliary variables 
dSdt   = 0.*T(inz,inx);  diss_h = 0.*T(inz,inx);
dCdt   = 0.*c(inz,inx);
dVdt   = 0.*v(inz,inx);
dFdt   = 0.*f(inz,inx);
dMdt   = 0.*m(inz,inx);
dXdt   = 0.*x(inz,inx);
dTEdt  = 0.*te(inz,inx);
dIRdt  = 0.*ir(inz,inx);

% initialise timing and iterative parameters
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumMdt = 0;
dsumSdt = 0;
dsumCdt = 0;
dsumVdt = 0;
dsumCdt_oxd = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','S','C','V','T','c','v','cm','cx','vm','vf','TE','IR','te','ir','dSdt','dCdt','dVdt','dFdt','dXdt','dTEdt','dIRdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm');
        name = [opdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        xq = x; fq = f;
        SOL = [W(:);U(:);P(:)];
        rhoo = rho; Div_rhoVo = Div_rhoV;
        update;
    else % continuation file does not exist, start from scratch
        fprintf('\n   !!! restart file does not exist !!! \n   => starting run from scratch %s \n\n',runID);
        update;
        fluidmech;
        history;
        output;
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
