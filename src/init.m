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

T = T0+273.15;
for i=1:5

% calculate dimensionless numbers characterising the system dynamics
[x0,cx0,cm0,f0,vf0,vm0] = equilibrium(0.01,v0/10,T-273.15,c0,v0,Ptop,cal,TINY);
m0 = 1-x0-f0;

% update oxide compositions
wt0 = (cal.perCm-cm0)./(cal.perCm-cal.cphs0);
wt1 = (cal.cphs1-cm0)./(cal.cphs1-cal.perCm);
cm0_cmp = (wt0(:) .* [1 0 0 0] + (1-wt0(:)) .* [0 0 1 0]) .* (cm0< cal.perCm) ...
        + (wt1(:) .* [0 0 1 0] + (1-wt1(:)) .* [0 0 0 1]) .* (cm0>=cal.perCm);

wt0 = (cal.perCx-cx0)./(cal.perCx-cal.cphs0);
wt1 = (cal.cphs1-cx0)./(cal.cphs1-cal.perCx);
cx0_cmp = (wt0(:) .* [1 0 0 0] + (1-wt0(:)) .* [0 1 0 0]) .* (cx0< cal.perCx) ...
        + (wt1(:) .* [0 1 0 0] + (1-wt1(:)) .* [0 0 0 1]) .* (cx0>=cal.perCx);

c0_cmp = (m0.*cm0_cmp + x0.*cx0_cmp)./(1-f0);

cm0_mem = cm0_cmp*cal.cmp_mem/100;
cx0_mem = cx0_cmp*cal.cmp_mem/100;
 c0_mem =  c0_cmp*cal.cmp_mem/100;

cm0_oxd = cm0_mem*cal.mem_oxd/100;
cx0_oxd = cx0_mem*cal.mem_oxd/100;
 c0_oxd =  c0_mem*cal.mem_oxd/100;

rhof0 = cal.rhof0;
rhox0 = sum(cx0_mem./cal.rhox0).^-1;

cm1_oxd = (0.9.*cm0_mem + 0.1.*cal.cmp_mem(1,:)/100)*cal.mem_oxd/100;
cm2_oxd = (0.9.*cm0_mem + 0.1.*cal.cmp_mem(4,:)/100)*cal.mem_oxd/100;

wtm   = [cm0_oxd.*100,100.*vm0];
etam0 = Giordano08(wtm,T0);
rhom0 = DensityX(wtm,T0,Ptop/1e8);

rho0  = (x0./rhox0 + m0./rhom0).^-1;
T     = (T0+273.15).*exp(cal.aT./rho0./cP.*Ptop);

end

fprintf('    initial T: %4.3f \n'  ,T0);
fprintf('    initial c: %4.3f \n'  ,c0);
fprintf('    initial v: %4.3f \n'  ,v0);
fprintf('    initial x: %4.3f \n'  ,x0);
fprintf('    initial f: %4.3f \n\n',f0);

wtm   = [cm1_oxd.*100,100.*vm0];
rhom1 = DensityX(wtm,T0,Ptop/1e8);

wtm   = [cm2_oxd.*100,100.*vm0];
rhom2 = DensityX(wtm,T0,Ptop/1e8);

DrhoT = rhom0.*cal.aT*max([abs(T0-Twall)/10,abs(T0-T1),T0/100]);
Drhoc = abs(rhom1-rhom2);
Drhox = 0.1*abs(rhox0-rhom0);
Drhof = 0.1*abs(cal.rhof0-rhom0) * (max([v0,v1,vwall])>TINY);
Drho0 = DrhoT + Drhoc + Drhox + Drhof;

uT    = DrhoT*g0*(D/10)^2/etam0/cnvreg;
uc    = Drhoc*g0*(D/10)^2/etam0/cnvreg;
ux    = Drhox*g0*(D/10)^2/etam0/cnvreg;
uf    = Drhof*g0*(D/10)^2/etam0/cnvreg * (max([v0,v1,vwall])>TINY);
u0    = Drho0*g0*(D/10)^2/etam0/cnvreg;

wx0   = abs(rhox0-rhom0)*g0*dx^2/etam0/sgrreg;
wf0   = abs(rhof0-rhom0)*g0*df^2/etam0/sgrreg * (max([v0,v1,vwall])>TINY);

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

Re    = u0*rhom0*(D/10)/etam0/cnvreg;
Rex   = wx0*rhom0*dx/etam0/sgrreg;
Ref   = wf0*rhom0*df/etam0/sgrreg;

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

bnd_S = zeros(size(topshape));
bnd_C = zeros(size(topshape));
bnd_V = zeros(size(topshape));

% set specified boundaries to no slip, else to free slip
if bndmode>=4;               sds = +1;      % no slip sides for 'all sides(4)'
else;                        sds = -1; end  % free slip sides for other types
if bndmode==1 || bndmode>=3; top = +1;      % no slip top for 'top only(1)', 'top/bot(3)', 'all sides(4)'
else;                        top = -1; end  % free slip for other types
if bndmode>=2;               bot = +1;      % no slip bot for 'bot only(2)', 'top/bot(3)', 'all sides(4)'
else;                        bot = -1; end  % free slip for other types
if bndmode==5;               top = -1; bot = -1; end % free slip top/bot for 'only walls(5)'

% initialise solution fields
Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_T))/2 + dTr.*rp + dTg.*gp;  % potential temperature [C]
c   =  c0 + (c1-c0) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_c))/2 + dcr.*rp + dcg.*gp;  % major component
v   =  v0 + (v1-v0) .* (1+erf((ZZ/D-zlay+rp*h*dlay)/wlay_c))/2 + dvr.*rp + dvg.*gp;  % volatile component

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

if any(topinit(:)) && ~isnan(cwall(1)); c  = c  + (cwall(1)-c ).*topinit; end
if any(botinit(:)) && ~isnan(cwall(2)); c  = c  + (cwall(2)-c ).*botinit; end
if any(sdsinit(:)) && ~isnan(cwall(3)); c  = c  + (cwall(3)-c ).*sdsinit; end
cin = c;

if any(topinit(:)) && ~isnan(vwall(1)); v  = v  + (vwall(1)-v ).*topinit; end
if any(botinit(:)) && ~isnan(vwall(2)); v  = v  + (vwall(2)-v ).*botinit; end
if any(sdsinit(:)) && ~isnan(vwall(3)); v  = v  + (vwall(3)-v ).*sdsinit; end
vin = v;

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
P   =  zeros(Nz+2,Nx+2);  Vel = 0.*Tp;
SOL = [W(:);U(:);P(:)];

% initialise auxiliary fields
Wf  = W;  Uf  = U; 
Wx  = W;  Ux  = U;
Wm  = W;  Um  = U;

eIIref = 1e-6;  
Div_V  = 0.*c;  advn_rho = 0.*c;  drhodt = 0.*c;  drhodto = drhodt;
exx    = 0.*c;  ezz = 0.*c;  exz = zeros(Nz-1,Nx-1);  eII = 0.*c;  
txx    = 0.*c;  tzz = 0.*c;  txz = zeros(Nz-1,Nx-1);  tII = 0.*c; 
VolSrc = 0.*c; 
rhom   = rhom0.*ones(size(c)); 
rhox   = rhox0.*ones(size(c));
rhof   = rhof0.*ones(size(c));
rho    = rhom0.*ones(size(c));
rhoref = mean(rho,'all');
Pt     = Ptop + rhoref.*g0.*ZZ;
xq     = ones(size(c))/10;  fq = v/2;
cmq    = c; cxq = c; 
vmq    = v; vfq = ones(size(v)); 

% get volume fractions and bulk density
step    = 0;
EQtime  = 0;
FMtime  = 0;
TCtime  = 0;
UDtime  = 0;
res  = 1;  tol = 1e-12;  
while res > tol
    Pti = Pt; xi = xq; fi = fq;
    
    rhoref = mean(rho,'all');
    Adbt   = cal.aT./rhoref;
    if Nz==1; Pt = Ptop.*ones(size(Tp)); else
        rhofz  = (rho(1:end-1,:)+rho(2:end,:))/2;
        Pt(1,:)     = repmat(mean(rhofz(1,:),2).*g0.*h/2,1,Nx) + Ptop;
        Pt(2:end,:) = Pt(1,:) + repmat(cumsum(mean(rhofz,2).*g0.*h),1,Nx);
    end

    T      =  (Tp+273.15).*exp(Adbt./cP.*Pt);

    eqtime = tic;

    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(xq,fq,T-273.15,c,v,Pt,cal,TINY);

    mq = 1-xq-fq;

    x  = xq;  f = fq;  m = mq;
    cm = cmq; cx = cxq;
    vm = vmq; vf = vfq;

    eqtime = toc(eqtime);
    EQtime = EQtime + eqtime;

    update;

    res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2) ...
         + norm((x(:)-xi(:)).*(x(:)>TINY^0.5),2)./(norm(x(:),2)+TINY) ...
         + norm((f(:)-fi(:)).*(f(:)>TINY^0.5),2)./(norm(f(:),2)+TINY);
end
rhoo = rho;
dto  = dt; 
ho   = h;

% get bulk enthalpy, silica, volatile content densities
S   = rho.*(cP.*log(T/(cal.Tphs1+273.15)) + x.*Dsx + f.*Dsf - Adbt.*(Pt-Ptop));  So = S;  res_S = 0.*S;
S0  = rho.*(cP.*log(cal.Tphs1+273.15) + x.*Dsx + f.*Dsf - Adbt.*Ptop);  
C   = rho.*(m.*cm + x.*cx); Co = C;  res_C = 0.*C;
V   = rho.*(m.*vm + f.*vf); Vo = V;  res_V = 0.*V;
X   = rho.*x; Xo = X;  res_X = 0.*X;
F   = rho.*f; Fo = F;  res_F = 0.*F;
M   = rho.*m; Mo = M;  res_M = 0.*M;
RHO = X+M+F;

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
dVdt   = 0.*v;  dVdto  = dVdt;
dFdt   = 0.*f;  dFdto  = dFdt;
dXdt   = 0.*x;  dXdto  = dXdt;
dMdt   = 0.*m;  dMdto  = dMdt;
dTEdt  = 0.*te; dTEdto = dTEdt;
dIRdt  = 0.*ir; dIRdto = dIRdt;

% initialise timing and iterative parameters
step    = 0;
time    = 0;
iter    = 0;
hist    = [];
dsumMdt = 0; dsumMdto = 0;
dsumSdt = 0; dsumSdto = 0;
dsumCdt = 0; dsumCdto = 0;
dsumVdt = 0; dsumVdto = 0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','M','S','C','V','T','c','v','cm','cx','vm','vf','TE','IR','te','ir','dSdt','dCdt','dVdt','dFdt','dXdt','dMdt','drhodt','dTEdt','dIRdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm');
        name = [opdir,'/',runID,'/',runID,'_hist'];
        load(name,'hist');

        xq  = x; fq = f; mq = 1-x-f;
        SOL = [W(:);U(:);P(:)];

        So = S;
        Co = C;
        Vo = V;
        Xo = X;
        Fo = F;
        Mo = M;
        rhoo = rho;
        TEo = TE;
        IRo = IR;
        dSdto = dSdt;
        dCdto = dCdt;
        dVdto = dVdt;
        dXdto = dXdt;
        dFdto = dFdt;
        dMdto = dMdt;
        drhodto = drhodt;
        dTEdto = dTEdt;
        dIRdto = dIRdt;

        update; output; restart = 0;
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
