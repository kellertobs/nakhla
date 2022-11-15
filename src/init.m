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
[x0,cx0,cm0,f0,vf0,vm0] = equilibrium(0.5,v0/2,T0,c0,v0,Ptop,cal,TINY);

fprintf('    initial T: %4.3f \n'  ,T0);
fprintf('    initial c: %4.3f \n'  ,c0);
fprintf('    initial v: %4.3f \n'  ,v0);
fprintf('    initial x: %4.3f \n'  ,x0);
fprintf('    initial f: %4.3f \n\n',f0);

wt0 = (cal.perCx-cm0)./(cal.perCx-cal.cphs0);
wt1 = (cal.perCm-cm0)./(cal.perCm-cal.perCx);
wt2 = (cal.cphs1-cm0)./(cal.cphs1-cal.perCm);
cm0_cmp = (wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(2,:)) .* (cm0(:)< cal.perCx) ...
        + (wt1(:) .* cal.cmp(2,:) + (1-wt1(:)) .* cal.cmp(3,:)) .* (cm0(:)>=cal.perCx & cm0(:)<=cal.perCm) ...
        + (wt2(:) .* cal.cmp(3,:) + (1-wt2(:)) .* cal.cmp(4,:)) .* (cm0(:)> cal.perCm);

cm0_oxd = cm0_cmp*cal.oxd./100;

wt0 = (cal.perCx-cx0)./(cal.perCx-cal.cphs0);
wt1 = (cal.perCm-cx0)./(cal.perCm-cal.perCx);
wt2 = (cal.cphs1-cx0)./(cal.cphs1-cal.perCm);
cx0_cmp = (wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(2,:)) .* (cx0(:)< cal.perCx) ...
        + (wt1(:) .* cal.cmp(2,:) + (1-wt1(:)) .* cal.cmp(3,:)) .* (cx0(:)>=cal.perCx & cx0(:)<=cal.perCm) ...
        + (wt2(:) .* cal.cmp(3,:) + (1-wt2(:)) .* cal.cmp(4,:)) .* (cx0(:)> cal.perCm);

cx0_oxd = cx0_cmp*cal.oxd./100;

wt0 = (cal.perCx-c0)./(cal.perCx-cal.cphs0);
wt1 = (cal.perCm-c0)./(cal.perCm-cal.perCx);
wt2 = (cal.cphs1-c0)./(cal.cphs1-cal.perCm);
c0_cmp = (wt0(:) .* cal.cmp(1,:) + (1-wt0(:)) .* cal.cmp(2,:)) .* (c0(:)< cal.perCx) ...
        + (wt1(:) .* cal.cmp(2,:) + (1-wt1(:)) .* cal.cmp(3,:)) .* (c0(:)>=cal.perCx & c0(:)<=cal.perCm) ...
        + (wt2(:) .* cal.cmp(3,:) + (1-wt2(:)) .* cal.cmp(4,:)) .* (c0(:)> cal.perCm);

c0_oxd = c0_cmp*cal.oxd./100;

rhof0 = cal.rhof0;
rhom0 = sum(cm0_oxd/100./cal.rhom0).^-1;
rhox0 = sum(cx0_oxd/100./cal.rhox0).^-1;

cm1_oxd = (0.999.*cm0_cmp + 0.001.*cal.cmp(1))*cal.oxd./100;
cm2_oxd = (0.999.*cm0_cmp + 0.001.*cal.cmp(4))*cal.oxd./100;

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

wx0   = abs(rhox0-rhom0)*g0*d0^2/etam0;
wf0   = abs(rhof0-rhom0)*g0*d0^2/etam0 * (max([v0,v1,vwall])>TINY);

ud0   = kT0/rhom0/cP/(D/10);

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
Rex   = wx0*rhom0*d0/etam0;
Ref   = wf0*rhom0*d0/etam0;

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
for k = 1:round(smth)
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
            bndinit = (1+erf( ( -ZZ+dw)/h))/2;
        case 2  % bot only
            bndinit = (1+erf(-(D-ZZ-dw)/h))/2;
        case 3  % top/bot only
            bndinit = (1+erf( ( -ZZ+dw)/h))/2 ...
                    + (1+erf(-(D-ZZ-dw)/h))/2;
        case 4 % all walls
            bndinit = (1+erf( ( -ZZ+dw)/h))/2 ...
                    + (1+erf(-(D-ZZ-dw)/h))/2 ...
                    + (1+erf( ( -XX+dw)/h))/2 ...
                    + (1+erf(-(L-XX-dw)/h))/2;
    end
    dw = h;
end
bndinit = max(0,min(1,bndinit));

switch bndmode
    case 0  % none
        bndshape = zeros(size(ZZ(inz,inx)));
    case 1  % top only
        bndshape = exp( ( -ZZ(inz,inx)+h/2)/dw);
    case 2  % bot only
        bndshape = exp(-(D-ZZ(inz,inx)-h/2)/dw);
    case 3  % top/bot only
        bndshape = exp( ( -ZZ(inz,inx)+h/2)/dw) ...
                 + exp(-(D-ZZ(inz,inx)-h/2)/dw);
    case 4 % all walls
        bndshape = exp( ( -ZZ(inz,inx)+h/2)/dw) ...
                 + exp(-(D-ZZ(inz,inx)-h/2)/dw) ...
                 + exp( ( -XX(inz,inx)+h/2)/dw) ...
                 + exp(-(L-XX(inz,inx)-h/2)/dw);
end
bndshape = max(0,min(1,bndshape));

bnd_S = zeros(size(bndshape));
bnd_C = zeros(size(bndshape));
bnd_V = zeros(size(bndshape));

% set specified boundaries to no slip, else to free slip
if bndmode==4;               sds = +1;      % no slip sides for 'all sides(4)'
else;                        sds = -1; end  % free slip sides for other types
if bndmode==1 || bndmode>=3; top = +1;      % no slip top for 'top only(1)', 'top/bot(3)', 'all sides(4)'
else;                        top = -1; end  % free slip for other types
if bndmode>=2;               bot = +1;      % no slip bot for 'bot only(2)', 'top/bot(3)', 'all sides(4)'
else;                        bot = -1; end  % free slip for other types

% initialise solution fields
Tp  =  T0 + (T1-T0) .* (1+erf((ZZ/D-zlay)/wlay_T))/2 + dT.*rp;  if any(bndinit(:)) && ~isnan(Twall); Tp = Tp + (Twall-Tp).*bndinit; end % potential temperature [C]
c   =  c0 + (c1-c0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dc.*rp;  if any(bndinit(:)) && ~isnan(cwall); c  = c  + (cwall-c ).*bndinit; end; cin = c;% major component
v   =  v0 + (v1-v0) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dv.*rp;  if any(bndinit(:)) && ~isnan(vwall); v  = v  + (vwall-v ).*bndinit; end % volatile component

te = zeros(Nz,Nx,length(te0));
for k = 1:length(te0)
    te(:,:,k)  =  te0(k) + (te1(k)-te0(k)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dte(k).*rp;  % trace elements
    if any(bndinit(:)) && ~isnan(tewall(k)); te(:,:,k)  = te(:,:,k) + (tewall(k)-te(:,:,k)).*bndinit; end 
end
ir = zeros(Nz,Nx,length(ir0));
for k = 1:length(ir0)
    ir(:,:,k)  =  ir0(k) + (ir1(k)-ir0(k)) .* (1+erf((ZZ/D-zlay)/wlay_c))/2 + dir(k).*rp;  % isotope ratios  
    if any(bndinit(:)) && ~isnan(irwall(k)); ir(:,:,k)  = ir(:,:,k) + (irwall(k)-ir(:,:,k)).*bndinit; end
end

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
res  = 1;  tol = 1e-12;  x = ones(size(Tp))./10;  f = v/2;
while res > tol
    Pti = Pt;
    
    rhoref =  mean(rho(inz,inx),'all');
    rhofz  = (rho(1:end-1,:)+rho(2:end,:))/2;
    Pt(2:end,:) = Ptop + repmat(cumsum(mean(rhofz,2).*g0.*h),1,Nx);
    Adbt   =  cal.aT./rhoref;
    if Nz<=10; Pt = Ptop.*ones(size(Tp)); end
    
    T    =  (Tp+273.15).*exp(Adbt./cP.*Pt);

    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T-273.15,c,v,Pt,cal,TINY);
    
    c = cin.*(1-f);

    x  = xq;  f = fq;  m = 1-x-f;
    cm = cmq; cx = cxq;
    vm = vmq; vf = vfq;
    Kc = cxq./cmq;
    Kf = vfq./vmq;

    update;

    res  = norm(Pt(:)-Pti(:),2)./norm(Pt(:),2);
end
rhoo = rho;
dto  = dt; 

% get bulk enthalpy, silica, volatile content densities
S  = rho.*(cP.*log(T/(T0+273.15)) + x.*Dsx + f.*Dsf - Adbt.*(Pt-Ptop));  
S0 = rho.*(cP.*log(T0+273.15) + x.*Dsx + f.*Dsf - Adbt.*Ptop);  
C  = rho.*(m.*cm + x.*cx);
V  = rho.*(m.*vm + f.*vf);
X  = rho.*x;
F  = rho.*f;

% get phase entropies
sm = S./rho - x.*Dsx - f.*Dsf;
sx = sm + Dsx;
sf = sm + Dsf;

% get trace element phase compositions
for k = 1:length(te0)
    tem(:,:,k)  = te(:,:,k) ./(m + x.*Kte(k) );
    tex(:,:,k)  = te(:,:,k) ./(m./Kte(k)  + x);
end

% get geochemical component densities
for k = 1:length(te0)
    TE(:,:,k)  = rho.*(m.*tem(:,:,k) + x.*tex(:,:,k));
end
for k = 1:length(ir0)
    IR(:,:,k)  = rho.*ir(:,:,k);
end

% initialise phase change rates
Gx = 0.*x(inz,inx);  Gf = 0.*f(inz,inx);

% initialise auxiliary variables 
dSdt   = 0.*T(inz,inx);  diss_h = 0.*T(inz,inx);
dCdt   = 0.*c(inz,inx);
dVdt   = 0.*v(inz,inx);
dFdt   = 0.*f(inz,inx);
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

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = [opdir,'/',runID,'/',runID,'_cont.mat'];
    elseif restart > 0  % restart from specified continuation frame
        name = [opdir,'/',runID,'/',runID,'_',num2str(restart),'.mat'];
    end
    if exist(name,'file')
        fprintf('\n   restart from %s \n\n',name);
        load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','S','C','V','T','c','v','cm','cx','vm','vf','TE','IR','te','ir','dSdt','dCdt','dVdt','dFdt','dXdt','dTEdt','dIRdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
        
        xq = x; fq = f;
        SOL = [W(:);U(:);P(:)];
        rhoo = rho; Div_rhoVo = Div_rhoV;
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
