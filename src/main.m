load ocean;  % load custom colormap

% get coordinate arrays
X     = -h/2:h:L+h/2;
Z     = -h/2:h:D+h/2;
[XX,ZZ] = meshgrid(X,Z);
Xf    = (X(1:end-1)+X(2:end))./2;
Zf    = (Z(1:end-1)+Z(2:end))./2;
[XXu,ZZu] = meshgrid(Xf,Z);
[XXw,ZZw] = meshgrid(X,Zf);

Nx = length(X);
Nz = length(Z);

% get smoothed initialisation field
rng(15);
rp = randn(Nz,Nx);
for i = 1:round(smth)
    rp(2:end-1,2:end-1) = rp(2:end-1,2:end-1) + diff(rp(:,2:end-1),2,1)./8 + diff(rp(2:end-1,:),2,2)./8;
    rp([1 end],:) = rp([2 end-1],:);
    rp(:,[1 end]) = rp(:,[2 end-1]);
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

% set initial solution fields
T   =  T0 + (T1-T0) .* (1+erf(40*(ZZ/D-zlay)))/2 + dT.*rp;  meanT0 = mean(mean(T(2:end-1,2:end-1)));
c   =  c0 + (c1-c0) .* (1+erf(40*(ZZ/D-zlay)))/2 + dc.*rp;  meanc0 = mean(mean(c(2:end-1,2:end-1)));
v   =  v0 + (v1-v0) .* (1+erf(40*(ZZ/D-zlay)))/2 + dv.*rp;  meanv0 = mean(mean(v(2:end-1,2:end-1))); 
x   =  0.*c + 1e-16;
f   =  0.*c + 1e-16;

si =  si0 + (si1-si0) .* (1+erf(50*(ZZ/D-zlay)))/2 + dsi.*rp;
it =  it0 + (it1-it0) .* (1+erf(50*(ZZ/D-zlay)))/2 + dit.*rp; 
ct =  ct0 + (ct1-ct0) .* (1+erf(50*(ZZ/D-zlay)))/2 + dct.*rp;

U   =  zeros(size((XX(:,1:end-1)+XX(:,2:end))));  Ui = U;  res_U = 0.*U;
W   =  zeros(size((XX(1:end-1,:)+XX(2:end,:))));  Wi = W;  res_W = 0.*W; wf = 0.*W; wc = 0.*W;
P   =  0.*c;  Pi = P;  res_P = 0.*P;  meanQ = 0; Pt = rhom.*g0.*ZZ + Ptop;
S   = [W(:);U(:);P(:)];

% get initial local phase equilibrium
[xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
mq = 1-xq-fq;

% set initial condition to local phase equilibrium
x  = xq;  f = fq;  m = mq;
cm = cmq;
cx = cxq;
vm = vmq;
vf = vfq;

% get volume fractions and bulk density
chi = x;  phi = f;  mu = m;
rho = mu.*rhom + chi.*rhox + phi.*rhof;
res = 1;  tol = 1e-12;
while res > tol
    chii = chi;  phii = phi;
    chi  = max(TINY,min(1-TINY, x.*rho./rhox ))./(chi+phi+mu);
    phi  = max(TINY,min(1-TINY, f.*rho./rhof ))./(chi+phi+mu);
    mu   = max(TINY,min(1-TINY, m.*rho./rhom ))./(chi+phi+mu);
    rho  = mu.*rhom + chi.*rhox + phi.*rhof;
    res  = (norm(chi(:)-chii(:),2) + norm(phi(:)-phii(:),2))./sqrt(2*length(chi(:)));
end
rhoref  = mean(mean(rho(2:end-1,2:end-1)));
Mass0   = sum(sum(rho(2:end-1,2:end-1)));
rhoo    = rho;
  
% get bulk enthalpy, silica, volatile content densities
rhoCp = mu*rhom*Cpm + chi*rhox*Cpx + phi*rhof*Cpf;
H = rhoCp.*T + chi.*rhox.*DLx + phi.*rhof.*DLf;  sumH0 = sum(sum(H(2:end-1,2:end-1)));
C = mu.*rhom.*cm + chi.*rhox.*cx;                sumC0 = sum(sum(C(2:end-1,2:end-1)));
V = mu.*rhom.*vm + phi.*rhof.*vf;                sumV0 = sum(sum(V(2:end-1,2:end-1)));

itm = it./(m + x.*KIT); itx = it./(m./KIT + x);
ctm = ct./(m + x.*KCT); ctx = ct./(m./KCT + x);

IT = mu.*rhom.*itm + chi.*rhox.*itx;
CT = mu.*rhom.*ctm + chi.*rhox.*ctx;

% initialise reaction rates
Gx = 0.*x;  Gf = 0.*f;

% initialise auxiliary variables 
dHdt   =  0.*T;  advn_H = 0.*H;  diff_T = 0.*T;
dCdt   =  0.*c;  advn_C = 0.*C;  diff_c = 0.*c;
dVdt   =  0.*v;  advn_V = 0.*V;  diff_v = 0.*v;  
dfdt   =  0.*f;  diff_f = 0.*f;  diff_it = 0.*it;  diff_ct = 0.*ct;
dxdt   =  0.*x;  diff_x = 0.*x;  diff_si = 0.*si;   
dsidt  =  0.*si; dITdt  = 0.*IT;  dCTdt = 0.*CT;

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_rhomVm = 0.*P;  Div_rhoxVx = 0.*P;  Div_rhofVf = 0.*P;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(size(f)-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(size(f)-1);  tII = 0.*P; 
dVoldt =  0;  VolSrc = 0.*P;  drhodt = 0.*P;  drhodto = 0.*P;

% initialise timing and iterative parameters
step   =  0;
time   =  0;
iter   =  0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = ['../out/',runID,'/',runID,'_cont'];
    elseif restart > 0  % restart from specified continuation frame
        name = ['../out/',runID,'/',runID,'_',num2str(restart)];
    end
    load(name,'U','W','P','Pt','f','x','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','it','ct','si','dHdt','dCdt','dVdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step');
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
    
    xq = x;  cmq = cm;  cxq = cx; 
    fq = f;  vmq = vm;  vfq = vf;
    time = time+dt;
    step = step+1;
end

update; output;

% physical time stepping loop
while time <= tend && step <= M
        
    fprintf(1,'\n\n\n*****  step %d;  dt = %4.4e;  time = %4.4e [hr]\n\n',step,dt./3600,time./3600);
    tic;
    
    % store previous solution
    Ho      = H;
    Co      = C;
    Vo      = V;
    To      = T;
    co      = c;
    vo      = v;
    xo      = x;
    fo      = f;
    sio     = si;
    ITo     = IT;
    CTo     = CT;
    rhoo    = rho;
    dHdto   = dHdt;
    dCdto   = dCdt;
    dVdto   = dVdt;
    dxdto   = dxdt;
    dfdto   = dfdt;
    dITdto  = dITdt;
    dCTdto  = dCTdt;
    dsidto  = dsidt;
    drhodto = drhodt;
    Pto     = Pt;
    dto     = dt;
    
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    iter     = 0;
    
    % non-linear iteration loop
    while resnorm/resnorm0 >= rtol && resnorm >= atol && iter <= maxit 
            
        % solve thermo-chemical equations
        thermochem;
        
        % update non-linear parameters and auxiliary variables
        update;
        
        % solve fluid-mechanics equations
        fluidmech;
        
        % update non-linear parameters and auxiliary variables
        update;
        
        % report convergence
        report;

        iter = iter+1;
    end
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);
    
    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)  ),mean(T(:)  ),max(T(:)  ));
    fprintf(1,'         min c   =  %1.4f;    mean c   = %1.4f;    max c   = %1.4f;   [wt]\n'   ,min(c(:)  ),mean(c(:)  ),max(c(:)  ));
    fprintf(1,'         min v   =  %1.4f;    mean v   = %1.4f;    max v   = %1.4f;   [wt]\n\n' ,min(v(:)  ),mean(v(:)  ),max(v(:)  ));
    
    fprintf(1,'         min x   =  %1.4f;    mean x   = %1.4f;    max x   = %1.4f;   [wt]\n'   ,min(x(:)  ),mean(x(:)  ),max(x(:)  ));
    fprintf(1,'         min f   =  %1.4f;    mean f   = %1.4f;    max f   = %1.4f;   [wt]\n'   ,min(f(:)  ),mean(f(:)  ),max(f(:)  ));
    fprintf(1,'         min m   =  %1.4f;    mean m   = %1.4f;    max m   = %1.4f;   [wt]\n\n' ,min(m(:)  ),mean(m(:)  ),max(m(:)  ));

    fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
    fprintf(1,'         min eta =  %1.2e;  mean eta = %1.2e;  max eta = %1.2e; [Pas]\n\n',min(eta(:)),geomean(eta(:)),max(eta(:)));

    fprintf(1,'         min U   = %1.4f;    mean U   = %1.4f;    max U   = %1.4f;   [m/s]\n'  ,min(U(:)  ),mean(U(:)  ),max(U(:)  ));
    fprintf(1,'         min W   = %1.4f;    mean W   = %1.4f;    max W   = %1.4f;   [m/s]\n'  ,min(-W(:) ),mean(-W(:) ),max(-W(:) ));
    fprintf(1,'         min P   = %2.4f;    mean P   = %2.4f;    max P   = %2.4f;  [kPa]\n\n',min(P(:)./1e3),mean(P(:)./1e3),max(P(:)./1e3));

    % plot results
    if ~mod(step,nop); output; end
    
    % increment time/step
    time = time+dt;
    step = step+1;
    
end

diary off
