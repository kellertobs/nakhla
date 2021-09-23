load ocean;  % load custom colormap

% get coordinate arrays
X     = -h/2:h:L+h/2;
Z     = -h/2:h:D+h/2;
[XX,ZZ] = meshgrid(X,Z);
Xc    = (X(1:end-1)+X(2:end))./2;
Zc    = (Z(1:end-1)+Z(2:end))./2;

Nx = length(X);
Nz = length(Z);

% get smoothed initialisation field
rng(15);
a = randn(Nz,Nx);
for i = 1:round(smth)
    a(2:end-1,2:end-1) = a(2:end-1,2:end-1) + diff(a(:,2:end-1),2,1)./8 + diff(a(2:end-1,:),2,2)./8;
    a([1 end],:) = a([end-1 2],:);
    a(:,[1 end]) = a(:,[2 end-1]);
end
a = a./max(abs(a(:)));

% get mapping arrays
NP =  Nz   * Nx   ;
NW = (Nz-1)* Nx   ;
NU =  Nz   *(Nx-1);
MapP = reshape(1:NP,Nz  ,Nx  );
MapW = reshape(1:NW,Nz-1,Nx  );
MapU = reshape(1:NU,Nz  ,Nx-1) + NW;

% set initial solution fields
T   =  T0   + a.*T1;  To   = T;  cool = 0.*T;
c   =  c0 + a.*c1;  co = c;
v   =  v0 + a.*v1;  vo = v;
x   =  0.*c + 1e-16;
f   =  0.*c + 1e-16;
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
    chi  = max(TINY,min(1-TINY, x.*rho./rhox ));
    phi  = max(TINY,min(1-TINY, f.*rho./rhof ));
    mu   = 1-chi-phi;
    rho  = mu.*rhom + chi.*rhox + phi.*rhof;
    res  = (norm(chi(:)-chii(:),2) + norm(phi(:)-phii(:),2))./sqrt(2*length(chi(:)));
end

% get bulk enthalpy, silica, volatile content densities
Cprhob = mu*rhom*Cpm + chi*rhox*Cpx + phi*rhof*Cpf;
H = Cprhob.*T + chi.*rhox.*DLx + phi.*rhof.*DLf;
C = mu.*rhom.*cm + chi.*rhox.*cx;
V = mu.*rhom.*vm + phi.*rhof.*vf;

% initialise reaction rates
Gx = 0.*x;  Gf = 0.*f;

% initialise auxiliary variables 
dHdt   =  0.*T;  diff_T = 0.*T; 
dCdt   =  0.*c;  diff_c = 0.*c;  
dVdt   =  0.*v;  diff_v = 0.*v;  
% dfdt   =  0.*f;  diff_f = 0.*f;  
% dxdt   =  0.*x;  diff_x = 0.*x; 

eIIref =  1e-6;  
Div_V  =  0.*P;  Div_mVm = 0.*P;  Div_xVx = 0.*P;  Div_fVf = 0.*P; dwfdz = 0.*P;  dwxdz = 0.*P;
exx    =  0.*P;  ezz = 0.*P;  exz = zeros(size(f)-1);  eII = 0.*P;  
txx    =  0.*P;  tzz = 0.*P;  txz = zeros(size(f)-1);  tII = 0.*P; 

% initialise timing and iterative parameters
step   =  0;
time   =  0;
it     =  0;

% overwrite fields from file if restarting run
if restart
    if     restart < 0  % restart from last continuation frame
        name = ['../out/',runID,'/',runID,'_cont'];
    elseif restart > 0  % restart from specified continuation frame
        name = ['../out/',runID,'/',runID,'_',num2str(restart)];
    end
    load(name,'U','W','P','Pt','f','x','phi','chi','T','c','v','cm','cx','vm','vf','dTdt','dcdt','dvdt','dfdt','dxdt','Gf','Gx','rho','eta','Div_V','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step');
    name = ['../out/',runID,'/',runID,'_par'];
    load(name);
    
    % update equilibrium
    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);

    time = time+dt;
    step = step+1;
end
        
% initialise material properties and plot initial condition
update;  output;

% physical time stepping loop
while time <= tend && step <= M
        
    fprintf(1,'\n\n\n*****  step %d;  dt = %4.4e;  time = %4.4e;\n\n',step,dt,time);
    tic;
    
    % store previous solution
    Ho      = H;
    Co      = C;
    Vo      = V;
    xo      = x;
    fo      = f;
    rhoo    = rho;
    dHdto   = dHdt;
    dCdto   = dCdt;
    dVdto   = dVdt;
    dxdto   = dxdt;
    dfdto   = dfdt;
    dto     = dt;
    
    % reset residuals and iteration count
    resnorm  = 1e3;
    resnorm0 = resnorm;
    it       = 0;
    
    % non-linear iteration loop
    startup = 2*double(step<=0) + double(step>0);
    while resnorm/resnorm0 >= rtol^startup && resnorm >= atol/startup && it <= maxit*startup
            
        % solve thermo-chemical equations
        thermochem;
        
        % update nonlinear coefficients & auxiliary fields
        update;
        
        % solve fluid-mechanics equations
        fluidmech;
        
        % report convergence
        report;

        it = it+1;
    end
    
    % print diagnostics
    fprintf(1,'\n         time to solution = %4.4f sec\n\n',toc);
    
    fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)  ),mean(T(:)  ),max(T(:)  ));
    fprintf(1,'         min c   =  %2.3f;    mean c   = %2.3f;    max c   = %2.3f;   [wt]\n'   ,min(c(:)  ),mean(c(:)  ),max(c(:)  ));
    fprintf(1,'         min v   =  %1.4f;    mean v   = %1.4f;    max v   = %1.4f;   [wt]\n'   ,min(v(:)  ),mean(v(:)  ),max(v(:)  ));
    fprintf(1,'         min x   =  %1.4f;    mean x   = %1.4f;    max x   = %1.4f;   [wt]\n'   ,min(x(:)  ),mean(x(:)  ),max(x(:)  ));
    fprintf(1,'         min f   =  %1.4f;    mean f   = %1.4f;    max f   = %1.4f;   [wt]\n\n' ,min(f(:)  ),mean(f(:)  ),max(f(:)  ));

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
