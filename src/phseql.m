%*** update phase equilibrium
eqtime = tic;

var.c      = reshape(c,Nx*Nz,cal.ncmp);   % component fractions [wt]
var.T      = reshape(T,Nx*Nz,1)-273.15;   % temperature [C]
var.P      = reshape(Pt,Nx*Nz,1)/1e9;     % pressure [GPa]
var.m      = reshape(mq,Nx*Nz,1);         % melt fraction [wt]
var.f      = reshape(fq,Nx*Nz,1);         % bubble fraction [wt]
var.H2O    = var.c(:,end);                % water concentration [wt]
var.X      = reshape(cm_oxd_all,Nz*Nx,9); % melt oxide fractions [wt %]
cal.H2Osat = fluidsat(var);               % water saturation [wt]

[var,cal]  = leappart(var,cal,'E');

Tsol   = reshape(cal.Tsol,Nz,Nx);
Tliq   = reshape(cal.Tliq,Nz,Nx);
H2Osat = reshape(cal.H2Osat,Nz,Nx);

mq = reshape(var.m.*(var.m>eps^0.5),Nz,Nx);
fq = reshape(var.f.*(var.f>eps^0.5),Nz,Nx);
xq = reshape(var.x.*(var.x>eps^0.5),Nz,Nx);
mq = mq./(mq+xq+fq);
xq = xq./(mq+xq+fq);
fq = fq./(mq+xq+fq);

% phase mass transfer rates
Gm  = (mq-m).*RHO/max(tau_r,5*dt);
Gx  = (xq-x).*RHO/max(tau_r,5*dt); 
Gf  = (fq-f).*RHO/max(tau_r,5*dt);

Gmc = (cmq.*mq-cm.*m).*RHO/max(tau_r,5*dt);
Gxc = (cxq.*xq-cx.*x).*RHO/max(tau_r,5*dt);
Gfc = (cfq.*fq-cf.*f).*RHO/max(tau_r,5*dt);

cxq = reshape(var.cx,Nz,Nx,cal.ncmp);
cmq = reshape(var.cm,Nz,Nx,cal.ncmp);

eqtime = toc(eqtime);
EQtime = EQtime + eqtime;