% create manufactured solution
load ./colmap/ocean.mat
clear x z
TINY = 1e-16;
syms U_mms(x,z) W_mms(x,z) P_mms(x,z) eta_mms(x,z) rho_mms(x,z) src_mms(x,z)

fprintf(1,'\n\n  ***  compose manufactured solution\n\n');

% compose manufactured solution variables
W_mms(x,z) = 5e-5.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L));
U_mms(x,z) = 4e-5.*(sin(4*(x)*pi/L).*cos(4*(z)*pi/L));
P_mms(x,z) =-3e+3.*(cos(4*(x)*pi/L).*cos(4*(z)*pi/L));

% compose manufactured material coefficients and volume source
eta_mms(x,z) = 1e+3-9e+2.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L));
rho_mms(x,z) = 3e+3-5e+1.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L)); rhoref = 3e+3;
src_mms(x,z) =     -1e-3.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L));

fprintf(1,'       W   = %s \n',char(W_mms));
fprintf(1,'       U   = %s \n',char(U_mms));
fprintf(1,'       P   = %s \n',char(P_mms));
fprintf(1,'       eta = %s \n',char(eta_mms));
fprintf(1,'       rho = %s \n',char(rho_mms));
fprintf(1,'       src = %s \n',char(src_mms));
fprintf(1,'       . ');

% update strain rates
DivV_mms(x,z)= (diff(W_mms,z) + diff(U_mms,x));
exx_mms(x,z) = diff(U_mms,x) - DivV_mms./3;         % x-normal strain rate
ezz_mms(x,z) = diff(W_mms,z) - DivV_mms./3;         % z-normal strain rate
exz_mms(x,z) = 1/2.*(diff(U_mms,z)+diff(W_mms,x));  % xz-shear strain rate
fprintf(1,' . ');

% update stresses
txx_mms(x,z) = eta_mms .* exx_mms;                  % x-normal stress
tzz_mms(x,z) = eta_mms .* ezz_mms;                  % z-normal stress
txz_mms(x,z) = eta_mms .* exz_mms;                  % xz-shear stress
fprintf(1,' . ');

% manufactured solution residuals
res_W_mms = -(diff(tzz_mms,z) + diff(txz_mms,x)) + diff(P_mms,z) - (rho_mms(x,z)-rhoref)*g0;
res_U_mms = -(diff(txx_mms,x) + diff(txz_mms,z)) + diff(P_mms,x);
res_P_mms =  DivV_mms - src_mms(x,z);
fprintf(1,' . ');

% plot manufactured solution
figure(15);
colormap(ocean);
subplot(2,3,1); fcontour( -W_mms*hr  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufcat. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); fcontour(  U_mms*hr  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); fcontour(  P_mms/1e3 ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
fprintf(1,' . ');
subplot(2,3,4); fcontour(      rho_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\rho$ [kg/m$^3$]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,5); fcontour(log10(eta_mms),[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\eta$ [log$_{10}$ Pas]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,6); fcontour(      src_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\dot{V} [1/s]$','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
drawnow;
fprintf(1,' . \n');

% evaluate mms source terms on appropriate coordinate grids
fprintf(1,'\n  ***  evaluate manufactured solution\n\n');
x_mms  = -h/2:h:L+h/2;
z_mms  = -h/2:h:L+h/2;
xu_mms = (x_mms(1:end-1)+x_mms(2:end))./2;
zw_mms = (z_mms(1:end-1)+z_mms(2:end))./2;

fprintf(1,'       Patience, my young Padawan!\n');
fprintf(1,'       . ');

[x,z] = meshgrid(x_mms,zw_mms);
src_W_mms = double(subs(res_W_mms)); fprintf(1,' . ');
[x,z] = meshgrid(xu_mms,z_mms);
src_U_mms = double(subs(res_U_mms)); fprintf(1,' . ');
[x,z] = meshgrid(x_mms,z_mms);
src_P_mms = double(subs(res_P_mms)); fprintf(1,' . ');

% plot manufactured residuals and evaluated source terms
figure(16);
colormap(ocean);
subplot(2,3,1); fcontour(-res_W_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $W$-res','Interpreter','latex');
subplot(2,3,2); fcontour(-res_U_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $U$-res','Interpreter','latex');
subplot(2,3,3); fcontour(-res_P_mms,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar; box on; title('manufactured $P$-res','Interpreter','latex');
subplot(2,3,4); imagesc(x_mms,zw_mms,-src_W_mms); axis ij equal tight; colorbar; box on; title('evaluated $W$-res','Interpreter','latex');
subplot(2,3,5); imagesc(xu_mms,z_mms,-src_U_mms); axis ij equal tight; colorbar; box on; title('evaluated $U$-res','Interpreter','latex');
subplot(2,3,6); imagesc(x_mms ,z_mms,-src_P_mms); axis ij equal tight; colorbar; box on; title('evaluated $P$-res','Interpreter','latex');
drawnow;

% evaluate analytical solution on appropriate coordinate grids
[x,z]  = meshgrid(x_mms,zw_mms);
W_mms  = double(subs(W_mms)); fprintf(1,' . ');
rhow   = double(subs(rho_mms)); fprintf(1,' . ');
rhow   = rhow(:,2:end-1);
[x,z]  = meshgrid(xu_mms,z_mms);
U_mms  = double(subs(U_mms)); fprintf(1,' . ');
rhou   = double(subs(rho_mms)); fprintf(1,' . ');
rhou   = rhou(2:end-1,:);
[x,z]  = meshgrid(x_mms,z_mms);
P_mms  = double(subs(P_mms)); fprintf(1,' . ');
eta    = double(subs(eta_mms)); fprintf(1,' . ');
eta    = eta(2:end-1,2:end-1);
dV     = double(subs(src_mms)); fprintf(1,' . ');
dV     = dV(2:end-1,2:end-1);
[x,z]  = meshgrid(xu_mms,zw_mms);
etaco  = double(subs(eta_mms)); fprintf(1,' . ');

Drho   = rhow-mean(rhow,2);
rhoWo  = zeros(size(rhow));
rhoWoo = zeros(size(rhow));
rhoUo  = zeros(size(rhou));
rhoUoo = zeros(size(rhou));
WBG    = 0.*W_mms;  W = WBG;
UBG    = 0.*U_mms;  U = UBG;
SOL    = [W_mms(:);U_mms(:);P_mms(:)];
dt     = 1e32;

% get mapping arrays
Nx = length(x_mms)-2;
Nz = length(z_mms)-2;

NP = (Nz+2) * (Nx+2);
NW = (Nz+1) * (Nx+2);
NU = (Nz+2) * (Nx+1);
MapP = reshape(1:NP,Nz+2,Nx+2);
MapW = reshape(1:NW,Nz+1,Nx+2);
MapU = reshape(1:NU,Nz+2,Nx+1) + NW;

% set time stepping parameters
a1 = 1; a2 = 1; a3 = 0;
b1 = 1; b2 = 0; b3 = 0;

% set boundary conditions to free slip
sds = -1;
top = -1;
bot = -1;
BCA = {'',''};

% set ghosted index arrays
if periodic
    icx = [Nx,1:Nx,1];
    icz = [Nz,1:Nz,1];
    ifx = [Nx,1:Nx+1,2];
    ifz = [Nz,1:Nz+1,2];
else
    icx = [1,1:Nx,Nx];
    icz = [1,1:Nz,Nz];
    ifx = [1,1:Nx+1,Nx+1];
    ifz = [1,1:Nz+1,Nz+1];
end

fprintf(1,' . \n');
   
% call fluid mechanics solver
FMtime = 0;
fluidmech;
