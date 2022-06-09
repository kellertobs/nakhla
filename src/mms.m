% create manufactured solution
clear x z;
syms U_mms(x,z) W_mms(x,z) P_mms(x,z) eta_mms(x,z) rho_mms(x,z) src_mms(x,z)

fprintf(1,'\n\n  ***  compose manufactured solution\n\n');

% compose manufactured solution variables
W_mms(x,z) = 5.00e-5.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L));
U_mms(x,z) = 1.00e-5.*(sin(4*(x)*pi/L).*cos(4*(z)*pi/L));
P_mms(x,z) = 2.00e+4.*(sin(4*(x)*pi/L).*sin(4*(z)*pi/L));

% compose manufactured material coefficients and volume source
eta_mms(x,z) = 1e+3-8e+2.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L));
rho_mms(x,z) = 3e+3-1e+1.*(cos(4*(x)*pi/L).*sin(4*(z)*pi/L)); rhoref = 3e+3;
src_mms(x,z) =     -1e-4.*(sin(4*(x)*pi/L).*sin(4*(z)*pi/L));

fprintf(1,'       W   = %s \n',char(W_mms));
fprintf(1,'       U   = %s \n',char(U_mms));
fprintf(1,'       P   = %s \n',char(P_mms));
fprintf(1,'       eta = %s \n',char(eta_mms));
fprintf(1,'       rho = %s \n',char(rho_mms));
fprintf(1,'       src = %s \n',char(src_mms));
fprintf(1,'       . ');

% update strain rates
exx_mms(x,z) = diff(U_mms,x) - (diff(W_mms,z) + diff(U_mms,x))./3;         % x-normal strain rate
ezz_mms(x,z) = diff(W_mms,z) - (diff(W_mms,z) + diff(U_mms,x))./3;         % z-normal strain rate
exz_mms(x,z) = 1/2.*(diff(U_mms,z)+diff(W_mms,x))                ;         % shear strain rate
fprintf(1,' . ');

% update stresses
txx_mms(x,z) = eta_mms .* exx_mms;                                         % x-normal stress
tzz_mms(x,z) = eta_mms .* ezz_mms;                                         % z-normal stress
txz_mms(x,z) = eta_mms .* exz_mms;                                         % xz-shear stress
fprintf(1,' . ');

% manufactured solution residuals
res_W_mms = (diff(tzz_mms,z) + diff(txz_mms,x)) - diff(P_mms,z) + (rho_mms(x,z)-rhoref)*g0;
res_U_mms = (diff(txx_mms,x) + diff(txz_mms,z)) - diff(P_mms,x);
res_P_mms =-(diff(  W_mms,z) + diff(  U_mms,x))                 + src_mms(x,z);
fprintf(1,' . ');

% plot manufactured solution
figure(15);
colormap(ocean);
subplot(2,3,1); fcontour( -W_mms*hr  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufcat. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,2); fcontour(  U_mms*hr  ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
subplot(2,3,3); fcontour(  P_mms/1e3 ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
fprintf(1,' . ');
subplot(2,3,4); fcontour(      rho_mms ,[0,L],'LineWidth',1.5); axis ij equal tight; colorbar('TicklabelInterpreter','latex'); box on; title('manufact. $\Delta \rho$ [kg/m$^3$]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
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
[x,z]  = meshgrid(xu_mms,z_mms);
U_mms  = double(subs(U_mms)); fprintf(1,' . ');
[x,z]  = meshgrid(x_mms,z_mms);
P_mms  = double(subs(P_mms)); fprintf(1,' . ');
rho    = double(subs(rho_mms)); fprintf(1,' . ');
eta    = double(subs(eta_mms)); fprintf(1,' . ');
VolSrc = double(subs(src_mms)); fprintf(1,' . ');
[x,z]  = meshgrid(xu_mms,zw_mms);
etaco  = double(subs(eta_mms)); fprintf(1,' . ');
theta  = 1;
WBG    = 0.*WBG;
UBG    = 0.*UBG;

fprintf(1,' . \n');
   
% call fluid mechanics solver
tic;
fluidmech;
solvetime = toc;
