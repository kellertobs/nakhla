% calibrate phase diagram
addpath('../src');
load ocean
TINY = 1e-16;
run('../usr/par_default');       % load default parameters


% calibration run options
runID     = 'test';              % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set phase diagram parameters
calID    =  'default';           % phase diagram calibration
run(['../cal/cal_',calID]);      % load melt model calibration

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

% set ranges for control variables T, c, v, P
T = linspace(1300,850,400).';    % temperature range [degC]
c = linspace(0.488,0.488,400).';   % major component range [wt SiO2]
v = linspace(0.00,0.00,400).';   % volatile component range [wt H2O]
P = linspace(125,125,400).'*1e6; % pressure range [Pa]

% equilibrium phase fractions and compositions
c0 = c; res = 1; x = zeros(size(T)); f = zeros(size(T));
while res>1e-13
    ci = c;
    [x,cx,cm,f,vf,vm]  =  equilibrium(x,f,T,c,v,P,cal,TINY);
    c = c0.*(1-f);
    m = 1-f-x;
    res = norm(c-ci,'fro')./norm(c,'fro');
end

Nz = length(T); Nx = 1; Ptop = min(P); Pt = P; etareg = 1; calibrt = 1; T = T+273.15;
update;
T = T-273.15;

wm =  mu.*Ksgr_m .* (rhom-rho)*g0; % melt segregation speed
wx = chi.*Ksgr_x .* (rhox-rho)*g0; % crystal segregation speed
wf = phi.*Ksgr_f .* (rhof-rho)*g0; % fluid segregation speed


% block out values below solidus and above liquidus
ind = x<1e-9 | m<1e-9;
c_oxd  = squeeze(c_oxd); c_oxd(ind,:) = [];
cm_oxd = squeeze(cm_oxd); cm_oxd(ind,:) = [];
cx_oxd = squeeze(cx_oxd); cx_oxd(ind,:) = [];
c_cmp  = squeeze(c_cmp); c_cmp(ind,:) = [];
cm_cmp = squeeze(cm_cmp); cm_cmp(ind,:) = [];
cx_cmp = squeeze(cx_cmp); cx_cmp(ind,:) = [];
c_mem  = squeeze(c_mem); c_mem(ind,:) = [];
cm_mem = squeeze(cm_mem); cm_mem(ind,:) = [];
cx_mem = squeeze(cx_mem); cx_mem(ind,:) = [];
cx_msy = squeeze(cx_msy); cx_msy(ind,:) = [];
cx (ind) = [];
cm (ind) = [];
vf (ind) = [];
vm (ind) = [];
c   (ind) = [];
v   (ind) = [];
T   (ind) = [];
P   (ind) = [];
rhox(ind) = [];
rhom(ind) = [];
rhof(ind) = [];
rho (ind) = [];
etam(ind) = [];
eta (ind) = [];
wm  (ind) = [];
wx  (ind) = [];
wf  (ind) = [];
f  (ind) = [];
x  (ind) = [];
m  (ind) = [];
phi(ind) = [];
chi(ind) = [];
mu (ind) = [];
Ptop = min(P); Pt = P;

if ~holdfig; close all; end

% plot phase diagram
figure(1); if ~holdfig; clf; end
TT = [linspace(cal.Tphs0+Ptop*cal.clap,cal.perT+Ptop*cal.clap,200),linspace(cal.perT+Ptop*cal.clap,cal.Tphs1+Ptop*cal.clap,200)];
cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2,200),linspace((cal.perCx+cal.perCm)/2,cal.cphs0,200)];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,Ptop*ones(size(TT)),cal,TINY);
plot(CCx,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k-','LineWidth',2);
perTs  = cal.perT;
Tphs0s = cal.Tphs0;
Tphs1s = cal.Tphs1;
vv = 0.10*ones(size(TT));
xx = 0.50*ones(size(TT));
ff = 0.05*ones(size(TT));
for i = 1:5
    TT = [linspace(Tphs0s+Ptop*cal.clap,perTs+Ptop*cal.clap,200),linspace(perTs+Ptop*cal.clap,Tphs1s+Ptop*cal.clap,200)];
    cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2,200),linspace((cal.perCx+cal.perCm)/2,cal.cphs0,200)];
    vmq_c0 = (4.7773e-7.*Ptop.^0.6 + 1e-11.*Ptop) .* exp(2565*(1./(TT+273.15)-1./(cal.perT+273.15))); % Katz et al., 2003; Moore et al., 1998
    vmq_c1 = (3.5494e-3.*Ptop.^0.5 + 9.623e-8.*Ptop - 1.5223e-11.*Ptop.^1.5)./(TT+273.15) + 1.2436e-14.*Ptop.^1.5; % Liu et al., 2015
    vmq0   = (1-cc).*vmq_c0 + cc.*vmq_c1;
    perTs  = cal.perT -cal.dTH2O(2).*vmq0(round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*200)).^0.75;
    Tphs0s = cal.Tphs0-cal.dTH2O(1).*vmq0(1  ).^0.75;
    Tphs1s = cal.Tphs1-cal.dTH2O(3).*vmq0(end).^0.75;
end
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,vmq0,Ptop*ones(size(TT)),cal,TINY);
plot(CCx,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k-','LineWidth',2);

Tplt = T - (Pt-Ptop)*cal.clap;
cplt = c./(1-f);
plot(cplt,Tplt,'k.','LineWidth',2,'MarkerSize',15);
plot(cx  ,Tplt,'b.','LineWidth',2,'MarkerSize',15);
plot(cm  ,Tplt,'r.','LineWidth',2,'MarkerSize',15);

title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,x.*100,'k',T,m.*100,'r',T,f.*1000,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',15,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Melting model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Phase fractions [wt\%]','Interpreter','latex','FontSize',15)

% plot major phase compositions
figure(3); if ~holdfig; clf; end
plot(T,cx.*100,'b',T,cm.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vf/10.*100,'b',T,vm.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('fluid $/10$','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Volatile component [wt\% H$_2$O]','Interpreter','latex','FontSize',15)

% plot phase densities
figure(5); if ~holdfig; clf; end
plot(T,rhox,'k',T,rhom,'r',T,rhof,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
plot(T,rho ,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','mixture','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Density model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Density [kg/m$^3$]','Interpreter','latex','FontSize',15)

% plot mixture rheology
figure(6); if ~holdfig; clf; end
semilogy(T,eta,'k',T,etam,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('mixture','melt','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Viscosity model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Viscosity [log$_{10}$ Pas]','Interpreter','latex','FontSize',15)

% plot phase segregation speeds
figure(7); if ~holdfig; clf; end
semilogy(T,max(1e-18,abs(chi.*wx)).*3600,'k',T,max(1e-18,abs(mu.*wm)).*3600,'r',T,max(1e-18,abs(phi.*wf)).*3600,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase segregation model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Segregation flux [m/hr]','Interpreter','latex','FontSize',15)

% plot oxide compositions
figure(14); clf;
subplot(2,1,1)
sgtitle('Phase Oxide Fractions','Interpreter','latex','FontSize',18)
for i=1:cal.noxd
    plot(T,cm_oxd(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.noxd)+1,:)); hold on; box on; axis tight;
end
legend(cal.oxdStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
for i=1:cal.noxd
    plot(T,cx_oxd(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.noxd)+1,:)); hold on; box on; axis tight;
end
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

% plot end-member component compositions
figure(15); clf;
subplot(2,1,1)
sgtitle('Phase Component Fractions','Interpreter','latex','FontSize',18)
for i=1:cal.ncmp
    plot(T,cm_cmp(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); hold on; box on; axis tight;
end
legend(cal.cmpStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
for i=1:cal.ncmp
    plot(T,cx_cmp(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); hold on; box on; axis tight;
end
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

% plot simplified mineral assemblage
figure(16); clf;
patch([T;flipud(T)],[zeros(size(T));flipud(cx_msy(:,1))],[0.6,0.8,0.5],'LineWidth',2); hold on; box on; axis tight;
patch([T;flipud(T)],[sum(cx_msy(:,1  ),2);flipud(sum(cx_msy(:,1:2),2))],0.7.*[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:2),2);flipud(sum(cx_msy(:,1:3),2))],[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:3),2);flipud(sum(cx_msy(:,1:4),2))],[0.9,0.9,0.9],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:4),2);flipud(sum(cx_msy(:,1:5),2))],[0.9,0.7,0.9],'LineWidth',2);
legend('~olivine','~oxide','~pyroxene','~feldspar','~quartz','Interpreter','latex','FontSize',15,'box','off','location','southeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Mineral assemblage','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Mineral fraction [wt]','Interpreter','latex','FontSize',15)


% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save output to file
if save_plot
    name = ['../out/',runID,'/',runID,'_phase_dgrm'];
    print(figure(1),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_melt_model'];
    print(figure(2),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_maj_compnt'];
    print(figure(3),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_vol_compnt'];
    print(figure(4),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_density'];
    print(figure(5),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_viscosity'];
    print(figure(6),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_segr_speed'];
    print(figure(7),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_oxides'];
    print(figure(8),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_components'];
    print(figure(9),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_modalcmp'];
    print(figure(10),name,'-dpng','-r300','-opengl');
end
