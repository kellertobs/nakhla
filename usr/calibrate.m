% calibrate phase diagram
clear; addpath('../src');
TINY = 1e-16;


% calibration run options
runID     = 'test';              % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set phase diagram parameters
calID    =  'krafla';            % phase diagram calibration

% set model rheology parameters
etaf0    =  0.1;                 % fluid viscosity [Pas]
etax0    =  1e16;                % crystal viscosity [Pas]
AA       = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
BB       = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
CC       = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% set model buoyancy parameters
rhom0    =  2750;                % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhox0    =  3050;                % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhof0    =  1000;                % bubble phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aT       =  4e-5;                % thermal expansivity [1/K]
gC       =  0.5;                 % compositional expansivity [1/wt]
bP       =  1e-8;                % mvp compressibility [1/Pa]
dx       =  5e-4;                % crystal size [m]
df       =  5e-4;                % bubble size [m]
dm       =  5e-4;                % melt film size [m]
g0       =  10.;                 % gravity [m/s2]

% set ranges for control variables T, c, v, P
T = linspace(500,1600,1e3).';    % temperature range [degC]
c = linspace(0.51,0.51,1e3).';   % major component range [wt SiO2]
v = linspace(0.04,0.04,1e3).';   % volatile component range [wt H2O]
P = linspace(125,125,1e3).'*1e6; % pressure range [Pa]

% equilibrium phase fractions and compositions
run(['../cal/cal_',calID]);  % load melt model calibration
[x,cx,cm,f,vf,vm]  =  equilibrium(ones(size(T)).*0.5,ones(size(T)).*0.0,T,c,v,P,cal,TINY);
m = 1-f-x;  

Nz = length(T); Nx = 1; Ptop = min(P); Pt = P; etareg = 1; calibrt = 1; T = T+273.15;
update;
T = T-273.15;

wm = mu .*Csgr_m .* (rhom-rho)*g0; % crystal segregation speed
wx = chi.*Csgr_x .* (rhox-rho)*g0; % crystal segregation speed
wf = phi.*Csgr_f .* (rhof-rho)*g0; % fluid segregation speed

% simplified mineral assemblage (ol+px+fs+qz)
olv0 = cal.perCx; olv1 = cal.cphs0; 
pxn0 = cal.cphs0; pxn1 = cal.perCx;  pxn2 = cal.cphs1;
plg0 = (cal.cphs0+cal.perCx)/2;  plg1 = cal.perCm;  plg2 = (cal.cphs1+1)/2;
kfs0 = cal.perCm;  kfs1 = cal.cphs1;  kfs2 = 1;
qtz0 = cal.perCx;  qtz1 = 1;

olv = max(0,min(1,(cx-olv0)./(olv1-olv0))).*(cx>=olv1 & cx<=olv0);
pxn = max(0,min(1,(cx-pxn0)./(pxn1-pxn0))).*(cx>=pxn0 & cx<=pxn1) + max(0,min(1,(cx-pxn2)./(pxn1-pxn2))).*(cx>=pxn1 & cx<=pxn2);
plg = max(0,min(1,(cx-plg0)./(plg1-plg0))).*(cx>=plg0 & cx<=plg1) + max(0,min(1,(cx-plg2)./(plg1-plg2))).*(cx>=plg1 & cx<=plg2);
kfs = max(0,min(1,(cx-kfs0)./(kfs1-kfs0))).*(cx>=kfs0 & cx<=kfs1) + max(0,min(1,(cx-kfs2)./(kfs1-kfs2))).*(cx>=kfs1 & cx<=kfs2);
qtz = max(0,min(1,(cx-qtz0)./(qtz1-qtz0))).*(cx>=qtz0 & cx<=qtz1);
all = olv+pxn+plg+kfs+qtz;
olv = olv./all;
pxn = pxn./all;
plg = plg./all;
kfs = kfs./all;
qtz = qtz./all;

% block out values below solidus and above liquidus
ind = x<1e-9 | m<1e-9;
c_oxds = squeeze(c_oxds); c_oxds(ind,:) = [];
cm_oxds = squeeze(cm_oxds); cm_oxds(ind,:) = [];
cx_oxds = squeeze(cx_oxds); cx_oxds(ind,:) = [];
cx (ind) = [];
cm (ind) = [];
vf (ind) = [];
vm (ind) = [];
olv (ind) = [];
pxn (ind) = [];
plg (ind) = [];
kfs (ind) = [];
qtz (ind) = [];
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
Ptop = min(P); Pt = P;

if ~holdfig; close all; end

% plot phase diagram
figure(1); if ~holdfig; clf; end
TT = linspace(cal.Tphs0+Ptop*cal.clap,cal.Tphs1+Ptop*cal.clap,500);
cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2, ceil((cal.perT-cal.Tphs0)./(cal.Tphs1-cal.Tphs0)*500)), ...
      linspace((cal.perCx+cal.perCm)/2,cal.cphs0,floor((cal.perT-cal.Tphs1)./(cal.Tphs0-cal.Tphs1)*500))];
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
    TT = linspace(Tphs0s+Ptop*cal.clap,Tphs1s+Ptop*cal.clap,500);
    cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2, ceil((perTs-Tphs0s)./(Tphs1s-Tphs0s)*500)), ...
          linspace((cal.perCx+cal.perCm)/2,cal.cphs0,floor((perTs-Tphs1s)./(Tphs0s-Tphs1s)*500))];
    vmq_c0 = (4.7773e-7.*Ptop.^0.6 + 1e-11.*Ptop) .* exp(2565*(1./(TT+273.15)-1./(cal.perT+273.15))); % Katz et al., 2003; Moore et al., 1998
    vmq_c1 = (3.5494e-3.*Ptop.^0.5 + 9.623e-8.*Ptop - 1.5223e-11.*Ptop.^1.5)./(TT+273.15) + 1.2436e-14.*Ptop.^1.5; % Liu et al., 2015
    vmq0   = (1-cc).*vmq_c0 + cc.*vmq_c1;
    perTs  = cal.perT -cal.dTH2O(2).*vmq0(round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*500)).^0.75;
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
semilogy(T,ones(size(T)).*etax0,'k',T,min(etam,eta),'r',T,ones(size(T)).*etaf0,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
semilogy(T,eta,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','mixture','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Viscosity model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Viscosity [log$_{10}$ Pas]','Interpreter','latex','FontSize',15)

% plot phase segregation speeds
figure(7); if ~holdfig; clf; end
semilogy(T,max(1e-15,abs(wx)).*3600,'k',T,max(1e-15,abs(wm)).*3600,'r',T,max(1e-15,abs(wf)).*3600,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase segregation model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Segregation flux [m/hr]','Interpreter','latex','FontSize',15)

% plot phase segregation speeds
figure(8); if ~holdfig; clf; end
subplot(2,1,1)
sgtitle('Phase Oxide Compositions','Interpreter','latex','FontSize',18)
plot(T,cm_oxds,'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend(cal.oxdStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
plot(T,cx_oxds,'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

% plot simplified mineral assemblage
figure(9); clf;
patch([T;flipud(T)],[zeros(size(T));flipud(olv)],[0.6,0.8,0.5],'LineWidth',2); hold on; box on; axis tight;
patch([T;flipud(T)],[olv;flipud(olv+pxn)],[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[olv+pxn;flipud(olv+pxn+plg)],[0.9,0.9,0.9],'LineWidth',2);
patch([T;flipud(T)],[olv+pxn+plg;flipud(olv+pxn+plg+qtz)],[0.9,0.7,1.0],'LineWidth',2);
patch([T;flipud(T)],[olv+pxn+plg+qtz;flipud(olv+pxn+plg+qtz+kfs)],[1.0,0.8,0.7],'LineWidth',2);
legend('~olivine','~pyroxene','~plagioclase','~quartz','~k-feldspar','Interpreter','latex','FontSize',15,'box','off','location','best')
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
    name = ['../out/',runID,'/',runID,'_densities'];
    print(figure(5),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_viscosity'];
    print(figure(6),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_segr_speed'];
    print(figure(7),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_minerals'];
    print(figure(7),name,'-dpng','-r300','-opengl');
end
