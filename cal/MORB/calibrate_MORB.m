% calibrate phase diagram
addpath('../../src');
TINY = 1e-16;


% calibration run options
runID     = 'test';              % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% load MAGEMin results
load MORB_fract_dry_P150_out.mat
MORB = OUT;
load KLB1_fract_dry_P150_out.mat
KLB1 = OUT;

% set phase diagram parameters
calID    =  'morb';             % phase diagram calibration

% set model buoyancy parameters
d0       =  1e-3;                % crystal size [m]
g0       =  10.;                 % gravity [m/s2]

% set ranges for control variables T, c, v, P
T = linspace(1100,1850,1e3).';    % temperature range [degC]
c = linspace(0.4486,0.4486,1e3).';   % major component range [wt SiO2]
v = linspace(0.0,0.0,1e3).';   % volatile component range [wt H2O]
P = linspace(150,150,1e3).'*1e6; % pressure range [Pa]

% equilibrium phase fractions and compositions
run(['./cal_',calID]);  % load melt model calibration
[x,cx,cm,f,vf,vm]  =  equilibrium(ones(size(T)).*0.5,v./10,T,c,v,P,cal,TINY);
m = 1-f-x;  

Nz = length(T); Nx = 1; Ptop = min(P); Pt = P; etareg = 1; calibrt = 1; T = T+273.15;
update;
T = T-273.15;

% block out values below solidus and above liquidus
ind = x<1e-9 | m<1e-9;
c_oxd  = squeeze(c_oxd); c_oxd(ind,:) = [];
cm_oxd = squeeze(cm_oxd); cm_oxd(ind,:) = [];
cx_oxd = squeeze(cx_oxd); cx_oxd(ind,:) = [];
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
f  (ind) = [];
x  (ind) = [];
m  (ind) = [];
Ptop = min(P); Pt = P;

if ~holdfig; close all; end

% plot phase diagram
figure(1); if ~holdfig; clf; end
TT = linspace(cal.Tphs0+Ptop*cal.clap,cal.Tphs1+Ptop*cal.clap,200);
cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2, ceil((cal.perT-cal.Tphs0)./(cal.Tphs1-cal.Tphs0)*200)), ...
      linspace((cal.perCx+cal.perCm)/2,cal.cphs0,floor((cal.perT-cal.Tphs1)./(cal.Tphs0-cal.Tphs1)*200))];
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
    TT = linspace(Tphs0s+Ptop*cal.clap,Tphs1s+Ptop*cal.clap,200);
    cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2, ceil((perTs-Tphs0s)./(Tphs1s-Tphs0s)*200)), ...
          linspace((cal.perCx+cal.perCm)/2,cal.cphs0,floor((perTs-Tphs1s)./(Tphs0s-Tphs1s)*200))];
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

plot(MORB.OxideFractions.sol(:,1),MORB.T,'bo')
plot(MORB.OxideFractions.liq(:,1),MORB.T,'ro')
plot(KLB1.OxideFractions.sol(:,1),KLB1.T,'bo')
plot(KLB1.OxideFractions.liq(:,1),KLB1.T,'ro')

axis([0.4 0.8,650 1900])
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,x.*100,'k',T,m.*100,'r',T,f.*1000,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',15,'box','off','location','east')

% plot(MORB.T,MORB.PhaseFractions.sol_wt*100,'ko')
% plot(MORB.T,MORB.PhaseFractions.liq_wt*100,'ro')
plot(KLB1.T,KLB1.PhaseFractions.sol_wt*100,'ko')
plot(KLB1.T,KLB1.PhaseFractions.liq_wt*100,'ro')

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

% plot oxide compositions
figure(7); if ~holdfig; clf; end
subplot(2,1,1)
sgtitle('Phase Oxide Compositions','Interpreter','latex','FontSize',18)
plot(T,cm_oxd,'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend(cal.oxdStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
plot(T,cx_oxd,'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

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
