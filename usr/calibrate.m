% calibrate phase diagram
clear; addpath('../src');

% calibration run options
runID     = 'test';    % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set parameters
cphs0    =  0.36;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.72;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  760;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  4.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.52;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.48;                % peritectic solidus  composition [wt SiO2]
perT     =  1100;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  1300;                % solidus shift from water content [degC/wt^0.75]

% calculate T-dependence of melting model
T = linspace(600,1600,1e3);    % temperature range [degC]
c = linspace(0.50,0.50,1e3);   % major component range [wt SiO2]
v = linspace(0.02,0.02,1e3);   % volatile component range [wt H2O]
P = linspace(200,200,1e3)*1e6; % pressure range [Pa]

% equilibrium phase fractions and compositions
[xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(ones(size(T)).*0.5,ones(size(T)).*0.0, ...
                                        T, c, v, P, ...
                                        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
mq = 1-fq-xq;  

% simplified mineral assemblage (ol+px+fs+qz)
olv0 = perCx;  olv1 = cphs0; 
pxn0 = cphs0; pxn1 = perCx;  pxn2 = cphs1;
plg0 = (cphs0+perCx)/2;  plg1 = perCm;  plg2 = (cphs1+1)/2;
kfs0 = perCm;  kfs1 = cphs1;  kfs2 = 1;
qtz0 = perCx;  qtz1 = 1;

olv = max(0,min(1,(cxq-olv0)./(olv1-olv0))).*(cxq>=olv1 & cxq<=olv0);
pxn = max(0,min(1,(cxq-pxn0)./(pxn1-pxn0))).*(cxq>=pxn0 & cxq<=pxn1) + max(0,min(1,(cxq-pxn2)./(pxn1-pxn2))).*(cxq>=pxn1 & cxq<=pxn2);
plg = max(0,min(1,(cxq-plg0)./(plg1-plg0))).*(cxq>=plg0 & cxq<=plg1) + max(0,min(1,(cxq-plg2)./(plg1-plg2))).*(cxq>=plg1 & cxq<=plg2);
kfs = max(0,min(1,(cxq-kfs0)./(kfs1-kfs0))).*(cxq>=kfs0 & cxq<=kfs1) + max(0,min(1,(cxq-kfs2)./(kfs1-kfs2))).*(cxq>=kfs1 & cxq<=kfs2);
qtz = max(0,min(1,(cxq-qtz0)./(qtz1-qtz0))).*(cxq>=qtz0 & cxq<=qtz1);
all = olv+pxn+plg+kfs+qtz;
olv = olv./all;
pxn = pxn./all;
plg = plg./all;
kfs = kfs./all;
qtz = qtz./all;

% block out values below solidus and above liquidus
cxq(xq<1e-12 | mq<1e-12) = [];
cmq(xq<1e-12 | mq<1e-12) = [];
vfq(xq<1e-12 | mq<1e-12) = [];
vmq(xq<1e-12 | mq<1e-12) = [];
olv(xq<1e-12 | mq<1e-12) = [];
pxn(xq<1e-12 | mq<1e-12) = [];
plg(xq<1e-12 | mq<1e-12) = [];
kfs(xq<1e-12 | mq<1e-12) = [];
qtz(xq<1e-12 | mq<1e-12) = [];
c  (xq<1e-12 | mq<1e-12) = [];
v  (xq<1e-12 | mq<1e-12) = [];
T  (xq<1e-12 | mq<1e-12) = [];
P  (xq<1e-12 | mq<1e-12) = [];
fq (xq<1e-12 | mq<1e-12) = [];
xq (xq<1e-12 | mq<1e-12) = [];
mq (mq>1-1e-12 | mq<1e-12) = [];


if ~holdfig; close all; end

% plot phase diagram
figure(1); if ~holdfig; clf; end
TT = linspace(Tphs0,Tphs1,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perT-Tphs0)./(Tphs1-Tphs0)*1e3)),linspace((perCx+perCm)/2,cphs0,floor((perT-Tphs1)./(Tphs0-Tphs1)*1e3))];
[~,CCx,CCm,FF,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,0*TT,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
plot(CCx,TT,'k','LineStyle',linestyle,'LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k','LineStyle',linestyle,'LineWidth',2);
plot(cxq,T-P.*clap + dTH2O*vmq.^0.75,'b','LineStyle',linestyle,'LineWidth',2);
plot(cmq,T-P.*clap + dTH2O*vmq.^0.75,'r','LineStyle',linestyle,'LineWidth',2);
plot(c./(1-fq+1e-16),T-P.*clap + dTH2O*vmq.^0.75,'k','LineStyle',':','LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,xq,'k',T,mq,'r',T,fq.*10,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',15,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Melting model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Phase fractions [wt]','Interpreter','latex','FontSize',15)

% plot major phase compositions
figure(3); if ~holdfig; clf; end
plot(T,cxq,'b',T,cmq,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Major component [wt SiO$_2$]','Interpreter','latex','FontSize',15)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vfq/10,'b',T,vmq,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('fluid $/10$','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Volatile component [wt H$_2$O]','Interpreter','latex','FontSize',15)

% plot simplified mineral assemblage
figure(5); if ~holdfig; clf; end
patch([T,fliplr(T)],[zeros(size(T)),fliplr(olv)],[0.6,0.8,0.5],'LineWidth',2); hold on; box on; axis tight;
patch([T,fliplr(T)],[olv,fliplr(olv+pxn)],[0.6,0.6,0.6],'LineWidth',2);
patch([T,fliplr(T)],[olv+pxn,fliplr(olv+pxn+plg)],[0.9,0.9,0.9],'LineWidth',2);
patch([T,fliplr(T)],[olv+pxn+plg,fliplr(olv+pxn+plg+qtz)],[0.9,0.7,1.0],'LineWidth',2);
patch([T,fliplr(T)],[olv+pxn+plg+qtz,fliplr(olv+pxn+plg+qtz+kfs)],[1.0,0.8,0.7],'LineWidth',2);
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
    name = ['../out/',runID,'/',runID,'_calibrate_phase_diagram'];
    print(figure(1),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_calibrate_melt_model'];
    print(figure(2),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_calibrate_maj_compnt'];
    print(figure(3),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_calibrate_vol_compnt'];
    print(figure(4),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_calibrate_minerals'];
    print(figure(5),name,'-dpng','-r300','-opengl');
end
