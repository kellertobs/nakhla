% calibrate phase diagram
clear;
addpath('../src');

% hold figures or plot new ones,choose line style
holdfig   = 0;    % set to 1 to hold figures, to 0 for new figures
linestyle = '-';

% set parameters
cphs0    =  0.35;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.70;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  4.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.55;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.50;                % peritectic solidus  composition [wt SiO2]
perT     =  1050;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  1300;                % solidus shift from water content [degC/wt^0.75]

% calculate T-dependence of melting model
T = linspace(700,1400,1e4);    % temperature [degC]
P = 1e8    .*ones(size(T));    % pressure [Pa]
c = 0.48   .*ones(size(T));    % major component [wt SiO2]
v = 0.01   .*ones(size(T));    % volatile component [wt H2O]

xq = ones(size(T)).*0.5;
fq = ones(size(T)).*0.0;

% equilibrium phase fractions and compositions
[xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T,c,v,P,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
mq = 1-fq-xq;  

% simplified mineral assemblage
ol = max(0,min(1,(cxq-perCx)./(cphs0-perCx))).*(cxq< perCx);
px = max(0,min(1,(cxq-cphs0)./(perCx-cphs0))).*(cxq< perCx) + max(0,min(1,(cxq-cphs1)./(perCx-cphs1))).*(cxq>=perCx);
fs = max(0,min(1,(cxq-perCx)./(cphs1-perCx))).*(cxq>=perCx);

% block out values below solidus and above liquidus
cxq(xq<1e-12 | mq<1e-12) = nan;
cmq(xq<1e-12 | mq<1e-12) = nan;
vfq(xq<1e-12 | mq<1e-12) = nan;
vmq(xq<1e-12 | mq<1e-12) = nan;
ol (xq<1e-12 | mq<1e-12) = nan;
px (xq<1e-12 | mq<1e-12) = nan;
fs (xq<1e-12 | mq<1e-12) = nan;
fq (xq<1e-12 | mq<1e-12) = nan;
xq (xq<1e-12 | mq<1e-12) = nan;
mq (mq>1-1e-12 | mq<1e-12) = nan;

% plot phase diagram
figure(1); if ~holdfig; clf; end
TT = linspace(Tphs0,Tphs1,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,(perT-Tphs0)./(Tphs1-Tphs0)*1e3),linspace((perCx+perCm)/2,cphs0,(perT-Tphs1)./(Tphs0-Tphs1)*1e3)];
[~,CCx,CCm,FF,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,0*TT,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
plot(CCx,TT,'k','LineStyle',linestyle,'LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k','LineStyle',linestyle,'LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase Diagram','Interpreter','latex','FontSize',22)
xlabel('Major component [wt SiO$_2$]','Interpreter','latex','FontSize',18)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)

% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,xq,'k',T,mq,'r',T,fq,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','Interpreter','latex','FontSize',18,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Melting model','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Phase fractions [wt]','Interpreter','latex','FontSize',18)

% plot major phase compositions
figure(3); if ~holdfig; clf; end
plot(T,cxq,'b',T,cmq,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',18,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase compositions','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Major component [wt SiO$_2$]','Interpreter','latex','FontSize',18)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vmq,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('melt','Interpreter','latex','FontSize',18,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase compositions','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Volatile component [wt H$_2$O]','Interpreter','latex','FontSize',18)

% plot simplified mineral assemblage
figure(5); if ~holdfig; clf; end
plot(T,ol,'k',T,px,'r',T,fs,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('olivine','pyroxene','feldspar','Interpreter','latex','FontSize',18,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Mineral assemblage','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Mineral fraction [wt]','Interpreter','latex','FontSize',18)
