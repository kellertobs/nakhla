% calibrate phase diagram
clear; addpath('../src');

% calibration run options
runID     = 'Krafla_074_003';    % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 1;                   % turn on (1) to save output file in /out directory

% set parameters
cphs0    =  0.36;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.72;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  4.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.54;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.50;                % peritectic solidus  composition [wt SiO2]
perT     =  1100;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  1300;                % solidus shift from water content [degC/wt^0.75]

% calculate T-dependence of melting model
T = linspace(500,1600,1e3);    % temperature [degC]
c = 0.55   .*ones(size(T));    % major component [wt SiO2]
v = 0.02   .*ones(size(T));    % volatile component [wt H2O]
P = 125e6  .*ones(size(T));    % pressure [Pa]

% equilibrium phase fractions and compositions
[xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(ones(size(T)).*0.5,ones(size(T)).*0.0, ...
                                        T, c, v, P, ...
                                        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
mq = 1-fq-xq;  

% simplified mineral assemblage (ol+px+fs+qz)
ol0 = perCx;  ol1 = cphs0; 
px0 = (cphs0+perCx)/2; px1 = perCx;  px2 = cphs1;
fs0 = perCx;  fs1 = cphs1;  fs2 = (1+cphs1)/2;
qz0 = perCm;  qz1 = 1;

ol = max(0,min(1,(cxq-ol0)./(ol1-ol0))).*(cxq>=ol1 & cxq<=ol0);
px = max(0,min(1,(cxq-px0)./(px1-px0))).*(cxq>=px0 & cxq<=px1) + max(0,min(1,(cxq-px2)./(px1-px2))).*(cxq>=px1 & cxq<=px2);
fs = max(0,min(1,(cxq-fs0)./(fs1-fs0))).*(cxq>=fs0 & cxq<=fs1) + max(0,min(1,(cxq-fs2)./(fs1-fs2))).*(cxq>=fs1 & cxq<=fs2);
qz = max(0,min(1,(cxq-qz0)./(qz1-qz0))).*(cxq>=qz0 & cxq<=qz1);
ol = ol./(ol+px+fs+qz);
px = px./(ol+px+fs+qz);
fs = fs./(ol+px+fs+qz);
qz = qz./(ol+px+fs+qz);

% block out values below solidus and above liquidus
cxq(xq<1e-12 | mq<1e-12) = nan;
cmq(xq<1e-12 | mq<1e-12) = nan;
vfq(xq<1e-12 | mq<1e-12) = nan;
vmq(xq<1e-12 | mq<1e-12) = nan;
ol (xq<1e-12 | mq<1e-12) = nan;
px (xq<1e-12 | mq<1e-12) = nan;
fs (xq<1e-12 | mq<1e-12) = nan;
qz (xq<1e-12 | mq<1e-12) = nan;
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
plot(cxq,T-P.*clap + dTH2O*vmq.^0.75,'b','LineStyle',linestyle,'LineWidth',2);
plot(cmq,T-P.*clap + dTH2O*vmq.^0.75,'r','LineStyle',linestyle,'LineWidth',2);
plot(c./(1-fq+1e-16),T-P.*clap + dTH2O*vmq.^0.75,'k','LineStyle',':','LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase Diagram','Interpreter','latex','FontSize',22)
xlabel('Major component [wt SiO$_2$]','Interpreter','latex','FontSize',18)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)

% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,xq,'k',T,mq,'r',T,fq.*10,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',18,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Melting model','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Phase fractions [wt]','Interpreter','latex','FontSize',18)

% plot major phase compositions
figure(3); if ~holdfig; clf; end
plot(T,cxq,'b',T,cmq,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',18,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase compositions','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Major component [wt SiO$_2$]','Interpreter','latex','FontSize',18)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vfq/10,'b',T,vmq,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('fluid $/10$','melt','Interpreter','latex','FontSize',18,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase compositions','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Volatile component [wt H$_2$O]','Interpreter','latex','FontSize',18)

% plot simplified mineral assemblage
figure(5); if ~holdfig; clf; end
plot(T,ol,'g',T,px,'k',T,fs,'b',T,qz,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('olivine','pyroxene','feldspar','quartz','Interpreter','latex','FontSize',18,'box','off','location','west')
set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Mineral assemblage','Interpreter','latex','FontSize',22)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',18)
ylabel('Mineral fraction [wt]','Interpreter','latex','FontSize',18)

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
