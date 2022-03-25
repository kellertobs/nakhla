% calibrate phase diagram
clear; addpath('../src');

% calibration run options
runID     = 'test';              % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set phase diagram parameters
cphs0    =  0.35;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.75;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  5.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.51;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.47;                % peritectic solidus  composition [wt SiO2]
perT     =  1100;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  [1200,1000,200];     % solidus shift from water content [degC/wt^0.75
beta     =  0.80;                % iterative lag parameter phase diagram [1]

% set model rheology parameters
etam0    =  300;                 % melt viscosity [Pas]
etaf0    =  0.1;                 % fluid viscosity [Pas]
etax0    =  1e15;                % crystal viscosity [Pas]
Fmc      =  1e+4;                % major component weakening factor of melt viscosity [1]
Fmv      =  0.4;                 % volatile component weakening factor of melt viscosity [1]
Em       =  150e3;               % activation energy melt viscosity [J/mol]
AA       = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
BB       = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
CC       = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% set model buoyancy parameters
rhom0    =  2700;                % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhox0    =  3100;                % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhof0    =  500;                 % bubble phase ref. density [kg/m3] (at T0,cphs0,Ptop)
aTm      =  3e-5;                % melt thermal expansivity [1/K]
aTx      =  1e-5;                % xtal thermal expansivity [1/K]
aTf      =  1e-4;                % mvp  thermal expansivity [1/K]
gCm      =  0.5;                 % melt compositional expansion [1/wt]
gCx      =  0.6;                 % xtal compositional expansion [1/wt]
bPf      =  1e-8;                % mvp compressibility [1/Pa]
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
dm       =  1e-4;                % melt film size [m]
g0       =  10.;                 % gravity [m/s2]

% set ranges for control variables T, c, v, P
T = linspace(500,1600,1e3);    % temperature range [degC]
c = linspace(0.49,0.49,1e3);   % major component range [wt SiO2]
v = linspace(0.02,0.02,1e3);   % volatile component range [wt H2O]
P = linspace(100,100,1e3)*1e6; % pressure range [Pa]

% equilibrium phase fractions and compositions
[xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(ones(size(T)).*0.5,ones(size(T)).*0.0, ...
                                        T, c, v, P, ...
                                        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
mq = 1-fq-xq;  

% phase densities and rheology
rhom = rhom0 .* (1 - aTm.*(T-perT) - gCm.*(cmq-(perCx+perCm)/2 ));
rhox = rhox0 .* (1 - aTx.*(T-perT) - gCx.*(cxq-(perCx+perCm)/2 ));
rhof = rhof0 .* (1 - aTf.*(T-perT) + bPf.*(P  -min(P)));

rho    = 1./(mq./rhom + xq./rhox + fq./rhof);
chi    = xq.*rho./rhox;
phi    = fq.*rho./rhof;
mu     = mq.*rho./rhom;

etam  = etam0 .* exp(Em./(8.3145.*(T+273.15))-Em./(8.3145.*(perT+273.15))) ...
              .* Fmc.^((cmq-(perCx+perCm)/2)./(cphs1-cphs0)) .* Fmv.^(vmq./0.01); % T-c-v-dep. melt viscosity
etaf  = etaf0.* ones(size(fq));                                  % constant volatile fluid viscosity
etax  = etax0.* ones(size(xq));                                  % constant volatile fluid viscosity

% get permission weights
kv = permute(cat(3,etax,etam,etaf),[3,1,2]);
Mv = squeeze(permute(repmat(kv,1,1,1,3),[4,1,2,3])./permute(repmat(kv,1,1,1,3),[1,4,2,3]));
kv = squeeze(kv);

ff = squeeze(permute(cat(3,chi,mu,phi),[3,1,2]));
FF = permute(repmat(ff,1,1,1,3),[4,1,2,3]);
Sf = (FF./BB).^(1./CC);  Sf = Sf./sum(Sf,2);
Xf = sum(AA.*Sf,2).*FF + (1-sum(AA.*Sf,2)).*Sf;

% get momentum and volume permissions
thtv = squeeze(prod(Mv.^Xf,2));

% get momentum and volume flux and transfer coefficients
Kv =    ff .*kv.*thtv;
Cv = (1-ff)./[dx;dm;df].^2.*Kv;

eta    = squeeze(sum(Kv,1));
Ksgr_x = max(1e-24,min(1e-6,chi./squeeze(Cv(1,:,:))));
Ksgr_m = max(1e-24,min(1e-6,mu ./squeeze(Cv(2,:,:))));
Ksgr_f = max(1e-24,min(1e-6,phi./squeeze(Cv(3,:,:))));

wx = chi.*Ksgr_x .* (rhox-rho)*g0; % crystal segregation speed
wf = phi.*Ksgr_f .* (rhof-rho)*g0; % fluid segregation speed
wm = mu .*Ksgr_m .* (rhom-rho)*g0; % melt segregation speed

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
ind = xq<1e-12 | mq<1e-12;
cxq (ind) = [];
cmq (ind) = [];
vfq (ind) = [];
vmq (ind) = [];
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
wx  (ind) = [];
wf  (ind) = [];
wm  (ind) = [];
fq  (ind) = [];
xq  (ind) = [];
mq  (ind) = [];


if ~holdfig; close all; end

% plot phase diagram
figure(1); if ~holdfig; clf;
TT = Tphs0+mean(P(:))*clap:1:Tphs1+mean(P(:))*clap;
[~,CCx,~,~,~,~] = equilibrium(0*TT,0*TT,TT,cphs1*ones(size(TT)),0*TT,mean(P(:))*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
[~,~,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cphs0*ones(size(TT)),0*TT,mean(P(:))*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
perTs  = perT;
Tphs0s = Tphs0;
Tphs1s = Tphs1;
vv = 0.10*ones(size(TT));
for i = 1:10
    perTs  = perT -dTH2O(2)*mean(vv(abs(TT-mean(P(:))*clap-perTs )<1)).^0.75;
    Tphs0s = Tphs0-dTH2O(1)*mean(vv(abs(TT-mean(P(:))*clap-Tphs0s)<1)).^0.75;
    Tphs1s = Tphs1-dTH2O(3)*mean(vv(abs(TT-mean(P(:))*clap-Tphs1s)<1)).^0.75;
    TTi = Tphs0s+mean(P(:))*clap:1:Tphs1s+mean(P(:))*clap;
    vv = interp1(TT,vv,TTi,'linear','extrap'); TT = TTi;
    cc = [linspace(cphs1,(perCx+perCm)/2,round((perTs-Tphs0s)./(Tphs1s-Tphs0s)*length(TT))),linspace((perCx+perCm)/2,cphs0,round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*length(TT)))];
    [~,CCx,CCm,~,~,vv] = equilibrium(0*TT,0*TT,TT,cc,vv,mean(P(:))*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
end
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
end
plot(cxq.*100,T-(P-mean(P(:))).*clap,'b','LineStyle',linestyle,'LineWidth',2);
plot(cmq.*100,T-(P-mean(P(:))).*clap,'r','LineStyle',linestyle,'LineWidth',2);
plot(c./(1-fq+1e-16).*100,T-(P-mean(P(:))).*clap,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,xq.*100,'k',T,mq.*100,'r',T,fq.*1000,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',15,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Melting model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Phase fractions [wt\%]','Interpreter','latex','FontSize',15)

% plot major phase compositions
figure(3); if ~holdfig; clf; end
plot(T,cxq.*100,'b',T,cmq.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vfq/10.*100,'b',T,vmq.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
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
semilogy(T,max(1e-12,abs(wx)).*3600,'k',T,max(1e-12,abs(wf)).*3600,'b',T,max(1e-12,abs(wm)).*3600,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','fluid','melt','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase segregation model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Segregation flux [m/hr]','Interpreter','latex','FontSize',15)

% plot simplified mineral assemblage
figure(8); clf;
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
