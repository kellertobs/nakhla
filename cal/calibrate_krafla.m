% calibrate phase diagram
clear; addpath('../src'); close all;

% calibration run options
runID     = 'test';     % run ID for output files; [system name_wt.% SiO2_wt.% H2O]
holdfig   = 0;                           % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                         % set line style for plots
save_plot = 0;                           % turn on (1) to save output file in /out directory

% set model rheology parameters
etam0    =  100;                 % melt viscosity [Pas]
etaf0    =  0.1;                 % fluid viscosity [Pas]
etax0    =  1e15;                % crystal viscosity [Pas]
Fmc      =  1e+4;                % major component weakening factor of melt viscosity [1]
Fmv      =  0.5;                 % volatile component weakening factor of melt viscosity [1]
Em       =  150e3;               % activation energy melt viscosity [J/mol]
AA       = [ 0.60, 0.25, 0.30; 0.20, 0.20, 0.20; 0.20, 0.20, 0.20; ];  % permission slopes
BB       = [ 0.30, 0.15, 0.55; 0.48, 0.02, 0.50; 0.80, 0.08, 0.12; ];  % permission step locations
CC       = [ 0.20, 0.20, 0.20; 0.60, 0.60, 0.12; 0.20, 0.25, 0.50; ];  % permission step widths

% set model buoyancy parameters
rhom0    =  2860;                % melt phase ref. density [kg/m3] (at T0,cphs0,Ptop)
rhox0    =  3300;                % crystal phase ref. density [kg/m3] (at T0,cphs0,Ptop)
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

% Load MELTS tables in csv format
krafla_mafic_liquid = load('./krafla_mafic/Liquid_comp_tbl.csv');   % output file of MELTS for Krafla mafic comp
krafla_mafic_solid  = load('./krafla_mafic/Solid_comp_tbl.csv');   % output file of MELTS for Krafla mafic comp
krafla_mafic_bulk   = load('./krafla_mafic/Bulk_comp_tbl.csv');   % output file of MELTS for Krafla mafic comp

% read in P,T,c,v variables from MELTS tables
P_mafic = krafla_mafic_bulk(:,1).'.*1e5;
T_mafic = krafla_mafic_bulk(:,2).';
c_mafic = krafla_mafic_bulk(:,4).'./100;
v_mafic = 0.*c_mafic+0.0000;

% read in phase fractions and compositions from MELTS tables
m_mafic  = krafla_mafic_liquid(:,3).'./krafla_mafic_bulk(:,3).';
x_mafic  = krafla_mafic_solid(:,3).'./krafla_mafic_bulk(:,3).';
cm_mafic = krafla_mafic_liquid(:,4).'./100;
cx_mafic = krafla_mafic_solid(:,4).'./100;
% cx_mafic(T_mafic>=1100) = cx_mafic(T_mafic>=1100) - 0.03;
% cm_mafic(T_mafic>=1100) = cm_mafic(T_mafic>=1100) + 0.03.*cx_mafic(T_mafic>=1100)./cm_mafic(T_mafic>=1100);

% Load MELTS tables in csv format
krafla_interm_liquid = load('./krafla_interm/Liquid_comp_tbl.csv');   % output file of MELTS for Krafla interm comp
krafla_interm_solid  = load('./krafla_interm/Solid_comp_tbl.csv');   % output file of MELTS for Krafla interm comp
krafla_interm_bulk   = load('./krafla_interm/Bulk_comp_tbl.csv');   % output file of MELTS for Krafla interm comp

% read in P,T,c,v variables from MELTS tables
P_interm = krafla_interm_bulk(:,1).'.*1e5;
T_interm = krafla_interm_bulk(:,2).';
c_interm = krafla_interm_bulk(:,4).'./100;
v_interm = 0.*c_interm+0.0000;

% read in phase fractions and compositions from MELTS tables
m_interm  = krafla_interm_liquid(:,3).'./krafla_interm_bulk(:,3).';
x_interm  = krafla_interm_solid(:,3).'./krafla_interm_bulk(:,3).';
cm_interm = krafla_interm_liquid(:,4).'./100;
cx_interm = krafla_interm_solid(:,4).'./100;
% cx_interm(T_interm>=940) = cx_interm(T_interm>=940) + 0.03;
% cm_interm(T_interm>=940) = cm_interm(T_interm>=940) - 0.03.*cx_interm(T_interm>=940)./cm_interm(T_interm>=940);

% Load MELTS tables in csv format
krafla_ultram_liquid = load('./krafla_ultram/Liquid_comp_tbl.csv');   % output file of MELTS for Krafla ultram comp
krafla_ultram_solid  = load('./krafla_ultram/Solid_comp_tbl.csv');   % output file of MELTS for Krafla ultram comp
krafla_ultram_bulk   = load('./krafla_ultram/Bulk_comp_tbl.csv');   % output file of MELTS for Krafla ultram comp

% read in P,T,c,v variables from MELTS tables
P_ultram = krafla_ultram_bulk(:,1).'.*1e5;
T_ultram = krafla_ultram_bulk(:,2).';
c_ultram = krafla_ultram_bulk(:,4).'./100;
v_ultram = 0.*c_ultram+0.0000;

% read in phase fractions and compositions from MELTS tables
m_ultram  = krafla_ultram_liquid(:,3).'./krafla_ultram_bulk(:,3).';
x_ultram  = krafla_ultram_solid(:,3).'./krafla_ultram_bulk(:,3).';
cm_ultram = krafla_ultram_liquid(:,4).'./100;
cx_ultram = krafla_ultram_solid(:,4).'./100;

cphs0_bst  =  0.428;               % phase diagram lower bound composition [wt SiO2]
cphs1_bst  =  0.783;               % phase diagram upper bound composition [wt SiO2]
Tphs0_bst  =  843.7;               % phase diagram lower bound temperature [degC]
Tphs1_bst  =  1867;                % phase diagram upper bound temperature [degC]
PhDg_bst   =  [7.0,4.2,1.0,0.93];  % Phase diagram curvature factor (> 1)
perCm_bst  =  0.517;               % peritectic liquidus composition [wt SiO2]
perCx_bst  =  0.475;               % peritectic solidus  composition [wt SiO2]
perT_bst   =  1147;                % peritectic temperature [degC]
clap       =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O      =  [1200,1000,100];     % solidus shift from water content [degC/wt^0.75]
beta       =  0.9;                 % iterative lag parameter phase diagram [1]

misfit  = 1e3; tol = 0.4;
bestfit = misfit;

time0 = tic;
timelimit = 60*60*0.25; % Limit 15min %%change back to 5min after completion of misfit

while misfit > tol
    % set phase diagram parameters
    cphs0    =  cphs0_bst * (1 + randn(1)*1e-3);          % phase diagram lower bound composition [wt SiO2]
    cphs1    =  cphs1_bst * (1 + randn(1)*1e-3);          % phase diagram upper bound composition [wt SiO2]
    Tphs0    =  Tphs0_bst * (1 + randn(1)*1e-3);          % phase diagram lower bound temperature [degC]
    Tphs1    =  Tphs1_bst * (1 + randn(1)*1e-3);          % phase diagram upper bound temperature [degC]
    PhDg     =  PhDg_bst .* (1 + randn(1,4)*1e-3);        % Phase diagram curvature factor (> 1)
    perCm    =  perCm_bst * (1 + randn(1)*3e-3);          % peritectic liquidus composition [wt SiO2]
    perCx    =  perCx_bst * (1 + randn(1)*3e-3);          % peritectic solidus  composition [wt SiO2]
    perT     =  perT_bst  * (1 + randn(1)*3e-3);          % peritectic temperature [degC]  
    
    % equilibrium phase fractions and compositions - Krafla mafic
    [xq_mafic,cxq_mafic,cmq_mafic,fq_mafic,~,~]  =  equilibrium(ones(size(T_mafic)).*0.5,ones(size(T_mafic)).*0.0, ...
        T_mafic, c_mafic, v_mafic, P_mafic, ...
        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
    mq_mafic = 1-fq_mafic-xq_mafic;

    % equilibrium phase fractions and compositions - Krafla intermediate
    [xq_interm,cxq_interm,cmq_interm,fq_interm,~,~]  =  equilibrium(ones(size(T_interm)).*0.5,ones(size(T_interm)).*0.0, ...
        T_interm, c_interm, v_interm, P_interm, ...
        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
    mq_interm = 1-fq_interm-xq_interm;
    
    % equilibrium phase fractions and compositions - Krafla ultramediate
    [xq_ultram,cxq_ultram,cmq_ultram,fq_ultram,~,~]  =  equilibrium(ones(size(T_ultram)).*0.5,ones(size(T_ultram)).*0.0, ...
        T_ultram, c_ultram, v_ultram, P_ultram, ...
        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
    mq_ultram = 1-fq_ultram-xq_ultram;
    
    misfit =          norm( m_mafic.*linspace(1.15,1,length(m_mafic))- mq_mafic,2)./norm( mq_mafic,2) ...
                    + norm( x_mafic- xq_mafic,2)./norm( xq_mafic,2) ...                    
                    + norm(cm_mafic-cmq_mafic,2)./norm(cmq_mafic,2) ...
                    + norm(cx_mafic-cxq_mafic,2)./norm(cxq_mafic,2)*0.8;
    misfit = misfit + norm( m_interm- mq_interm,2)./norm( mq_interm,2) ...
                    + norm( x_interm- xq_interm,2)./norm( xq_interm,2) ...                    
                    + norm(cm_interm-cmq_interm,2)./norm(cmq_interm,2) ...
                    + norm(cx_interm-cxq_interm,2)./norm(cxq_interm,2)*0.8;
   misfit = misfit + norm(  m_ultram- mq_ultram,2)./norm( mq_ultram,2) ...
                    + norm( x_ultram- xq_ultram,2)./norm( xq_ultram,2) ...                    
                    + norm(cm_ultram-cmq_ultram,2)./norm(cmq_ultram,2) ...
                    + norm(cx_ultram-cxq_ultram,2)./norm(cxq_ultram,2); 

    if misfit < 1.0005*bestfit         % if misfit < bestfit set bestfit values as misfit
        fprintf(1,'   misfit = %1.4e \n',misfit);
        bestfit = misfit;       %set bestfit values of parameter to model parameters
        cphs0_bst = cphs0;
        cphs1_bst = cphs1;
        Tphs0_bst = Tphs0;
        Tphs1_bst = Tphs1;
        PhDg_bst  = PhDg;
        perCm_bst = perCm;
        perCx_bst = perCx;
        perT_bst  = perT;

        % plot phase fractions
        figure(1); if ~holdfig; clf; end
        plot(T_mafic,xq_mafic.*100,'k-' ,T_mafic,mq_mafic.*100,'r-' ,T_mafic,fq_mafic.*1000,'b','LineWidth',2); hold on; box on; axis tight;
        plot(T_interm,xq_interm.*100,'k-' ,T_interm,mq_interm.*100,'r-' ,T_interm,fq_interm.*1000,'b','LineWidth',2);
        plot(T_ultram,xq_ultram.*100,'k-' ,T_ultram,mq_ultram.*100,'r-' ,T_ultram,fq_ultram.*1000,'b','LineWidth',2);
        plot(T_mafic,x_mafic*100,'kd',T_mafic,m_mafic*100,'rd','LineWidth',2,'MarkerSize',5);
        plot(T_interm,x_interm*100,'kv',T_interm,m_interm*100,'rv','LineWidth',2,'MarkerSize',5);
        plot(T_ultram,x_ultram*100,'kv',T_ultram,m_ultram*100,'rv','LineWidth',2,'MarkerSize',5);
        set(gca,'TickLabelInterpreter','latex','FontSize',13)
        title('Melting model','Interpreter','latex','FontSize',18)
        xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
        ylabel('Phase fractions [wt\%]','Interpreter','latex','FontSize',15)

        % plot major phase compositions
        figure(2); if ~holdfig; clf; end
        plot(T_mafic,cxq_mafic.*100,'b-' ,T_mafic,cmq_mafic.*100,'r-' ,'LineWidth',2); hold on; box on; axis tight;
        plot(T_interm,cxq_interm.*100,'b-' ,T_interm,cmq_interm.*100,'r-' ,'LineWidth',2);
        plot(T_ultram,cxq_ultram.*100,'b-' ,T_ultram,cmq_ultram.*100,'r-' ,'LineWidth',2);
        plot(T_mafic,cx_mafic.*100,'bd' ,T_mafic,cm_mafic.*100,'rd' ,'LineWidth',2);
        plot(T_interm,cx_interm.*100,'bv' ,T_interm,cm_interm.*100,'rv' ,'LineWidth',2);
        plot(T_ultram,cx_ultram.*100,'bv' ,T_ultram,cm_ultram.*100,'rv' ,'LineWidth',2);
        set(gca,'TickLabelInterpreter','latex','FontSize',13)
        title('Phase compositions','Interpreter','latex','FontSize',18)
        xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
        ylabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
        
        drawnow
    end

    if toc(time0)>timelimit  % if loop takes too long, break after time limit
        break
    end
end

cphs0 = cphs0_bst
cphs1 = cphs1_bst
Tphs0 = Tphs0_bst
Tphs1 = Tphs1_bst
PhDg  = PhDg_bst
perCm = perCm_bst
perCx = perCx_bst
perT  = perT_bst

% set ranges for control variables T, c, v, P - Krafla
T = linspace(min(T_mafic)-100,max(T_mafic)+100,1e3);  % temperature range [degC]
c = 0.514 * ones(size(T));                       % major component range [wt SiO2]
v = 0.000 * ones(size(T));                       % volatile component range [wt H2O]
P = 100e6 * ones(size(T));                       % pressure range [Pa]

    % equilibrium phase fractions and compositions - Krafla
    [xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(ones(size(T)).*0.5,ones(size(T)).*0.0, ...
        T, c, v, P, ...
        Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
    mq = 1-fq-xq;
    
% plot phase diagram
figure(3); if ~holdfig; clf; end
vv = (4.8e-5.*mean(P(:)).^0.6 + 1e-9.*mean(P(:)))./100;
TT = linspace(Tphs0+mean(P(:))*clap,Tphs1+mean(P(:))*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perT-Tphs0)./(Tphs1-Tphs0)*1e3)),linspace((perCx+perCm)/2,cphs0,floor((perT-Tphs1)./(Tphs0-Tphs1)*1e3))];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,mean(P(:))*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
Tphs0s = Tphs0-dTH2O(1)*vv^0.75;
Tphs1s = Tphs1-dTH2O(3)*vv^0.75;
perTs  = perT-dTH2O(2)*vv^0.75;
TT = linspace(Tphs0s+mean(P(:))*clap,Tphs1s+mean(P(:))*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,round((perTs-Tphs0s)./(Tphs1s-Tphs0s)*1e3)),linspace((perCx+perCm)/2,cphs0,round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*1e3))];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,vv*ones(size(TT)),mean(P(:))*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
plot(cxq.*100,T-(P-mean(P(:))).*clap,'b','LineStyle',linestyle,'LineWidth',2);
plot(cmq.*100,T-(P-mean(P(:))).*clap,'r','LineStyle',linestyle,'LineWidth',2);
plot(c./(1-fq+1e-16).*100,T-(P-mean(P(:))).*clap,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2);
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vfq/10.*100,'b',T,vmq.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('fluid $/10$','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Volatile component [wt\% H$_2$O]','Interpreter','latex','FontSize',15)



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
    
   % save ('name','cphs0','cphs1','Tphs0','Tphs1','PhDg','perCm','perCx','perT','clap','dTH2O' 'beta', 'bestfit', '-ascii');
end

% for i = bestfit         %not yet done to generate an output file
%     X = []
%    fid = sprintf('output_matlab_%d.txt', i);
%     
%     fprintf(fid,'%4.4f\n',[cphs0 cphs1 Tphs0 Tphs1 PhDg perCm perCx perT]);
%     fclose(fid, 'all');
% end