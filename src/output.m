
if plot_op
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',16};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',12};
    UN = {'Units','Centimeters'};
    
    % set axis and border dimensions
    axh = 8.00; axw = axh*L/D;
    ahs = 0.50; avs = 0.25;
    axb = 1.50; axt = 0.75;
    axl = 1.50; axr = 0.75;
    
    % initialize figures and axes
    fh1 = figure(1); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh1,UN{:},'Position',[3 3 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off');
    set(fh1,'Resize','off');
    ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(13) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(14) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    fh2 = figure(2); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh2,UN{:},'Position',[6 6 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off');
    set(fh2,'Resize','off');
    ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(23) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(24) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    
    fh3 = figure(3); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh3,UN{:},'Position',[9 9 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','off');
    ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(33) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(34) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    
    fh4 = figure(4); clf; colormap(ocean);
    fh = axb + 1*axh + 0*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh4,UN{:},'Position',[12 12 fw fh]);
    set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh4,'Color','w','InvertHardcopy','off');
    set(fh4,'Resize','off');
    ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(43) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
        
    if plot_cv
        fh5 = figure(5); clf; colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh5,UN{:},'Position',[15 15 fw fh]);
        set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh5,'Color','w','InvertHardcopy','off');
        set(fh5,'Resize','off');
        ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(53) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end
    
    % plot velocity-pressure solution in Fig. 1
    figure(1);
    axes(ax(11));
    imagesc(X(2:end-1),Z(2:end-1),-W(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/s]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(12));
    imagesc(X(2:end-1),Z(2:end-1), U(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [m/s]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(13));
    imagesc(X(2:end-1),Z(2:end-1), P(2:end-1,2:end-1)./1e3); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [kPa]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(14));
    imagesc(X(2:end-1),Z(2:end-1),T(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T \ [^\circ$C]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);

    
    % plot thermo-chemical solution in Fig. 2
    figure(2);
    axes(ax(21));
    imagesc(X(2:end-1),Z(2:end-1),chi(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(22));
    imagesc(X(2:end-1),Z(2:end-1),phi(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\phi$ [vol]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(23));
    imagesc(X(2:end-1),Z(2:end-1),c(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$c$ [wt\% SiO$_2$]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(24));
    imagesc(X(2:end-1),Z(2:end-1),v(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$v$ [wt\% H$_2$O]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    % plot density and rheology fields in Fig. 3
    figure(3);
    axes(ax(31));
    imagesc(X(2:end-1),Z(2:end-1),      rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\rho$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(32));
    imagesc(X(2:end-1),Z(2:end-1),log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\eta$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(33));
    imagesc(X(2:end-1),Z(2:end-1),Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_\chi/\bar{\rho}$ [1/s]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(34));
    imagesc(X(2:end-1),Z(2:end-1),Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_\phi/\bar{\rho}$ [1/s]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
%     axes(ax(33));
%     imagesc(X(2:end-1),Z(2:end-1),log10(eII(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
%     set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\varepsilon_{II}$ [log$_{10}$ 1/s]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
%     axes(ax(34));
%     imagesc(X(2:end-1),Z(2:end-1),log10(tII(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
%     set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\tau_{II}$ [log$_{10}$ Pa]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);

    figure(4);
    axes(ax(41));
    imagesc(X(2:end-1),Z(2:end-1),it(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['incomp. trace'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(42));
    imagesc(X(2:end-1),Z(2:end-1),ct(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['comp. trace'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(43));
    imagesc(X(2:end-1),Z(2:end-1),si(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['stable isotope'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    if plot_cv && iter > 0
        % plot residual fields in Fig. 4
        figure(5);
        axes(ax(51));
        imagesc(X(2:end-1),Z(2:end-1),-res_W(:      ,2:end-1)./(1e-16+norm(RR(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $W$'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
        axes(ax(52));
        imagesc(X(2:end-1),Z(2:end-1), res_U(2:end-1,:      )./(1e-16+norm(RR(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $U$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
        axes(ax(53));
        imagesc(X(2:end-1),Z(2:end-1), res_P(2:end-1,2:end-1)./(1e-16+norm(RR(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $P$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); 
    end
end

% plot phase diagram
fh6 = figure(6); clf;
TT = linspace(Tphs0,Tphs1,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,(perT-Tphs0)./(Tphs1-Tphs0)*1e3),linspace((perCx+perCm)/2,cphs0,(perT-Tphs1)./(Tphs0-Tphs1)*1e3)];
[~,CCx,CCm,FF,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,0*TT,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
plot(CCx,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k-','LineWidth',2);

% plot([perCx,cphs1],[Tphs0,Tphs0],'k-','LineWidth',1.5)
% plot([cphs0,perCx],[perT,perT],'k-','LineWidth',1.5)
% plot([perCx,perCm],[perT,perT],'k-','LineWidth',1)
% plot([perCx,perCx],[Tphs0,perT],'k-','LineWidth',1.5)

Tplt = T - Pt*clap + dTH2O*vm.^0.75;
cplt = c./(1-f);
plot(cplt(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'k.',cx(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'b.',cm(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'r.','LineWidth',2,'MarkerSize',15);

set(gca,'TickLabelInterpreter','latex','FontSize',15)
% text(perCx/2,(perT-perT/2)/2,'sol 1 + sol 2','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(perCx/2,perT+0.4*(1-perT),'sol 1 + liq','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
% text((perCx+1)/2,perT*0.4,'sol 2 + liq','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
% text((perCx+1)/2,-perT/4,'sol 2 + sol 3','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
% text(perCm-0.2,perT+0.75*(1-perT),'liq','Interpreter','latex','FontSize',18,'HorizontalAlignment','center','VerticalAlignment','middle')
title('Phase Diagram','Interpreter','latex','FontSize',22)
xlabel('Composition','Interpreter','latex','FontSize',18)
ylabel('Temperature','Interpreter','latex','FontSize',18)
        
% plot model history
fh7 = figure(7);
if iter > 0 % don't plot before solver has run first time
    subplot(4,1,1);
    Qcool = sum(sum(rhoCp(2:end-1,2:end-1).*cool(2:end-1,2:end-1).*h^3));
    plot(time./3600,Qcool./1e6,'bo','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('$Q_{h,c}$ [MW]',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(4,1,2);
    meanT = mean(mean(T(2:end-1,2:end-1)));
    plot(time./3600,meanT,'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('$T$ [$^\circ$C]',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(4,1,3);
    rmsV = sqrt(sum(sum(W(:,2:end-1).^2))+sum(sum(U(2:end-1,:).^2)))./sqrt((N-2)*(N-1));
    maxV = max(max(max(abs(W(:,2:end-1)))),max(max(abs(U(2:end-1,:)))));
    plot(time./3600,rmsV,'ko',time./3600,maxV,'ro','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('$\mathbf{V}$ [m/s]',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(4,1,4);
    meanf = mean(mean(f(2:end-1,2:end-1)));
    meanx = mean(mean(x(2:end-1,2:end-1)));
    plot(time./3600,meanf,'ro',time./3600,meanx,'bo','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    xlabel('Time [hr]',TX{:},FS{:});
    ylabel('$\phi$, $\chi$ [vol]',TX{:},FS{:});
end

fh8 = figure(8);
if iter > 0 % don't plot before solver has run first time
    if step==0
        Mass0 = sum(sum(rho(2:end-1,2:end-1)));
        sumH0 = sum(sum(H(2:end-1,2:end-1)));
        sumC0 = sum(sum(C(2:end-1,2:end-1)));
        sumV0 = sum(sum(V(2:end-1,2:end-1)));
    end
    Mass = sum(sum(rho(2:end-1,2:end-1)));
    subplot(3,1,1);
    plot(time./3600,sum(sum(H(2:end-1,2:end-1)))./Mass-(sumH0./Mass0),'bo','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $H$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(3,1,2);
    plot(time./3600,sum(sum(C(2:end-1,2:end-1)))./Mass-(sumC0./Mass0),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $C$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
    subplot(3,1,3);
    plot(time./3600,sum(sum(V(2:end-1,2:end-1)))./Mass-(sumV0./Mass0),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
    ylabel('consv. $V$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
end

drawnow

% save output to file
if save_op
    name = ['../out/',runID,'/',runID,'_vp_',num2str(floor(step/nop))];
    print(fh1,name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
    print(fh2,name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_mat_',num2str(floor(step/nop))];
    print(fh3,name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_chm',num2str(floor(step/nop))];
    print(fh4,name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
    print(fh6,name,'-dpng','-r300','-opengl');
    
    name = ['../out/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','f','x','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','dHdt','dCdt','dVdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step');
    name = ['../out/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','f','x','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','dHdt','dCdt','dVdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step');
    
    if step == 1
        logfile = ['../out/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end
    