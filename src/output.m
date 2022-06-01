
if plot_op
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
    LW = {'LineWidth',2};
    
    if Nx <= 10 && Nz <= 10  % create 0D plots
            
        fh1 = figure(1); clf;
        subplot(4,1,1)
        plot(hist.time/hr,hist.T(:,2)-273.15,CL{[1,2]},LW{:}); axis xy tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        plot(hist.time/hr,hist.cx(:,2)*100,CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
        plot(hist.time/hr,hist.cm(:,2)*100,CL{[1,3]},LW{:});
        plot(hist.time/hr,hist.c(:,2)./(1-hist.f(:,2))*100,CL{[1,2]},LW{:});
        title('$\bar{c}/(1-f)$ [wt\% SiO$_2$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        semilogy(hist.time/hr,hist.vf(:,2)*100,CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
        semilogy(hist.time/hr,hist.vm(:,2)*100,CL{[1,3]},LW{:});
        semilogy(hist.time/hr,hist.v(:,2)./(1-hist.x(:,2))*100,CL{[1,2]},LW{:});
        title('$\bar{v}/(1-x)$ [wt\% H$_2$O]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(hist.time/hr,hist.mu (:,2)*100.*(hist.mu (:,2)>1e-9),CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
        plot(hist.time/hr,hist.phi(:,2)*100.*(hist.phi(:,2)>1e-9),CL{[1,5]},LW{:});
        plot(hist.time/hr,hist.chi(:,2)*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:});
        title(['$\mu$, $\chi$, $\phi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
        fh2 = figure(2); clf;
        subplot(4,1,1)
        plot(hist.time/hr,hist.rhom(:,2),'-',CL{[1,3]},LW{:}); axis xy tight; box on; hold on
        plot(hist.time/hr,hist.rhox(:,2),'-',CL{[1,4]},LW{:});
        plot(hist.time/hr,hist.rho (:,2),'-',CL{[1,2]},LW{:});
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,2)
        semilogy(hist.time/hr,hist.eta(:,2),'k-',LW{:}); axis xy tight; box on;
        title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,3)
        plot(hist.time/hr,    hist.Gx(:,2)./hist.rho(:,2)*hr.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
        plot(hist.time/hr,10.*hist.Gf(:,2)./hist.rho(:,2)*hr.*(hist.phi(:,2)>1e-9),CL{[1,5]},LW{:});
        title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt \%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(4,1,4)
        plot(hist.time/hr,hist.dV(:,2)*hr,'k-',LW{:}); axis xy tight; box on;
        title('$\dot{V}$ [\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
    elseif Nx <= 10  % create 1D plots
        
        fh1 = figure(1); clf;
        subplot(1,5,1)
        plot(mean(T(2:end-1,2:end-1),2)-273.15,Z(2:end-1).',CL{[1,2]},LW{:}); axis ij tight; box on;
        title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,2)
        plot(mean(cx(2:end-1,2:end-1),2)*100,Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        plot(mean(cm(2:end-1,2:end-1),2)*100,Z(2:end-1).',CL{[1,3]},LW{:});
        plot(mean(c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)),2)*100,Z(2:end-1).',CL{[1,2]},LW{:});
        title('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        semilogx(mean(max(1e-6,vf(2:end-1,2:end-1)*100),2),Z(2:end-1).',CL{[1,5]},LW{:}); axis ij tight; box on; hold on;
        semilogx(mean(max(1e-6,vm(2:end-1,2:end-1)*100),2),Z(2:end-1).',CL{[1,3]},LW{:});
        semilogx(mean(max(1e-6,v(2:end-1,2:end-1)./(1-x(2:end-1,2:end-1))*100),2),Z(2:end-1).',CL{[1,2]},LW{:});
        title('$\bar{v}/(1-x)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        plot(mean(mu (2:end-1,2:end-1),2)*100.*(mean(mu (2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
        plot(mean(phi(2:end-1,2:end-1),2)*100.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',CL{[1,5]},LW{:});
        plot(mean(chi(2:end-1,2:end-1),2)*100.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',CL{[1,4]},LW{:});
        title(['$\mu$, $\phi$, $\chi$ [vol\%]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,5)
        plot(mean(-W(:,2:end-1),2)*hr,Zfc.',CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
        plot(mean(-(f(1:end-1,2:end-1)+f(2:end,2:end-1))/2.*wf(:,2:end-1),2)*hr,Zfc.',CL{[1,5]},LW{:});
        plot(mean(-(x(1:end-1,2:end-1)+x(2:end,2:end-1))/2.*wx(:,2:end-1),2)*hr,Zfc.',CL{[1,4]},LW{:});
        title('$W$, $w_\Delta^f$, $w_\Delta^x$ [m/hr]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh2 = figure(3); clf;
        subplot(1,5,1)
        plot(mean(rhox(2:end-1,2:end-1),2),Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        plot(mean(rhom(2:end-1,2:end-1),2),Z(2:end-1).',CL{[1,3]},LW{:});
        plot(mean(rho (2:end-1,2:end-1),2),Z(2:end-1).',CL{[1,2]},LW{:});
        title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,2)
        semilogx(geomean(min(eta(2:end-1,2:end-1),etam(2:end-1,2:end-1)),2),Z(2:end-1).',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
        semilogx(geomean(eta(2:end-1,2:end-1),2),Z(2:end-1).',CL{[1,2]},LW{:});
        title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        plot(mean(    Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*hr.*(mean(chi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        plot(mean(10.*Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1),2)*hr.*(mean(phi(2:end-1,2:end-1),2)>1e-9),Z(2:end-1).',CL{[1,5]},LW{:}); 
        title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        plot(mean(VolSrc(2:end-1,2:end-1),2),Z(2:end-1).',CL{[1,2]},LW{:}); axis ij tight; box on;
        title('$\dot{V}$ [1/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,5)
        plot(mean(P(2:end-1,2:end-1),2)/1e3,Z(2:end-1).',CL{[1,2]},LW{:}); axis ij tight; box on;
        title('$P$ [kPa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        
        fh3 = figure(4); clf;
        subplot(1,5,1)
        semilogx(mean(max(1e-3,min(1e3,itx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(mean(max(1e-3,min(1e3,itm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,3]},LW{:});
        semilogx(mean(max(1e-3,min(1e3,it (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',CL{[1,2]},LW{:});
        title('incomp. trace',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,2)
        semilogx(mean(max(1e-3,min(1e3,ctx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(mean(max(1e-3,min(1e3,ctm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,3]},LW{:});
        semilogx(mean(max(1e-3,min(1e3,ct (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',CL{[1,2]},LW{:});
        title('comp. trace',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,3)
        plot(mean(si(2:end-1,2:end-1),2),Z(2:end-1).',CL{[1,2]},LW{:}); axis ij tight; box on;
        title('stable isotope',TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,4)
        semilogx(mean(max(1e-3,min(1e3,ripx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(mean(max(1e-3,min(1e3,ripm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,3]},LW{:});
        semilogx(mean(max(1e-3,min(1e3,rip (2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',CL{[1,2]},LW{:});
        title(['radiogenic parent'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        subplot(1,5,5)
        semilogx(mean(max(1e-3,min(1e3,ridx(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
        semilogx(mean(max(1e-3,min(1e3,ridm(2:end-1,2:end-1))),2)*100,Z(2:end-1).',CL{[1,3]},LW{:});
        semilogx(mean(max(1e-3,min(1e3,rid(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)))),2)*100,Z(2:end-1).',CL{[1,2]},LW{:});
        title(['radiogenic daughter'],TX{:},FS{:}); set(gca,TL{:},TS{:});
        
    else % create 2D plots
        
    % set axis and border dimensions
    axh = 6.00; axw = axh*L/D;
    ahs = 0.40; avs = 0.2;
    axb = 1.00; axt = 0.4;
    axl = 1.20; axr = 0.4;
    
    % initialize figures and axes
    fh1 = figure(1); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh1,UN{:},'Position',[1 1 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off');
    set(fh1,'Resize','off');
    ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(13) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(14) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    fh2 = figure(2); clf; colormap(ocean);
    fh = axb + 1*axh + 0*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh2,UN{:},'Position',[3 3 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off');
    set(fh2,'Resize','off');
    ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(23) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    fh3 = figure(3); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh3,UN{:},'Position',[5 5 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off');
    set(fh3,'Resize','off');
    ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(33) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(34) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    
    fh4 = figure(4); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh4,UN{:},'Position',[7 7 fw fh]);
    set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh4,'Color','w','InvertHardcopy','off');
    set(fh4,'Resize','off');
    ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(43) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(44) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    
    fh5 = figure(5); clf; colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh5,UN{:},'Position',[9 9 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off');
    set(fh5,'Resize','off');
    ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(53) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(54) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(55) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(56) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    
    if plot_cv
        fh6 = figure(6); clf; colormap(ocean);
        fh = axb + 1*axh + 0*avs + axt;
        fw = axl + 3*axw + 2*ahs + axr;
        set(fh6,UN{:},'Position',[11 11 fw fh]);
        set(fh6,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
        set(fh6,'Color','w','InvertHardcopy','off');
        set(fh6,'Resize','off');
        ax(61) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
        ax(62) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
        ax(63) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);
    end
    
    % plot velocity-pressure solution in Fig. 1
    figure(1);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(11));
    imagesc(X(2:end-1),Z(2:end-1),-W(:      ,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(12));
    imagesc(X(2:end-1),Z(2:end-1), U(2:end-1,:      ).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [m/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(13));
    imagesc(X(2:end-1),Z(2:end-1), P(2:end-1,2:end-1)./1e3); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [kPa]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(14));
    imagesc(X(2:end-1),Z(2:end-1),Div_V(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \mathbf{v}$ [1/s]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    % plot temperature and composition in Fig. 2
    figure(2);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(21));
    imagesc(X(2:end-1),Z(2:end-1),T(2:end-1,2:end-1)-273.15); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T [^\circ$C]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(22));
    imagesc(X(2:end-1),Z(2:end-1),c(2:end-1,2:end-1)./(1-f(2:end-1,2:end-1)).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{c}/(1-f)$ [wt\% SiO$_2$]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(23));
    imagesc(X(2:end-1),Z(2:end-1),v(2:end-1,2:end-1).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{v}$ [wt\% H$_2$O]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot phase fractions and reaction rates in Fig. 3
    figure(3);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(31));
    imagesc(X(2:end-1),Z(2:end-1),chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(32));
    imagesc(X(2:end-1),Z(2:end-1),phi(2:end-1,2:end-1).*100.*(phi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\phi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(33));
    imagesc(X(2:end-1),Z(2:end-1),Gx(2:end-1,2:end-1)./rho(2:end-1,2:end-1)*hr*100.*(chi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(34));
    imagesc(X(2:end-1),Z(2:end-1),Gf(2:end-1,2:end-1)./rho(2:end-1,2:end-1)*hr*100.*(phi(2:end-1,2:end-1)>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_f/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot density, rheology, and segregation speeds in Fig. 4
    figure(4);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(41));
    imagesc(X(2:end-1),Z(2:end-1),      rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(42));
    imagesc(X(2:end-1),Z(2:end-1),log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(43));
    imagesc(X(2:end-1),Z(2:end-1),-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(44));
    imagesc(X(2:end-1),Z(2:end-1),-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^f$ [m/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot geochemical variables in Fig. 5
    figure(5);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    axes(ax(51));
    imagesc(X(2:end-1),Z(2:end-1),it(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['incomp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(52));
    imagesc(X(2:end-1),Z(2:end-1),ct(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['comp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(53));
    imagesc(X(2:end-1),Z(2:end-1),si(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['stable isotope'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(54));
    imagesc(X(2:end-1),Z(2:end-1),rip(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. parent'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(55));
    imagesc(X(2:end-1),Z(2:end-1),rid(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. daughter'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(56));
    imagesc(X(2:end-1),Z(2:end-1),(dcy_rip(2:end-1,2:end-1)-dcy_rid(2:end-1,2:end-1))./(dcy_rip(2:end-1,2:end-1)+dcy_rid(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['radiogen. disequilibrium'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    
    if plot_cv && iter > 0
        % plot residual fields in Fig. 4
        figure(6);
        axes(ax(61));
        imagesc(X(2:end-1),Z(2:end-1),-res_W(:      ,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $W$'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
        axes(ax(62));
        imagesc(X(2:end-1),Z(2:end-1), res_U(2:end-1,:      )); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $U$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
        axes(ax(63));
        imagesc(X(2:end-1),Z(2:end-1), res_P(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $P$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); 
    end
    
    end
    
    % plot phase diagram
    fh7 = figure(7); clf;
    TT = linspace(Tphs0+Ptop*clap,Tphs1+Ptop*clap,500);
    cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perT-Tphs0)./(Tphs1-Tphs0)*500)),linspace((perCx+perCm)/2,cphs0,floor((perT-Tphs1)./(Tphs0-Tphs1)*500))];
    [~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,TINY);
    plot(CCx,TT,'k-',LW{:}); axis tight; hold on; box on;
    plot(CCm,TT,'k-',LW{:});
    perTs  = perT;
    Tphs0s = Tphs0;
    Tphs1s = Tphs1;
    vv = 0.10*ones(size(TT));
    xx = 0.50*ones(size(TT));
    ff = 0.05*ones(size(TT));
    for i = 1:5
        TT = linspace(Tphs0s+Ptop*clap,Tphs1s+Ptop*clap,500);
        cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perTs-Tphs0s)./(Tphs1s-Tphs0s)*500)),linspace((perCx+perCm)/2,cphs0,floor((perTs-Tphs1s)./(Tphs0s-Tphs1s)*500))];
        vmq_c0 = (4.7773e-7.*Ptop.^0.6 + 1e-11.*Ptop) .* exp(2565*(1./(TT+273.15)-1./(perT+273.15))); % Katz et al., 2003; Moore et al., 1998
        vmq_c1 = (3.5494e-3.*Ptop.^0.5 + 9.623e-8.*Ptop - 1.5223e-11.*Ptop.^1.5)./(TT+273.15) + 1.2436e-14.*Ptop.^1.5; % Liu et al., 2015
        vmq0   = (1-cc).*vmq_c0 + cc.*vmq_c1;
        perTs  = perT -dTH2O(2).*vmq0(round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*500)).^0.75;
        Tphs0s = Tphs0-dTH2O(1).*vmq0(1).^0.75;
        Tphs1s = Tphs1-dTH2O(3).*vmq0(end).^0.75;
    end
    [~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,vmq0,Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,TINY);
    plot(CCx,TT,'k-',LW{:}); axis tight; hold on; box on;
    plot(CCm,TT,'k-',LW{:});

    Tplt = T-273.15 - (Pt-Ptop)*clap;
    cplt = c./(1-f);
    plot(cplt(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'.',CL{[1,2]},LW{:},'MarkerSize',15);
    plot(cx  (2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'.',CL{[1,4]},LW{:},'MarkerSize',15);
    plot(cm  (2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'.',CL{[1,3]},LW{:},'MarkerSize',15);
    
    set(gca,'TickLabelInterpreter','latex','FontSize',15)
    title('Phase Diagram','Interpreter','latex','FontSize',22)
    xlabel('Composition','Interpreter','latex','FontSize',18)
    ylabel('Temperature','Interpreter','latex','FontSize',18)
    
    % plot model history
    if plot_cv
        figure(8); clf;
        subplot(4,1,1);
        plot(hist.time/hr,hist.DM./hist.sumM,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(4,1,2);
        plot(hist.time/hr,hist.DH./hist.sumH,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $H$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(4,1,3);
        plot(hist.time/hr,hist.DC./hist.sumC,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(4,1,4);
        plot(hist.time/hr,hist.DV./hist.sumV,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('consv. $V$',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
        
        figure(9); clf;
        subplot(4,1,1);
        plot(hist.time/hr,hist.EM,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('error $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(4,1,2);
        plot(hist.time/hr,hist.EH,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('error $H$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(4,1,3);
        plot(hist.time/hr,hist.EC,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('error $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
        subplot(4,1,4);
        plot(hist.time/hr,hist.EV,'k-',LW{:}); hold on; axis tight; box on;
        ylabel('error $V$',TX{:},FS{:}); set(gca,TL{:},TS{:});
        xlabel('Time [hr]',TX{:},FS{:});
    end
    
    drawnow

end

% save output to file
if save_op
    if plot_op
        if Nx <= 10 && Nz <= 10  % print 0D plots
            name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
            print(fh1,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
            print(fh2,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
            print(fh7,name,'-dpng','-r300','-opengl');
        elseif Nx <= 10  % create 1D plots
            name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
            print(fh1,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
            print(fh2,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_gch_',num2str(floor(step/nop))];
            print(fh3,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
            print(fh7,name,'-dpng','-r300','-opengl');
        else
            name = [opdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
            print(fh1,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
            print(fh2,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
            print(fh3,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
            print(fh4,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
            print(fh5,name,'-dpng','-r300','-opengl');
            name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
            print(fh7,name,'-dpng','-r300','-opengl');
        end
    end
    
    name = [opdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
    name = [opdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SI','RIP','RID','it','ct','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSIdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist','VolSrc','wf','wx');
    
    if step == 0
        logfile = [opdir,'/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

clear fh1 fh2 fh3 fh4 fh5 fh6 fh7
    