
% prepare for plotting
TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
UN = {'Units','Centimeters'};
CL = {'Color',[0.0 0.0 0.0],[0.80 0.15 0.10],[0.10 0.15 0.65],[0.45 0.60 0.95]};
LW = {'LineWidth',2};
if plot_op
    VIS = {'Visible','on'};
else
    VIS = {'Visible','off'};
end

if Nx <= 10 && Nz <= 10  % create 0D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
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
    plot(hist.time/hr,hist.chi(:,2)*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:});
    plot(hist.time/hr,hist.phi(:,2)*100.*(hist.phi(:,2)>1e-9),CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end
    subplot(4,1,1)
    plot(hist.time/hr,hist.rhom(:,2),'-',CL{[1,3]},LW{:}); axis xy tight; box on; hold on
    plot(hist.time/hr,hist.rhox(:,2),'-',CL{[1,4]},LW{:});
    plot(hist.time/hr,hist.rho (:,2),'-',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    semilogy(hist.time/hr,hist.eta (:,2),'k-',LW{:}); axis xy tight; box on; hold on
    semilogy(hist.time/hr,hist.etam(:,2),'r-',LW{:});
    title('$\eta^m,\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,3)
    plot(hist.time/hr,    hist.Gx(:,2)./hist.rho(:,2)*hr.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/hr,10.*hist.Gf(:,2)./hist.rho(:,2)*hr.*(hist.phi(:,2)>1e-9),CL{[1,5]},LW{:});
    title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt \%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,4)
    plot(hist.time/hr,hist.dV(:,2)*hr,'k-',LW{:}); axis xy tight; box on;
    title('$\dot{V}$ [\%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end
    subplot(4,1,1)
    for i=1:cal.nc
        plot(hist.time/hr,squeeze(hist.cx_cmp(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis xy tight; box on; hold on
    end
    title('Xtal cmps [wt\%]',TX{:},FS{:}); legend(cal.cmpStr,TX{:},FS{1},8,'Location','northeast'); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    for i=1:cal.nc
        plot(hist.time/hr,squeeze(hist.cm_cmp(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis xy tight; box on; hold on
    end
    title('Melt cmps [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,3)
    for i=1:cal.nc
        plot(hist.time/hr,squeeze(hist.cx_oxd(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis xy tight; box on; hold on
    end
    title('Xtal oxds [wt\%]',TX{:},FS{:}); legend(cal.oxdStr,TX{:},FS{1},8,'Location','northeast'); set(gca,TL{:},TS{:});
    subplot(4,1,4)
    for i=1:cal.nc
        plot(hist.time/hr,squeeze(hist.cm_oxd(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis xy tight; box on; hold on
    end
    title('Melt oxds [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});

elseif Nx <= 10  % create 1D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    subplot(1,5,1)
    plot(T-273.15,Zc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    plot(cx*100,Zc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(cm*100,Zc.',CL{[1,3]},LW{:});
    plot(c./(1-f)*100,Zc.',CL{[1,2]},LW{:});
    title('$\bar{c}/(1-f)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,3)
    semilogx(max(1e-6,vf*100).*any(v(:)>10*TINY),Zc.',CL{[1,5]},LW{:}); axis ij tight; box on; hold on;
    semilogx(max(1e-6,vm*100).*any(v(:)>10*TINY),Zc.',CL{[1,3]},LW{:});
    semilogx(max(1e-6,v ./(1-x)*100).*any(v(:)>10*TINY),Zc.',CL{[1,2]},LW{:});
    title('$\bar{v}/(1-x)$ [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,4)
    plot(mu *100.*(mu >1e-9),Zc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(chi*100.*(chi>1e-9),Zc.',CL{[1,4]},LW{:});
    plot(phi*100.*(phi>1e-9),Zc.',CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,5)
    plot(-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1)*hr,Zf.',CL{[1,5]},LW{:}); axis ij tight; box on; hold on;
    plot(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)*hr,Zf.',CL{[1,4]},LW{:});
    plot(-(mu ([1,1:end],:)+mu ([1:end,end],:))/2.*wm(:,2:end-1)*hr,Zf.',CL{[1,3]},LW{:});
    plot(-                                         W (:,2:end-1)*hr,Zf.',CL{[1,2]},LW{:});
    title('$W$, $w_\Delta^f$, $w_\Delta^x$ [m/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end 
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    subplot(1,5,1)
    plot(rhox,Zc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(rhom,Zc.',CL{[1,3]},LW{:});
    plot(rho ,Zc.',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    semilogx(min(eta,etam),Zc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    semilogx(eta,Zc.',CL{[1,2]},LW{:});
    title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,3)
    plot(    Gx./rho*hr.*(chi>1e-9),Zc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(10.*Gf./rho*hr.*(phi>1e-9),Zc.',CL{[1,5]},LW{:});
    title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,4)
    plot(VolSrc,Zc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title('$\dot{V}$ [1/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,5)
    plot(P(2:end-1,2:end-1),Zc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title('$P$ [Pa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
 
    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh4); clf;
    end 
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    for i=1:cal.nc
        plot(squeeze( c_oxd(:,:,i)),Zc.',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis ij tight; box on; hold on;
    end
    title('Bulk oxds [wt\%]',TX{:},FS{:});ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    for i=1:cal.nc
        plot(squeeze(cm_oxd(:,:,i)),Zc.',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis ij tight; box on; hold on;
    end
    title('Melt oxds [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    for i=1:cal.nc
        plot(squeeze(cx_oxd(:,:,i)),Zc.',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis ij tight; box on; hold on;
    end
    title('Xtal oxds [wt\%]',TX{:},FS{:}); legend(cal.oxdStr,TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});
    subplot(1,4,4)
    for i=1:cal.nc
        plot(squeeze(x.*cx_cmp(:,:,i)),Zc.',LW{:},'color',ocean(round((i-1)*213/cal.nc)+1,:)); axis ij tight; box on; hold on;
    end
    title('Xtal cmps [wt\%]',TX{:},FS{:}); legend(cal.cmpStr,TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});

else % create 2D plots

    % set axis and border dimensions
    axh = 6.00*sqrt(D/L); axw = 6.00*sqrt(L/D)+1.50;
    ahs = 0.60; avs = 0.80;
    axb = 1.20; axt = 1.50;
    axl = 1.20; axr = 0.60;

    % initialize figures and axes
    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh1,UN{:},'Position',[1 1 fw fh]);
    set(fh1,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh1,'Color','w','InvertHardcopy','off','Resize','off');
    ax(11) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(12) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(13) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(14) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end 
    colormap(ocean);
    fh = axb + 1*axh + 0*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh2,UN{:},'Position',[3 3 fw fh]);
    set(fh2,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh2,'Color','w','InvertHardcopy','off','Resize','off');
    ax(21) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(22) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(23) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh3,UN{:},'Position',[5 5 fw fh]);
    set(fh3,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh3,'Color','w','InvertHardcopy','off','Resize','off');
    ax(31) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(32) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(33) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(34) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh4','var'); fh4 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh4); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh4,UN{:},'Position',[7 7 fw fh]);
    set(fh4,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh4,'Color','w','InvertHardcopy','off','Resize','off');
    ax(41) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(42) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(43) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(44) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh5','var'); fh5 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh5); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh5,UN{:},'Position',[9 9 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off','Resize','off');
    ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(53) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(54) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(55) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(56) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    % plot velocity-pressure solution in Fig. 1
    set(0,'CurrentFigure',fh1)
    set(fh1,'CurrentAxes',ax(11));
    imagesc(Xc,Zc,-W(:      ,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh1,'CurrentAxes',ax(12));
    imagesc(Xc,Zc, U(2:end-1,:      ).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [m/hr]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh1,'CurrentAxes',ax(13));
    imagesc(Xc,Zc, P); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [Pa]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    set(fh1,'CurrentAxes',ax(14));
    imagesc(Xc,Zc,Div_V); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \mathbf{v}$ [1/s]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot temperature and composition in Fig. 2
    set(0,'CurrentFigure',fh2)
    set(fh2,'CurrentAxes',ax(21));
    imagesc(Xc,Zc,T-273.15); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T [^\circ$C]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh2,'CurrentAxes',ax(22));
    imagesc(Xc,Zc,c./(1-f).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{c}/(1-f)$ [wt\% SiO$_2$]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    set(fh2,'CurrentAxes',ax(23));
    imagesc(Xc,Zc,v./(1-x).*100.*(v>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{v}/(1-x)$ [wt\% H$_2$O]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot phase fractions and reaction rates in Fig. 3
    set(0,'CurrentFigure',fh3)
    set(fh3,'CurrentAxes',ax(31));
    imagesc(Xc,Zc,chi.*100.*(chi>1e-9) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh3,'CurrentAxes',ax(32));
    imagesc(Xc,Zc,phi.*100.*(phi>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\phi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh3,'CurrentAxes',ax(33));
    imagesc(Xc,Zc,Gx./rho*hr*100.*(chi>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(34));
    imagesc(Xc,Zc,Gf./rho*hr*100.*(phi>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_f/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot density, rheology, and segregation speeds in Fig. 4
    set(0,'CurrentFigure',fh4)
    set(fh4,'CurrentAxes',ax(41));
    imagesc(Xc,Zc,      rho ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\rho}$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh4,'CurrentAxes',ax(42));
    imagesc(Xc,Zc,log10(eta)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh4,'CurrentAxes',ax(43));
    imagesc(Xc,Zc,-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(44));
    imagesc(Xc,Zc,-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1).*hr.*(any(phi>1e-9))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^f$ [m/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot geochemical variables in Fig. 5
    set(0,'CurrentFigure',fh5)
    set(fh5,'CurrentAxes',ax(51));
    imagesc(Xc,Zc,te(:,:,1)./(1-f)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 1'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh5,'CurrentAxes',ax(52));
    imagesc(Xc,Zc,te(:,:,2)./(1-f)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 2'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(53));
    imagesc(Xc,Zc,ir(:,:,1)./(1-f)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['isotope ratio 1'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(54));
    imagesc(Xc,Zc,te(:,:,3)./(1-f)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 3'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(55));
    imagesc(Xc,Zc,te(:,:,4)./(1-f)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 4'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(56));
    imagesc(Xc,Zc,ir(:,:,2)./(1-f)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['isotope ratio 2'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

end

% plot phase diagram
if ~exist('fh7','var'); fh7 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh7); clf;
end
TT = [linspace(cal.Tphs0+Ptop*cal.clap,cal.perT+Ptop*cal.clap,200),linspace(cal.perT+Ptop*cal.clap,cal.Tphs1+Ptop*cal.clap,200)];
cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2,200),linspace((cal.perCx+cal.perCm)/2,cal.cphs0,200)];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,Ptop*ones(size(TT)),cal,TINY);
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);
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
plot(CCx.*100,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm.*100,TT,'k-','LineWidth',2);

Tplt = T-273.15 - (Pt-Ptop)*cal.clap;
cplt = c./(1-f);
plot(cplt.*100,Tplt,'.',CL{[1,2]},LW{:},'MarkerSize',15);
plot(cx  .*100,Tplt,'.',CL{[1,4]},LW{:},'MarkerSize',15);
plot(cm  .*100,Tplt,'.',CL{[1,3]},LW{:},'MarkerSize',15);

set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot model history
if plot_cv
    if ~exist('fh8','var'); fh8 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh8); clf;
    end 
    subplot(4,1,1);
    plot(hist.time/hr,hist.DM./hist.sumM,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    plot(hist.time/hr,hist.DS./hist.sumS,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    plot(hist.time/hr,hist.DC./hist.sumC,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    plot(hist.time/hr,hist.DV./hist.sumV,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $V$',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});

    if ~exist('fh9','var'); fh9 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh9); clf;
    end 
    subplot(4,1,1);
    plot(hist.time/hr,hist.EM,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,2);
    plot(hist.time/hr,hist.ES,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,3);
    plot(hist.time/hr,hist.EC,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(4,1,4);
    plot(hist.time/hr,hist.EV,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $V$',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});
end

drawnow

% save output to file
if save_op && ~restart
    if Nx <= 10 && Nz <= 10  % print 0D plots
        name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_cmp',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
        print(fh7,name,'-dpng','-r300','-image');
    elseif Nx <= 10  % create 1D plots
        name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
        print(fh7,name,'-dpng','-r300','-image');
    else
        name = [opdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
        print(fh4,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
        print(fh5,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
        print(fh7,name,'-dpng','-r300','-image');
    end

    name = [opdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','M','S','C','V','T','c','v','cm','cx','vm','vf','TE','IR','te','ir','dSdt','dCdt','dVdt','dFdt','dXdt','dMdt','drhodt','dTEdt','dIRdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm');
    name = [opdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','f','x','m','phi','chi','mu','X','F','M','S','C','V','T','c','v','cm','cx','vm','vf','TE','IR','te','ir','dSdt','dCdt','dVdt','dFdt','dXdt','dMdt','drhodt','dTEdt','dIRdt','Gf','Gx','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm');
    name = [opdir,'/',runID,'/',runID,'_hist'];
    save(name,'hist');

    if step == 0
        logfile = [opdir,'/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end
    