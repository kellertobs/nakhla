
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

if Nx <= 1 && Nz <= 1  % create 0D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    subplot(3,1,1)
    plot(hist.time/hr,hist.T(:,2)-273.15,CL{[1,2]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/hr,hist.Tliq(:,2),CL{[1,3]},LW{:});
    plot(hist.time/hr,hist.Tsol(:,2),CL{[1,4]},LW{:});
    title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(3,1,2)
    for i=1:cal.ncmp
        plot(hist.time/hr,squeeze(hist.cm(:,2,i)./sum(hist.cm(:,2,1:end-1),3)).*100,'--',LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis xy tight; box on; hold on
        plot(hist.time/hr,squeeze(hist.cx(:,2,i)./sum(hist.cx(:,2,1:end-1),3)).*100,'-.',LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis xy tight; box on; hold on
        pcmp(i) = plot(hist.time/hr,squeeze(hist. c(:,2,i)./sum(hist. c(:,2,1:end-1),3))./(1-hist.f(:,2)).*100,'-' ,LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis xy tight; box on; hold on
    end
    title('cmps [wt\%]',TX{:},FS{:}); legend(pcmp,cal.cmpStr(1:end),TX{:},FS{1},8,'Location','northeast'); set(gca,TL{:},TS{:});
    subplot(3,1,3)
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
    subplot(3,1,1)
    for i=1:cal.noxd-1
        plot(hist.time/hr,squeeze(hist.cm_oxd(:,2,i)./sum(hist.cm_oxd(:,2,1:end-1),3)).*100,'-',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis xy tight; box on; hold on
    end
    title('Melt oxds [wt\%]',TX{:},FS{:}); legend(cal.oxdStr(1:end-1),TX{:},FS{1},8,'Location','northeast'); set(gca,TL{:},TS{:});
    subplot(3,1,2)
    for i=1:cal.noxd-1
        plot(hist.time/hr,squeeze(hist.cx_oxd(:,2,i)./sum(hist.cx_oxd(:,2,1:end-1),3)).*100,'-',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis xy tight; box on; hold on
    end
    title('Xtal oxds [wt\%]',TX{:},FS{:});
    subplot(3,1,3)
    for i=1:cal.nmsy
        plot(hist.time/hr,squeeze(hist.x(:,2).*hist.cx_msy(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.nmsy)+1,:)); axis xy tight; box on; hold on
    end
    title('Xtal msys [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); legend(cal.msyStr,TX{:},FS{1},8,'Location','northeast'); set(gca,TL{:},TS{:});
    xlabel('Time [hr]',TX{:},FS{:});

elseif Nx <= 1  % create 1D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    plot(T-273.15,Zc.',CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
    plot(cal.Tliq,Zc.',CL{[1,3]},LW{:});
    plot(cal.Tsol,Zc.',CL{[1,4]},LW{:});
    title('$T [^\circ$C]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    for i=1:cal.ncmp
        plot(squeeze(c(:,:,i)).*100.*(1+9.*cfq(:,:,i)),Zc.',LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis ij tight; box on; hold on;
    end
    title('Bulk cmps [wt\%]',TX{:},FS{:}); legend(cal.cmpStr,TX{:},FS{:},'Location','west'); ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    plot(mu *100.*(mu >1e-9),Zc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(chi*100.*(chi>1e-9),Zc.',CL{[1,4]},LW{:});
    plot(phi*100.*(phi>1e-9),Zc.',CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,4)
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
    else; set(0, 'CurrentFigure', fh3); clf;
    end 
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    for i=1:cal.noxd-1
        plot(squeeze( c_oxd(:,:,i)./sum(c_oxd(:,:,1:end-1),3)).*100,Zc.',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis ij tight; box on; hold on;
    end
    title('Bulk oxds [wt\%]',TX{:},FS{:});ylabel('Depth [m]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    for i=1:cal.noxd-1
        plot(squeeze(cm_oxd(:,:,i)./sum(cm_oxd(:,:,1:end-1),3)).*100,Zc.',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis ij tight; box on; hold on;
    end
    title('Melt oxds [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    for i=1:cal.noxd-1
        plot(squeeze(cx_oxd(:,:,i)./sum(cx_oxd(:,:,1:end-1),3)).*100,Zc.',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis ij tight; box on; hold on;
    end
    title('Xtal oxds [wt\%]',TX{:},FS{:}); legend(cal.oxdStr(1:end-1),TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});
    subplot(1,4,4)
    for i=1:cal.nmsy
        plot(squeeze(x.*cx_msy(:,:,i)),Zc.',LW{:},'color',ocean(round((i-1)*213/cal.nmsy)+1,:)); axis ij tight; box on; hold on;
    end
    title('Xtal msys [wt\%]',TX{:},FS{:}); legend(cal.msyStr,TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});

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
    fw = axl + 2*axw + 1*ahs + axr;
    set(fh5,UN{:},'Position',[9 9 fw fh]);
    set(fh5,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh5,'Color','w','InvertHardcopy','off','Resize','off');
    ax(51) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(52) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(53) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(54) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh6','var'); fh6 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh6); clf;
    end 
    colormap(ocean);
    fh = axb + 3*axh + 2*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh6,UN{:},'Position',[11 11 fw fh]);
    set(fh6,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh6,'Color','w','InvertHardcopy','off','Resize','off');
    ax(61) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+2*axh+2*avs axw axh]);
    ax(62) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+2*axh+2*avs axw axh]);
    ax(63) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+2*axh+2*avs axw axh]);
    ax(64) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(65) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(66) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(67) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(68) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(69) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh7','var'); fh7 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh7); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh7,UN{:},'Position',[13 13 fw fh]);
    set(fh7,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh7,'Color','w','InvertHardcopy','off','Resize','off');
    ax(71) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(72) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(73) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(74) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(75) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(76) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    if ~exist('fh8','var'); fh8 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh8); clf;
    end 
    colormap(ocean);
    fh = axb + 2*axh + 1*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh8,UN{:},'Position',[15 15 fw fh]);
    set(fh8,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh8,'Color','w','InvertHardcopy','off','Resize','off');
    ax(81) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(82) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(83) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(84) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(85) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(86) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

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
    imagesc(Xc,Zc,squeeze(c_oxd(:,:,cal.Si)./sum(c_oxd(:,:,1:end-1),3).*100)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{c}$ [wt\% SiO$_2$]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    set(fh2,'CurrentAxes',ax(23));
    imagesc(Xc,Zc,squeeze(c_oxd(:,:,cal.H))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{v}$ [wt\% H$_2$O]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
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

    % plot pseudo-component composition in Fig. 5
    set(0,'CurrentFigure',fh5)
    set(fh5,'CurrentAxes',ax(51));
    imagesc(Xc./1e3,Zc./1e3,c(:,:,1).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(52));
    imagesc(Xc./1e3,Zc./1e3,c(:,:,2).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(53));
    imagesc(Xc./1e3,Zc./1e3,c(:,:,3).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{3},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [km]',TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(54));
    imagesc(Xc./1e3,Zc./1e3,c(:,:,4).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{4},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot major oxide composition in Fig. 6
    set(0,'CurrentFigure',fh6)
    set(fh6,'CurrentAxes',ax(61));
    sumanh = sum(c_oxd(:,:,1:end-1),3);
    imagesc(Xc,Zc,c_oxd(:,:,1)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]);
    set(fh6,'CurrentAxes',ax(62));
    imagesc(Xc,Zc,c_oxd(:,:,2)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(63));
    imagesc(Xc,Zc,c_oxd(:,:,3)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{3},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(64));
    imagesc(Xc,Zc,c_oxd(:,:,4)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{4},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [km]',TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(65));
    imagesc(Xc,Zc,c_oxd(:,:,5)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{5},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]); 
    set(fh6,'CurrentAxes',ax(66));
    imagesc(Xc,Zc,c_oxd(:,:,6)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{6},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(67));
    imagesc(Xc,Zc,c_oxd(:,:,7)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{7},' [wt\%]'],TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(68));
    imagesc(Xc,Zc,c_oxd(:,:,8)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{8},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [km]',TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(69));
    imagesc(Xc,Zc,c_oxd(:,:,9)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{9},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot mineral assemblage in Fig. 7
    set(0,'CurrentFigure',fh7)
    set(fh7,'CurrentAxes',ax(71));
    imagesc(Xc,Zc,cx_msy(:,:,1).*x); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{1},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh7,'CurrentAxes',ax(72));
    imagesc(Xc,Zc,cx_msy(:,:,2).*x); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{2},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh7,'CurrentAxes',ax(73));
    imagesc(Xc,Zc,cx_msy(:,:,3).*x); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{3},' [wt\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh7,'CurrentAxes',ax(74));
    imagesc(Xc,Zc,cx_msy(:,:,4).*x); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{4},' [wt\%]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    set(fh7,'CurrentAxes',ax(75));
    imagesc(Xc,Zc,cx_msy(:,:,5).*x); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{5},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    set(fh7,'CurrentAxes',ax(76));
    imagesc(Xc,Zc,cx_msy(:,:,6).*x); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{6},' [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');

    % plot geochemical variables in Fig. 7
    set(0,'CurrentFigure',fh8)
    set(fh8,'CurrentAxes',ax(81));
    imagesc(Xc,Zc,te(:,:,1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 1'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:}); 
    set(fh8,'CurrentAxes',ax(82));
    imagesc(Xc,Zc,te(:,:,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 2'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh8,'CurrentAxes',ax(83));
    imagesc(Xc,Zc,ir(:,:,1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['isotope ratio 1'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh8,'CurrentAxes',ax(84));
    imagesc(Xc,Zc,te(:,:,3)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 3'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    set(fh8,'CurrentAxes',ax(85));
    imagesc(Xc,Zc,te(:,:,4)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['trace element 4'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    set(fh8,'CurrentAxes',ax(86));
    imagesc(Xc,Zc,ir(:,:,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['isotope ratio 2'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    sgtitle(['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k');


end

% plot phase diagram
if ~exist('fh9','var'); fh9 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh9);
end
if Nz>1 || step==0; clf; end

plot( c_oxd(:,:,cal.Si)./sum( c_oxd(:,:,1:end-1),3).*100,T-273.15,'.',CL{[1,2]},LW{:},'MarkerSize',15); axis tight; hold on
plot(cx_oxd(:,:,cal.Si)./sum(cx_oxd(:,:,1:end-1),3).*100,T-273.15,'.',CL{[1,4]},LW{:},'MarkerSize',15);
plot(cm_oxd(:,:,cal.Si)./sum(cm_oxd(:,:,1:end-1),3).*100,T-273.15,'.',CL{[1,3]},LW{:},'MarkerSize',15);

set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

% plot composition on TAS, AFM diagrams
if ~exist('fh10','var'); fh10 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh10);
end
if Nz>1 || step==0; clf;
TAS; axis tight; box on; hold on;
end
cxSi = cx_oxd(:,:,cal.Si)./sum(cx_oxd(:,:,1:end-1),3).*100;
cmSi = cm_oxd(:,:,cal.Si)./sum(cm_oxd(:,:,1:end-1),3).*100;
 cSi =  c_oxd(:,:,cal.Si)./sum( c_oxd(:,:,1:end-1),3).*100;
cxNK = sum(cx_oxd(:,:,[cal.Na,cal.K]),3)./sum(cx_oxd(:,:,1:end-1),3).*100;
cmNK = sum(cm_oxd(:,:,[cal.Na,cal.K]),3)./sum(cm_oxd(:,:,1:end-1),3).*100;
 cNK = sum( c_oxd(:,:,[cal.Na,cal.K]),3)./sum( c_oxd(:,:,1:end-1),3).*100;
scatter(cxSi(:),cxNK(:),50,T(:)-273.15,'filled','^'); colormap(ocean); cb = colorbar;
scatter(cmSi(:),cmNK(:),50,T(:)-273.15,'filled','o');
scatter( cSi(:), cNK(:),80,T(:)-273.15,'filled','s');
set(cb,TL{:},'FontSize',12); set(gca,TL{:},'FontSize',15); xlabel('SiO$_2$ [wt \%]',TX{:},'FontSize',15); ylabel('Na$_2$O + K$_2$O [wt \%]',TX{:},'FontSize',15);

if ~exist('fh11','var'); fh11 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh11);
end
if Nz>1 || step==0; clf;
AFM; axis tight; box on; hold on;
end
[A,B] = terncoords(cx_oxd(:,:, cal.Mg          )./sum(cx_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3), ...
                   cx_oxd(:,:, cal.Fe          )./sum(cx_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3), ...
               sum(cx_oxd(:,:,[cal.Na,cal.K]),3)./sum(cx_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3));
scatter(A(:),B(:),50,T(:)-273.15,'filled','^'); colormap(ocean); cb = colorbar;
[A,B] = terncoords(cm_oxd(:,:, cal.Mg          )./sum(cm_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3), ...
                   cm_oxd(:,:, cal.Fe          )./sum(cm_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3), ...
               sum(cm_oxd(:,:,[cal.Na,cal.K]),3)./sum(cm_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3));
scatter(A(:),B(:),50,T(:)-273.15,'filled','o'); colormap(ocean);
[A,B] = terncoords(c_oxd(:,:, cal.Mg          )./(sum(c_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3)), ...
                   c_oxd(:,:, cal.Fe          )./(sum(c_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3)), ...
               sum(c_oxd(:,:,[cal.Na,cal.K]),3)./(sum(c_oxd(:,:,[cal.Fe,cal.Mg,cal.Na,cal.K]),3)));
scatter(A(:),B(:),80,T(:)-273.15,'filled','s'); colormap(ocean);
set(cb,TL{:},'FontSize',12); set(gca,TL{:},'FontSize',15); xlabel('SiO$_2$ [wt \%]',TX{:},'FontSize',15); ylabel('Na$_2$O + K$_2$O [wt \%]',TX{:},'FontSize',15);

% plot model history
if plot_cv
    if ~exist('fh12','var'); fh12 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh12); clf;
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

    if ~exist('fh13','var'); fh13 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh13); clf;
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
        name = [opdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
        print(fh10,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
        print(fh11,name,'-dpng','-r300','-image');
    elseif Nx <= 10  % create 1D plots
        name = [opdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
        print(fh10,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
        print(fh11,name,'-dpng','-r300','-image');
    else
        name = [opdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
        print(fh4,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_cmp',num2str(floor(step/nop))];
        print(fh5,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_oxd',num2str(floor(step/nop))];
        print(fh6,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_mnr',num2str(floor(step/nop))];
        print(fh7,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
        print(fh8,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
        print(fh10,name,'-dpng','-r300','-image');
        name = [opdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
        print(fh11,name,'-dpng','-r300','-image');
    end

    name = [opdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','c','cm','cx','cf','TE','IR','te','ir','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTEdt','dIRdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
    name = [opdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','c','cm','cx','cf','TE','IR','te','ir','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTEdt','dIRdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
    name = [opdir,'/',runID,'/',runID,'_hist'];
    save(name,'hist');

    if step == 0
        logfile = [opdir,'/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end
    