
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

% adjust scales and units for intuitive visualisation
if time < 1e3*hr
    TimeScale = hr;
    TimeUnits = 'hr';
elseif time >= 1e3*hr && time < 1e2*yr
    TimeScale = yr;
    TimeUnits = 'yr';
elseif time >= 1e2*yr
    TimeScale = 1e3*yr;
    TimeUnits = 'kyr';
end
if D < 1e3
    SpaceScale = 1;
    SpaceUnits = 'm';
elseif D >= 1e3
    SpaceScale = 1e3;
    SpaceUnits = 'km';
end
if max(Vel(:)) < 1000/yr
    SpeedScale = 1/yr;
    SpeedUnits = 'm/yr';
elseif max(Vel(:)) >= 1000/yr && max(Vel(:)) < 1000/hr
    SpeedScale = 1/hr;
    SpeedUnits = 'm/hr';
elseif max(Vel(:)) >= 1000/hr
    SpeedScale = 1;
    SpeedUnits = 'm/s';
end
Xsc = Xc./SpaceScale;
Zsc = Zc./SpaceScale;

if Nx <= 1 && Nz <= 1  % create 0D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    subplot(4,1,1)
    plot(hist.time/TimeScale,hist.T(:,2)-273.15,CL{[1,2]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/TimeScale,hist.Tliq(:,2),CL{[1,3]},LW{:});
    plot(hist.time/TimeScale,hist.Tsol(:,2),CL{[1,4]},LW{:});
    title('$T [^\circ$C]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    for i=1:cal.ncmp
        pcmp(i) = plot(hist.time/TimeScale,squeeze(hist.cm(:,2,i)).*100,'-',LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis xy tight; box on; hold on
    end
    title('melt cmps [wt\%]',TX{:},FS{:}); legend(pcmp,cal.cmpStr(1:end),TX{:},FS{1},8,'Location','northwest'); set(gca,TL{:},TS{:});
    subplot(4,1,3)
    for i=1:cal.ncmp
        plot(hist.time/TimeScale,squeeze(hist.cx(:,2,i)).*100,'-',LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis xy tight; box on; hold on
    end
    title('xtal cmps [wt\%]',TX{:},FS{:});
    subplot(4,1,4)
    plot(hist.time/TimeScale,hist.mu (:,2)*100.*(hist.mu (:,2)>1e-9),CL{[1,3]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/TimeScale,hist.chi(:,2)*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:});
    plot(hist.time/TimeScale,hist.phi(:,2)*100.*(hist.phi(:,2)>1e-9),CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end
    subplot(4,1,1)
    plot(hist.time/TimeScale,hist.rhom(:,2),'-',CL{[1,3]},LW{:}); axis xy tight; box on; hold on
    plot(hist.time/TimeScale,hist.rhox(:,2),'-',CL{[1,4]},LW{:});
    plot(hist.time/TimeScale,hist.rho (:,2),'-',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    semilogy(hist.time/TimeScale,hist.eta (:,2),'k-',LW{:}); axis xy tight; box on; hold on
    semilogy(hist.time/TimeScale,hist.etam(:,2),'r-',LW{:});
    title('$\eta^m,\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,3)
    plot(hist.time/TimeScale,    hist.Gx(:,2)./hist.rho(:,2)*hr*100.*(hist.chi(:,2)>1e-9),CL{[1,4]},LW{:}); axis xy tight; box on; hold on;
    plot(hist.time/TimeScale,10.*hist.Gf(:,2)./hist.rho(:,2)*hr*100.*(hist.phi(:,2)>1e-9),CL{[1,5]},LW{:});
    title('$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt \%/hr]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,4)
    plot(hist.time/TimeScale,hist.dV(:,2)*100*TimeScale,'k-',LW{:}); axis xy tight; box on;
    title(['$\dot{V}$ [\%/',TimeUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});

    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end
    subplot(4,1,1)
    for i=1:cal.noxd
        plot(hist.time/TimeScale,squeeze(hist.c_oxd(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis xy tight; box on; hold on
    end
    title('Bulk oxds [wt\%]',TX{:},FS{:}); legend(cal.oxdStr,TX{:},FS{1},8,'Location','northwest'); set(gca,TL{:},TS{:});
    subplot(4,1,2)
    for i=1:cal.noxd
        plot(hist.time/TimeScale,squeeze(hist.cm_oxd(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis xy tight; box on; hold on
    end
    title('Melt oxds [wt\%]',TX{:},FS{:});
    subplot(4,1,3)
    for i=1:cal.noxd
        plot(hist.time/TimeScale,squeeze(hist.cx_oxd(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis xy tight; box on; hold on
    end
    title('Xtal oxds [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(4,1,4)
    for i=1:cal.nmsy
        plot(hist.time/TimeScale,squeeze(hist.cx_msy(:,2,i)),'-',LW{:},'color',ocean(round((i-1)*213/cal.nmsy)+1,:)); axis xy tight; box on; hold on
    end
    title('Xtal msys [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:}); legend(cal.msyStr,TX{:},FS{1},8,'Location','northwest'); set(gca,TL{:},TS{:});
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});

elseif Nx <= 1  % create 1D plots

    if ~exist('fh1','var'); fh1 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh1); clf;
    end
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    plot(T-273.15,Zsc.',CL{[1,2]},LW{:}); axis ij tight; box on; hold on;
    plot(cal.Tliq,Zsc.',CL{[1,3]},LW{:});
    plot(cal.Tsol,Zsc.',CL{[1,4]},LW{:});
    title('$T [^\circ$C]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    for i=1:cal.ncmp
        plot(squeeze(c(:,:,i)).*100,Zsc.',LW{:},'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); axis ij tight; box on; hold on;
    end
    title('Bulk cmps [wt\%]',TX{:},FS{:}); legend(cal.cmpStr,TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    plot(mu *100.*(mu >1e-9),Zsc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    plot(chi*100.*(chi>1e-9),Zsc.',CL{[1,4]},LW{:});
    plot(phi*100.*(phi>1e-9),Zsc.',CL{[1,5]},LW{:});
    title('$\mu$, $\chi$, $\phi$ [vol\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,4)
    plot(-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1)/SpeedScale,Zf.',CL{[1,5]},LW{:}); axis ij tight; box on; hold on;
    plot(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)/SpeedScale,Zf.',CL{[1,4]},LW{:});
    plot(-(mu ([1,1:end],:)+mu ([1:end,end],:))/2.*wm(:,2:end-1)/SpeedScale,Zf.',CL{[1,3]},LW{:});
    plot(-                                         W (:,2:end-1)/SpeedScale,Zf.',CL{[1,2]},LW{:});
    title(['$W$, $w_\Delta^f$, $w_\Delta^x$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});

    if ~exist('fh2','var'); fh2 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh2); clf;
    end 
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,5,1)
    plot(rhox,Zsc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(rhom,Zsc.',CL{[1,3]},LW{:});
    plot(rho ,Zsc.',CL{[1,2]},LW{:});
    title('$\bar{\rho}$ [kg/m$^3$]',TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,2)
    semilogx(min(eta,etam),Zsc.',CL{[1,3]},LW{:}); axis ij tight; box on; hold on;
    semilogx(eta,Zsc.',CL{[1,2]},LW{:});
    title('$\bar{\eta}$ [Pas]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,3)
    plot(    Gx./rho*100*hr.*(chi>1e-9),Zsc.',CL{[1,4]},LW{:}); axis ij tight; box on; hold on;
    plot(10.*Gf./rho*100*hr.*(phi>1e-9),Zsc.',CL{[1,5]},LW{:});
    title(['$10 \times \Gamma_f/\bar{\rho}$, $\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,4)
    plot(VolSrc*hr,Zsc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title(['$\dot{V}$ [1/hr]'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,5,5)
    plot(P(2:end-1,2:end-1),Zsc.',CL{[1,2]},LW{:}); axis ij tight; box on;
    title('$P$ [Pa]',TX{:},FS{:}); set(gca,TL{:},TS{:});
 
    if ~exist('fh3','var'); fh3 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh3); clf;
    end 
    sgtitle(['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k');
    subplot(1,4,1)
    for i=1:cal.noxd
        plot(squeeze( c_oxd(:,:,i)),Zsc.',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis ij tight; box on; hold on;
    end
    title('Bulk oxds [wt\%]',TX{:},FS{:});ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,2)
    for i=1:cal.noxd
        plot(squeeze(cm_oxd(:,:,i)),Zsc.',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis ij tight; box on; hold on;
    end
    title('Melt oxds [wt\%]',TX{:},FS{:}); set(gca,TL{:},TS{:});
    subplot(1,4,3)
    for i=1:cal.noxd
        plot(squeeze(cx_oxd(:,:,i)),Zsc.',LW{:},'color',ocean(round((i-1)*213/cal.noxd)+1,:)); axis ij tight; box on; hold on;
    end
    title('Xtal oxds [wt\%]',TX{:},FS{:}); legend(cal.oxdStr,TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});
    subplot(1,4,4)
    for i=1:cal.nmsy
        plot(squeeze(cx_msy(:,:,i)),Zsc.',LW{:},'color',ocean(round((i-1)*213/cal.nmsy)+1,:)); axis ij tight; box on; hold on;
    end
    title('Xtal msys [wt\%]',TX{:},FS{:}); legend(cal.msyStr,TX{:},FS{:},'Location','west'); set(gca,TL{:},TS{:});

else % create 2D plots

    % set axis and border dimensions
    axh = 6.00*sqrt(D/L); axw = 6.00*sqrt(L/D)+1.50;
    ahs = 0.60; avs = 0.80;
    axb = 1.20; axt = 1.50;
    axl = 1.50; axr = 0.50;

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
    fh = axb + 3*axh + 2*avs + axt;
    fw = axl + 3*axw + 2*ahs + axr;
    set(fh8,UN{:},'Position',[11 11 fw fh]);
    set(fh8,'PaperUnits','Centimeters','PaperPosition',[0 0 fw fh],'PaperSize',[fw fh]);
    set(fh8,'Color','w','InvertHardcopy','off','Resize','off');
    ax(81) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+2*axh+2*avs axw axh]);
    ax(82) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+2*axh+2*avs axw axh]);
    ax(83) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+2*axh+2*avs axw axh]);
    ax(84) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+1*axh+1*avs axw axh]);
    ax(85) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+1*axh+1*avs axw axh]);
    ax(86) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+1*axh+1*avs axw axh]);
    ax(87) = axes(UN{:},'position',[axl+0*axw+0*ahs axb+0*axh+0*avs axw axh]);
    ax(88) = axes(UN{:},'position',[axl+1*axw+1*ahs axb+0*axh+0*avs axw axh]);
    ax(89) = axes(UN{:},'position',[axl+2*axw+2*ahs axb+0*axh+0*avs axw axh]);

    % plot velocity-pressure solution in Fig. 1
    set(0,'CurrentFigure',fh1)
    set(fh1,'CurrentAxes',ax(11));
    imagesc(Xsc,Zsc,-W(:      ,2:end-1)./SpeedScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh1,'CurrentAxes',ax(12));
    imagesc(Xsc,Zsc, U(2:end-1,:      )./SpeedScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$U$ [',SpeedUnits,']'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh1,'CurrentAxes',ax(13));
    imagesc(Xsc,Zsc, P); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$P$ [Pa]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh1,'CurrentAxes',ax(14));
    imagesc(Xsc,Zsc,Div_V*TimeScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\nabla \cdot \mathbf{v}$ [1/',TimeUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot temperature and composition in Fig. 2
    set(0,'CurrentFigure',fh2)
    set(fh2,'CurrentAxes',ax(21));
    imagesc(Xsc,Zsc,Tp-273.15); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T_p$ [$^\circ$C]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh2,'CurrentAxes',ax(22));
    imagesc(Xsc,Zsc,squeeze(c_oxd(:,:,cal.Si)./sum(c_oxd(:,:,1:end-1),3).*100)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['SiO$_2$ [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    text(0.5,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh2,'CurrentAxes',ax(23));
    if any(c(:,:,end),'all')
        imagesc(Xsc,Zsc,squeeze(c_oxd(:,:,cal.H))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['H$_2$O [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    else
        imagesc(Xsc,Zsc,squeeze(c_oxd(:,:,cal.Fe)./sum(c_oxd(:,:,1:end-1),3).*100)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['FeO [wt\%]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    end

    % plot phase fractions and reaction rates in Fig. 3
    set(0,'CurrentFigure',fh3)
    set(fh3,'CurrentAxes',ax(31));
    imagesc(Xsc,Zsc,chi.*100.*(chi>1e-9) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh3,'CurrentAxes',ax(32));
    imagesc(Xsc,Zsc,phi.*100.*(phi>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\phi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh3,'CurrentAxes',ax(33));
    imagesc(Xsc,Zsc,Gx./rho*hr*100.*(chi>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_x/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh3,'CurrentAxes',ax(34));
    imagesc(Xsc,Zsc,Gf./rho*hr*100.*(phi>1e-9)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Gamma_f/\bar{\rho}$ [wt\%/hr]'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot density, rheology, and segregation speeds in Fig. 4
    set(0,'CurrentFigure',fh4)
    set(fh4,'CurrentAxes',ax(41));
    imagesc(Xsc,Zsc,rho-mean(rho,2)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\Delta \bar{\rho}_h$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh4,'CurrentAxes',ax(42));
    imagesc(Xsc,Zsc,log10(eta)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{\eta}$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(-0.1,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh4,'CurrentAxes',ax(43));
    imagesc(Xsc,Zsc,-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1).*1e3/SpeedScale); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m',SpeedUnits,']'],TX{:},FS{:}); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh4,'CurrentAxes',ax(44));
    if any(phi(:)>1e-9)
        imagesc(Xsc,Zsc,-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1).*1e3/SpeedScale.*(any(phi>1e-9))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^f$ [m',SpeedUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    else
        imagesc(Xsc,Zsc,-(mu ([1,1:end],:)+mu ([1:end,end],:))/2.*wm(:,2:end-1).*1e3/SpeedScale.*(any(mu>1e-9))); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^m$ [m',SpeedUnits,']'],TX{:},FS{:}); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:}); set(gca,'YTickLabel',[]);
    end

    % plot pseudo-component composition in Fig. 5
    set(0,'CurrentFigure',fh5)
    sumanh = sum(c(:,:,1:end-1),3);
    for i = 1:cal.ncmp-1
        set(fh5,'CurrentAxes',ax(50+i));
        if i<cal.ncmp
            imagesc(Xsc,Zsc,c(:,:,i)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
        else
            imagesc(Xsc,Zsc,c(:,:,i).*100); axis ij equal tight; box on; cb = colorbar;
        end
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.cmpStr{i},' [wt\%]'],TX{:},FS{:});
    end
    set(fh5,'CurrentAxes',ax(51));
    set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh5,'CurrentAxes',ax(52));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(0.5,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh5,'CurrentAxes',ax(53));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh5,'CurrentAxes',ax(54));
    ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(55));
    set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh5,'CurrentAxes',ax(56));
    set(gca,'YTickLabel',[]);

    % plot major oxide composition in Fig. 6
    set(0,'CurrentFigure',fh6)
    sumanh = sum(c_oxd(:,:,1:end-1),3);
    for i = 1:cal.noxd
        set(fh6,'CurrentAxes',ax(60+i));
        if i<cal.noxd
            imagesc(Xsc,Zsc,c_oxd(:,:,i)./sumanh.*100); axis ij equal tight; box on; cb = colorbar;
        else
            imagesc(Xsc,Zsc,c_oxd(:,:,i)); axis ij equal tight; box on; cb = colorbar;
        end
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.oxdStr{i},' [wt\%]'],TX{:},FS{:});
    end
    set(fh6,'CurrentAxes',ax(61));
    set(gca,'XTickLabel',[]);
    set(fh6,'CurrentAxes',ax(62));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(0.5,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh6,'CurrentAxes',ax(63));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(64));
    set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(65));
    set(gca,'XTickLabel',[],'YTickLabel',[]); 
    set(fh6,'CurrentAxes',ax(66));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh6,'CurrentAxes',ax(67));
    set(fh6,'CurrentAxes',ax(68));
    set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh6,'CurrentAxes',ax(69));
    set(gca,'YTickLabel',[]);

    % plot mineral assemblage in Fig. 7
    set(0,'CurrentFigure',fh7)
    for i = 1:cal.nmsy
        set(fh7,'CurrentAxes',ax(70+i));
        imagesc(Xsc,Zsc,cx_msy(:,:,i)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.msyStr{i},' [wt\%]'],TX{:},FS{:});
    end
    set(fh7,'CurrentAxes',ax(71));
    set(gca,TL{:},TS{:}); set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:}); 
    set(fh7,'CurrentAxes',ax(72));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(0.5,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh7,'CurrentAxes',ax(73));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh7,'CurrentAxes',ax(74));
    ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
    set(fh7,'CurrentAxes',ax(75));
    set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh7,'CurrentAxes',ax(76));
    set(gca,'YTickLabel',[]);

    % plot geochemical variables in Fig. 7
    set(0,'CurrentFigure',fh8)
    for i = 1:cal.ntrc
        set(fh8,'CurrentAxes',ax(80+i));
        imagesc(Xsc,Zsc,trc(:,:,i)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title([cal.trcStr{i},' [1]'],TX{:},FS{:});
    end
    set(fh8,'CurrentAxes',ax(81));
    set(gca,'XTickLabel',[]);
    set(fh8,'CurrentAxes',ax(82));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    text(0.5,1.1,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
    set(fh8,'CurrentAxes',ax(83));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh8,'CurrentAxes',ax(84));
    set(gca,'XTickLabel',[]); ylabel(['Depth [',SpaceUnits,']'],TX{:},FS{:});
    set(fh8,'CurrentAxes',ax(85));
    set(gca,'XTickLabel',[],'YTickLabel',[]); 
    set(fh8,'CurrentAxes',ax(86));
    set(gca,'XTickLabel',[],'YTickLabel',[]);
    set(fh8,'CurrentAxes',ax(87));
    set(fh8,'CurrentAxes',ax(88));
    set(gca,'YTickLabel',[]); xlabel(['Width [',SpaceUnits,']'],TX{:},FS{:});
    set(fh8,'CurrentAxes',ax(89));
    set(gca,'YTickLabel',[]);

end


% plot T-X phase diagram in component space
if ~exist('fh9','var'); fh9 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh9);
end
if Nz>1; clf; end
for i=1:min(6,cal.ncmp)
subplot(3,2,i)
if i<cal.ncmp
    plot( c(:,:,i)./sum( c(:,:,1:end-1),3)*100,T-273.15,'k.','LineWidth',2,'MarkerSize',15); axis tight; box on; hold on;
    plot(cx(:,:,i)./sum(cx(:,:,1:end-1),3)*100,T-273.15,'b.','LineWidth',2,'MarkerSize',15);
    plot(cm(:,:,i)./sum(cm(:,:,1:end-1),3)*100,T-273.15,'r.','LineWidth',2,'MarkerSize',15);
else
    plot( c(:,:,i)*100,T-273.15,'k.','LineWidth',2,'MarkerSize',15); axis tight; box on; hold on;
    plot(cx(:,:,i)*100,T-273.15,'b.','LineWidth',2,'MarkerSize',15);
    plot(cm(:,:,i)*100,T-273.15,'r.','LineWidth',2,'MarkerSize',15);
end
xlabel([cal.cmpStr{i},' [wt \%]'],'Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',10)
end

if Nz>1
subplot(3,2,2)
text(-0.5,1.2,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
end
subplot(3,2,3)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)


% plot T-X phase diagram in oxide space
if ~exist('fh10','var'); fh10 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh10);
end
if Nz>1; clf; end
for i=1:cal.noxd
subplot(3,3,i)
if i<cal.ncmp
    plot( c_oxd(:,:,i)./sum( c_oxd(:,:,1:end-1),3)*100,T-273.15,'k.','LineWidth',2,'MarkerSize',15); axis tight; box on; hold on;
    plot(cx_oxd(:,:,i)./sum(cx_oxd(:,:,1:end-1),3)*100,T-273.15,'b.','LineWidth',2,'MarkerSize',15);
    plot(cm_oxd(:,:,i)./sum(cm_oxd(:,:,1:end-1),3)*100,T-273.15,'r.','LineWidth',2,'MarkerSize',15);
else
    plot( c_oxd(:,:,i),T-273.15,'k.','LineWidth',2,'MarkerSize',15); axis tight; box on; hold on;
    plot(cx_oxd(:,:,i),T-273.15,'b.','LineWidth',2,'MarkerSize',15);
    plot(cm_oxd(:,:,i),T-273.15,'r.','LineWidth',2,'MarkerSize',15);
end
xlabel([cal.oxdStr{i},' [wt \%]'],'Interpreter','latex','FontSize',12)
set(gca,'TickLabelInterpreter','latex','FontSize',10)
end

if Nz>1
subplot(3,3,2)
text(0.5,1.2,['time = ',num2str(time/TimeScale,3),' [',TimeUnits,']'],TX{:},FS{:},'Color','k','HorizontalAlignment','center','Units','normalized');
end
subplot(3,3,4)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)


% plot composition on TAS, AFM diagrams
if ~exist('fh11','var'); fh11 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh11);
end
if Nz>1 || step==0 || frst; clf;
TAS; axis tight; box on; hold on;
end
cxSi = cx_oxd_all(:,:,1)./sum(cx_oxd_all(:,:,1:end-1),3).*100;
cmSi = cm_oxd_all(:,:,1)./sum(cm_oxd_all(:,:,1:end-1),3).*100;
 cSi =  c_oxd_all(:,:,1)./sum( c_oxd_all(:,:,1:end-1),3).*100;
cxNK = sum(cx_oxd_all(:,:,[7,8]),3)./sum(cx_oxd_all(:,:,1:end-1),3).*100;
cmNK = sum(cm_oxd_all(:,:,[7,8]),3)./sum(cm_oxd_all(:,:,1:end-1),3).*100;
 cNK = sum( c_oxd_all(:,:,[7,8]),3)./sum( c_oxd_all(:,:,1:end-1),3).*100;
scatter(cxSi(:),cxNK(:),50,T(:)-273.15,'filled','^'); colormap(ocean); cb = colorbar;
scatter(cmSi(:),cmNK(:),50,T(:)-273.15,'filled','o');
scatter( cSi(:), cNK(:),80,T(:)-273.15,'filled','s');
set(cb,TL{:},'FontSize',12); set(gca,TL{:},'FontSize',15); xlabel('SiO$_2$ [wt \%]',TX{:},'FontSize',15); ylabel('Na$_2$O + K$_2$O [wt \%]',TX{:},'FontSize',15);


if ~exist('fh12','var'); fh12 = figure(VIS{:});
else; set(0, 'CurrentFigure', fh12);
end
if Nz>1 || step==0 || frst; clf;
AFM; axis tight; box on; hold on;
end
[A,B] = terncoords(cx_oxd_all(:,:, 5      )./sum(cx_oxd_all(:,:,[5,4,7,8]),3), ...
                   cx_oxd_all(:,:, 4      )./sum(cx_oxd_all(:,:,[5,4,7,8]),3), ...
               sum(cx_oxd_all(:,:,[7,8]),3)./sum(cx_oxd_all(:,:,[5,4,7,8]),3));
scatter(A(:),B(:),50,T(:)-273.15,'filled','^'); colormap(ocean); cb = colorbar;
[A,B] = terncoords(cm_oxd_all(:,:, 5      )./sum(cm_oxd_all(:,:,[5,4,7,8]),3), ...
                   cm_oxd_all(:,:, 4      )./sum(cm_oxd_all(:,:,[5,4,7,8]),3), ...
               sum(cm_oxd_all(:,:,[7,8]),3)./sum(cm_oxd_all(:,:,[5,4,7,8]),3));
scatter(A(:),B(:),50,T(:)-273.15,'filled','o'); colormap(ocean);
[A,B] = terncoords(c_oxd_all(:,:, 5      )./(sum(c_oxd_all(:,:,[5,4,7,8]),3)), ...
                   c_oxd_all(:,:, 4      )./(sum(c_oxd_all(:,:,[5,4,7,8]),3)), ...
               sum(c_oxd_all(:,:,[7,8]),3)./(sum(c_oxd_all(:,:,[5,4,7,8]),3)));
scatter(A(:),B(:),80,T(:)-273.15,'filled','s'); colormap(ocean);
set(cb,TL{:},'FontSize',12); set(gca,TL{:},'FontSize',15); xlabel('SiO$_2$ [wt \%]',TX{:},'FontSize',15); ylabel('Na$_2$O + K$_2$O [wt \%]',TX{:},'FontSize',15);

% plot model history
if plot_cv
    if ~exist('fh13','var'); fh13 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh13); clf;
    end 
    subplot(3,1,1);
    plot(hist.time/TimeScale,hist.DM./hist.sumM,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,2);
    plot(hist.time/TimeScale,hist.DS./hist.sumS,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,3);
    plot(hist.time/TimeScale,hist.DC./hist.sumC,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('consv. $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});

    if ~exist('fh14','var'); fh14 = figure(VIS{:});
    else; set(0, 'CurrentFigure', fh14); clf;
    end 
    subplot(3,1,1);
    plot(hist.time/TimeScale,hist.EM,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $M$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,2);
    plot(hist.time/TimeScale,hist.ES,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $S$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    subplot(3,1,3);
    plot(hist.time/TimeScale,hist.EC,'k-',LW{:}); hold on; axis tight; box on;
    ylabel('error $C$',TX{:},FS{:}); set(gca,TL{:},TS{:},'XTickLabel',[]);
    xlabel(['Time [',TimeUnits,']'],TX{:},FS{:});
end

drawnow

% save output to file
if save_op && ~restart
    if Nx <= 10 && Nz <= 10  % print 0D plots
        name = [outdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_cmp',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
        print(fh11,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
        print(fh12,name,'-dpng','-r300','-image');
    elseif Nx <= 10  % create 1D plots
        name = [outdir,'/',runID,'/',runID,'_sol_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_aux_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_cmp_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
        print(fh11,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
        print(fh12,name,'-dpng','-r300','-image');
    else
        name = [outdir,'/',runID,'/',runID,'_vep_',num2str(floor(step/nop))];
        print(fh1,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_tch_',num2str(floor(step/nop))];
        print(fh2,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
        print(fh3,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
        print(fh4,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_cmp',num2str(floor(step/nop))];
        print(fh5,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_oxd',num2str(floor(step/nop))];
        print(fh6,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_mnr',num2str(floor(step/nop))];
        print(fh7,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_gch',num2str(floor(step/nop))];
        print(fh8,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_TAS',num2str(floor(step/nop))];
        print(fh11,name,'-dpng','-r300','-image');
        name = [outdir,'/',runID,'/',runID,'_AFM',num2str(floor(step/nop))];
        print(fh12,name,'-dpng','-r300','-image');
    end

    name = [outdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','Pchmb','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','Tp','c','cm','cx','cf','TRC','trc','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTRCdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
    name = [outdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','Pchmb','f','x','m','fq','xq','mq','phi','chi','mu','X','F','M','S','C','T','Tp','c','cm','cx','cf','TRC','trc','dSdt','dCdt','dFdt','dXdt','dMdt','drhodt','dTRCdt','Gf','Gx','Gm','rho','eta','eII','tII','dt','time','step','VolSrc','wf','wx','wm','cal');
    name = [outdir,'/',runID,'/',runID,'_hist'];
    save(name,'hist');

end

if save_op && (step==0 || restart)
    logfile = [outdir,'/',runID,'/',runID,'.log'];
    if exist(logfile,'file') && step==0; delete(logfile); end
    diary(logfile)
end
    