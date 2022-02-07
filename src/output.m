
if plot_op
    % prepare for plotting
    TX = {'Interpreter','Latex'}; FS = {'FontSize',12};
    TL = {'TickLabelInterpreter','Latex'}; TS = {'FontSize',10};
    UN = {'Units','Centimeters'};
    
    % set axis and border dimensions
    axh = 6.00; axw = axh*L/D;
    ahs = 0.40; avs = 0.2;
    axb = 1.20; axt = 0.4;
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
    axes(ax(11));
    imagesc(X(2:end-1),Z(2:end-1),-W(:      ,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$W$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    text(L/2,0.9*D,['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k','VerticalAlignment','middle','HorizontalAlignment','center');
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
    axes(ax(21));
    imagesc(X(2:end-1),Z(2:end-1),T(2:end-1,2:end-1)     ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$T [^\circ$C]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
    axes(ax(22));
    imagesc(X(2:end-1),Z(2:end-1),c(2:end-1,2:end-1).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{c}$ [wt\% SiO$_2$]'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(23));
    imagesc(X(2:end-1),Z(2:end-1),v(2:end-1,2:end-1).*100); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\bar{v}$ [wt\% H$_2$O]'],TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot phase fractions and reaction rates in Fig. 3
    figure(3);
    axes(ax(31));
    imagesc(X(2:end-1),Z(2:end-1),chi(2:end-1,2:end-1).*100.*(chi(2:end-1,2:end-1)>1e-9) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\chi$ [vol\%]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    text(L/2,0.9*D,['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','w','VerticalAlignment','middle','HorizontalAlignment','center');
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
    axes(ax(41));
    imagesc(X(2:end-1),Z(2:end-1),      rho(2:end-1,2:end-1) ); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\rho$ [kg/m$^3$]'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    text(L/2,0.9*D,['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','w','VerticalAlignment','middle','HorizontalAlignment','center');
    axes(ax(42));
    imagesc(X(2:end-1),Z(2:end-1),log10(eta(2:end-1,2:end-1))); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$\eta$ [log$_{10}$ Pas]'],TX{:},FS{:}); set(gca,'XTickLabel',[],'YTickLabel',[]);
    axes(ax(43));
    imagesc(X(2:end-1),Z(2:end-1),-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^x$ [m/hr]'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:});
    axes(ax(44));
    imagesc(X(2:end-1),Z(2:end-1),-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1).*hr); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['$w_\Delta^f$ [m/hr]'],TX{:},FS{:}); xlabel('Width [m]',TX{:},FS{:}); set(gca,'YTickLabel',[]);

    % plot geochemical variables in Fig. 5
    figure(5);
    axes(ax(51));
    imagesc(X(2:end-1),Z(2:end-1),it(2:end-1,2:end-1)); axis ij equal tight; box on; cb = colorbar;
    set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['incomp. trace'],TX{:},FS{:}); set(gca,'XTickLabel',[]); ylabel('Depth [m]',TX{:},FS{:});
    text(L/2,0.9*D,['time = ',num2str(time/hr,3),' [hr]'],TX{:},FS{:},'Color','k','VerticalAlignment','middle','HorizontalAlignment','center');
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
        imagesc(X(2:end-1),Z(2:end-1),-res_W(:      ,2:end-1)./(1e-16+norm(RR(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $W$'],TX{:},FS{:}); ylabel('Depth [m]',TX{:},FS{:});
        axes(ax(62));
        imagesc(X(2:end-1),Z(2:end-1), res_U(2:end-1,:      )./(1e-16+norm(RR(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $U$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); xlabel('Width [m]',TX{:},FS{:});
        axes(ax(63));
        imagesc(X(2:end-1),Z(2:end-1), res_P(2:end-1,2:end-1)./(1e-16+norm(RR(:),2)./N)); axis ij equal tight; box on; cb = colorbar;
        set(cb,TL{:},TS{:}); set(gca,TL{:},TS{:}); title(['res. $P$'],TX{:},FS{:}); set(gca,'YTickLabel',[]); 
    end
end

% plot phase diagram
fh7 = figure(7); clf;
vv = (4.8e-5.*Ptop.^0.6 + 1e-9.*Ptop)./100;
TT = linspace(Tphs0+Ptop*clap,Tphs1+Ptop*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,ceil((perT-Tphs0)./(Tphs1-Tphs0)*1e3)),linspace((perCx+perCm)/2,cphs0,floor((perT-Tphs1)./(Tphs0-Tphs1)*1e3))];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
plot(CCx,TT,'k:','LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k:','LineWidth',2);
Tphs0s = Tphs0-dTH2O(1)*v0^0.75;
Tphs1s = Tphs1-dTH2O(3)*v0^0.75;
perTs  = perT-dTH2O(2)*v0^0.75;
TT = linspace(Tphs0s+Ptop*clap,Tphs1s+Ptop*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,round((perTs-Tphs0s)./(Tphs1s-Tphs0s)*1e3)),linspace((perCx+perCm)/2,cphs0,round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*1e3))];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,v0*ones(size(TT)),Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
plot(CCx,TT,'k-','LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k-','LineWidth',2);
Tphs0s = Tphs0-dTH2O(1)*vv^0.75;
Tphs1s = Tphs1-dTH2O(3)*vv^0.75;
perTs  = perT-dTH2O(2)*vv^0.75;
TT = linspace(Tphs0s+Ptop*clap,Tphs1s+Ptop*clap,1e3);
cc = [linspace(cphs1,(perCx+perCm)/2,round((perTs-Tphs0s)./(Tphs1s-Tphs0s)*1e3)),linspace((perCx+perCm)/2,cphs0,round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*1e3))];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,vv*ones(size(TT)),Ptop*ones(size(TT)),Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
plot(CCx,TT,'k:','LineWidth',2); axis tight; hold on; box on;
plot(CCm,TT,'k:','LineWidth',2);

Tplt = T - (Pt-Ptop)*clap;
cplt = c./(1-f);
plot(cplt(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'k.',cx(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'b.',cm(2:end-1,2:end-1),Tplt(2:end-1,2:end-1),'r.','LineWidth',2,'MarkerSize',15);

set(gca,'TickLabelInterpreter','latex','FontSize',15)
title('Phase Diagram','Interpreter','latex','FontSize',22)
xlabel('Composition','Interpreter','latex','FontSize',18)
ylabel('Temperature','Interpreter','latex','FontSize',18)
        
% plot model history
fh8 = figure(8);
if step>0
subplot(3,1,1);
plot(hist.time(1:nop:end)./3600,hist.T(1:nop:end,2),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('$T$ [$^\circ$C]',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(3,1,2);
histV = sqrt(hist.W.^2 + hist.U.^2);
plot(hist.time(1:nop:end)./3600,log10(histV(1:nop:end,2)),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('log$_{10}$ $|\mathbf{v}|$ [m/s]',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(3,1,3);
plot(hist.time(1:nop:end)./3600,hist.chi(1:nop:end,2),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
plot(hist.time(1:nop:end)./3600,hist.phi(1:nop:end,2),'bo','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
plot(hist.time(1:nop:end)./3600,hist.mu (1:nop:end,2),'ro','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
xlabel('Time [hr]',TX{:},FS{:});
ylabel('$\mu$, $\phi$, $\chi$ [vol]',TX{:},FS{:});
end

fh9 = figure(9); clf;
if step>0
subplot(4,1,1);
plot(hist.time(1:nop:end)./3600,hist.DM(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('consv. $M$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(4,1,2);
plot(hist.time(1:nop:end)./3600,hist.DH(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('consv. $H$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(4,1,3);
plot(hist.time(1:nop:end)./3600,hist.DC(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('consv. $C$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(4,1,4);
plot(hist.time(1:nop:end)./3600,hist.DV(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('consv. $V$',TX{:},FS{:});
xlabel('Time [hr]',TX{:},FS{:});
end

fh10 = figure(10); clf;
if step>0
subplot(4,1,1);
plot(hist.time(1:nop:end)./3600,hist.EM(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('error $M$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(4,1,2);
plot(hist.time(1:nop:end)./3600,hist.EH(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('error $H$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(4,1,3);
plot(hist.time(1:nop:end)./3600,hist.EC(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('error $C$',TX{:},FS{:}); set(gca,'XTickLabel',[]);
subplot(4,1,4);
plot(hist.time(1:nop:end)./3600,hist.EV(1:nop:end),'ko','MarkerSize',5,'LineWidth',2); hold on; axis tight; box on;
ylabel('error $V$',TX{:},FS{:});
xlabel('Time [hr]',TX{:},FS{:});
end

drawnow

% save output to file
if save_op
    name = [opdir,'/',runID,'/',runID,'_vp_',num2str(floor(step/nop))];
    print(fh1,name,'-dpng','-r300','-opengl');
    name = [opdir,'/',runID,'/',runID,'_tc_',num2str(floor(step/nop))];
    print(fh2,name,'-dpng','-r300','-opengl');
    name = [opdir,'/',runID,'/',runID,'_phs_',num2str(floor(step/nop))];
    print(fh3,name,'-dpng','-r300','-opengl');
    name = [opdir,'/',runID,'/',runID,'_sgr_',num2str(floor(step/nop))];
    print(fh4,name,'-dpng','-r300','-opengl');
    name = [opdir,'/',runID,'/',runID,'_chm',num2str(floor(step/nop))];
    print(fh5,name,'-dpng','-r300','-opengl');
    name = [opdir,'/',runID,'/',runID,'_eql',num2str(floor(step/nop))];
    print(fh7,name,'-dpng','-r300','-opengl');
    
    name = [opdir,'/',runID,'/',runID,'_',num2str(floor(step/nop))];
    save(name,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SIm','SIx','SI','RIP','RID','it','ct','sim','six','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSImdt','dSIxdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist');
    name = [opdir,'/',runID,'/',runID,'_cont'];
    save(name,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SIm','SIx','SI','RIP','RID','it','ct','sim','six','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSImdt','dSIxdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist');
    
    if step == 0
        logfile = [opdir,'/',runID,'/',runID,'.log'];
        if exist(logfile,'file'); delete(logfile); end
        diary(logfile)
    end
end

clear fh1 fh2 fh3 fh4 fh5 fh6 fh7 fh8 fh9 fh10
    