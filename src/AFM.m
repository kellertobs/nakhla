% prepare for plotting ternary plots
mm = 5;
grids = linspace(0,1,mm+1);
grids = grids(1:end-1);
labels = num2str(grids(2:end)');
[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);
nn = mm-1;

% plot ternary axes
hold on;
set(gca,'visible','off');
maxz = 0;
plot3([0,1,0.5,0],[0,0,sin(1/3*pi),0],maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,'handlevisibility','off');
axis equal tight;

% plot grid on ternary axes
for i = 1:nn
    plot3([x1(i+1),x2(nn-i+2)],[y1(i+1),y2(nn-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x2(i+1),x3(nn-i+2)],[y2(i+1),y3(nn-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x3(i+1),x1(nn-i+2)],[y3(i+1),y1(nn-i+2)],maxz.*ones(2,1),'w','handlevisibility','off');
    plot3([x1(i+1),x2(nn-i+2)],[y1(i+1),y2(nn-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x2(i+1),x3(nn-i+2)],[y2(i+1),y3(nn-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
    plot3([x3(i+1),x1(nn-i+2)],[y3(i+1),y1(nn-i+2)],maxz.*ones(2,1),'k:','handlevisibility','off');
end

text(x3(2:end)+0.02,y3(2:end)-0.02,labels,'FontSize',16,'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom');
text(x2(2:end)-0.00,y2(2:end)-0.02,labels,'FontSize',16,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','top');
text(x1(2:end)-0.02,y1(2:end)-0.02,labels,'FontSize',16,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','bottom');
text( 0.99,-0.04,'MgO','FontSize',18,'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','middle');
text(0.450,0.87,'FeO$_\mathrm{tot}$','FontSize',18,'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom');
text(0.11,-0.04,'Na$_2$O+K$_2$O','FontSize',18,'Interpreter','latex','HorizontalAlignment','right','VerticalAlignment','middle');
text(1.0,0.87,'T [$^\circ$C]','FontSize',16,'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom');
text(0.50,0.60,'tholeiitic','Color',[60 5 15]/100,'FontSize',16,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
text(0.45,0.20,'calc-alcaline','Color',[60 5 15]/100,'FontSize',16,'Interpreter','latex','HorizontalAlignment','center','VerticalAlignment','middle');
[A,B] = terncoords([0.5,0.34,0.24,0.19,0.11,0.05],[0.36,0.49,0.56,0.56,0.46,0.35]);
plot(A,B,'--','Color',[60 5 15]/100,'LineWidth',3)
