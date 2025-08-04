% plot ternary axes
hold on;
HVS = {'handlevisibility','off'};
set(gca,'visible','off');

% prepare for plotting ternary plots
mm       = 5;
nn       = mm-1;
sin60    = sin(1/3*pi);
xshift   = 0.55;
zshift   = 2*xshift*sin60;
zshiftm  = 0.025/sin60;
maxz     = 0;
plot3([0,1,0.5,0],[0,0,sin60,0]+zshift,maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
plot3([0.5,1,0,0.5],[0,sin60,sin60,0]+zshiftm,maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
plot3([0,1,0.5,0]-xshift,[0,0,sin60,0],maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
plot3([0,1,0.5,0]+xshift,[0,0,sin60,0],maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
axis equal tight;

grids    = linspace(0,1,mm+1);
grids    = grids(1:end-1);
labels   = num2str(grids(2:end)');

[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);

% middle ternary
x1m = x1;    y1m = sin60 - y1 + zshiftm;
x2m = x2;    y2m = sin60 - y2 + zshiftm;
x3m = x3;    y3m = sin60 - y3 + zshiftm;

% plot grid on ternary axes
for i = 1:nn
    plot3([x1m(i+1),x2m(nn-i+2)],[y1m(i+1),y2m(nn-i+2)],maxz.*ones(2,1),'w',HVS{:});
    plot3([x2m(i+1),x3m(nn-i+2)],[y2m(i+1),y3m(nn-i+2)],maxz.*ones(2,1),'w',HVS{:});
    plot3([x3m(i+1),x1m(nn-i+2)],[y3m(i+1),y1m(nn-i+2)],maxz.*ones(2,1),'w',HVS{:});
    plot3([x1m(i+1),x2m(nn-i+2)],[y1m(i+1),y2m(nn-i+2)],maxz.*ones(2,1),'k:',HVS{:});
    plot3([x2m(i+1),x3m(nn-i+2)],[y2m(i+1),y3m(nn-i+2)],maxz.*ones(2,1),'k:',HVS{:});
    plot3([x3m(i+1),x1m(nn-i+2)],[y3m(i+1),y1m(nn-i+2)],maxz.*ones(2,1),'k:',HVS{:});
end

% text(x3m(2:end)+0.03,y3m(2:end)-0.03,labels,'FontSize',15,TX{:},'HorizontalAlignment','left','VerticalAlignment','bottom');
% text(x2m(2:end)-0.00,y2m(2:end)+0.07,labels,'FontSize',15,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');
% text(x1m(2:end)-0.03,y1m(2:end)-0.03,labels,'FontSize',15,TX{:},'HorizontalAlignment','right','VerticalAlignment','bottom');

% upper ternary
x1u = x1;    y1u = y1 + zshift;
x2u = x2;    y2u = y2 + zshift;
x3u = x3;    y3u = y3 + zshift;

% plot grid on ternary axes
for i = 1:nn
    plot3([x1u(i+1),x2u(nn-i+2)],[y1u(i+1),y2u(nn-i+2)],maxz.*ones(2,1),'w',HVS{:});
    plot3([x2u(i+1),x3u(nn-i+2)],[y2u(i+1),y3u(nn-i+2)],maxz.*ones(2,1),'w',HVS{:});
    plot3([x3u(i+1),x1u(nn-i+2)],[y3u(i+1),y1u(nn-i+2)],maxz.*ones(2,1),'w',HVS{:});
    plot3([x1u(i+1),x2u(nn-i+2)],[y1u(i+1),y2u(nn-i+2)],maxz.*ones(2,1),'k:',HVS{:});
    plot3([x2u(i+1),x3u(nn-i+2)],[y2u(i+1),y3u(nn-i+2)],maxz.*ones(2,1),'k:',HVS{:});
    plot3([x3u(i+1),x1u(nn-i+2)],[y3u(i+1),y1u(nn-i+2)],maxz.*ones(2,1),'k:',HVS{:});
end

text(x3u(2:end)+0.03,y3u(2:end)-0.03,labels,'FontSize',15,TX{:},'HorizontalAlignment','left','VerticalAlignment','bottom');
% text(x2u(2:end)-0.00,y2u(2:end)+0.07,labels,'FontSize',15,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');
text(x1u(2:end)-0.03,y1u(2:end)-0.03,labels,'FontSize',15,TX{:},'HorizontalAlignment','right','VerticalAlignment','bottom');


% transform grid coordinates for left ternary
x1l = x1 - xshift;    y1l = y1;
x2l = x2 - xshift;    y2l = y2;
x3l = x3 - xshift;    y3l = y3;

% plot grid lines for left ternary
for i = 1:nn
    plot3([x1l(i+1), x2l(nn-i+2)], [y1l(i+1), y2l(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x2l(i+1), x3l(nn-i+2)], [y2l(i+1), y3l(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x3l(i+1), x1l(nn-i+2)], [y3l(i+1), y1l(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x1l(i+1), x2l(nn-i+2)], [y1l(i+1), y2l(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x2l(i+1), x3l(nn-i+2)], [y2l(i+1), y3l(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x3l(i+1), x1l(nn-i+2)], [y3l(i+1), y1l(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
end

text(x2l(2:end)-0.00,y2l(2:end)-0.02,labels,'FontSize',15,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');
text(x1l(2:end)-0.02,y1l(2:end)-0.02,labels,'FontSize',15,TX{:},'HorizontalAlignment','right','VerticalAlignment','bottom');


% transform grid coordinates for right ternary
x1r = x1 + xshift;    y1r = y1;
x2r = x2 + xshift;    y2r = y2;
x3r = x3 + xshift;    y3r = y3;

% plot grid lines for right ternary
for i = 1:nn
    plot3([x1r(i+1), x2r(nn-i+2)], [y1r(i+1), y2r(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x2r(i+1), x3r(nn-i+2)], [y2r(i+1), y3r(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x3r(i+1), x1r(nn-i+2)], [y3r(i+1), y1r(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x1r(i+1), x2r(nn-i+2)], [y1r(i+1), y2r(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x2r(i+1), x3r(nn-i+2)], [y2r(i+1), y3r(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x3r(i+1), x1r(nn-i+2)], [y3r(i+1), y1r(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
end

text(x3r(2:end)+0.02,y3r(2:end)-0.02,labels,'FontSize',15,TX{:},'HorizontalAlignment','left','VerticalAlignment','bottom');
text(x2r(2:end)-0.00,y2r(2:end)-0.02,labels,'FontSize',15,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');

% vertex labels
text(1.05,0.94,LBL1,'FontSize',20,TX{:},'HorizontalAlignment','left','VerticalAlignment','middle');
text(0.50,-0.1,LBL2,'FontSize',20,TX{:},'HorizontalAlignment','center','VerticalAlignment','bottom');
text(-0.05,0.94,LBL3,'FontSize',20,TX{:},'HorizontalAlignment','right','VerticalAlignment','middle');
text(0.50-1.1,-0.1,LBL4,'FontSize',20,TX{:},'HorizontalAlignment','center','VerticalAlignment','bottom');
text(0.50+1.1,-0.1,LBL4,'FontSize',20,TX{:},'HorizontalAlignment','center','VerticalAlignment','bottom');
text(0.50,zshift+0.93,LBL4,'FontSize',20,TX{:},'HorizontalAlignment','center','VerticalAlignment','middle');
