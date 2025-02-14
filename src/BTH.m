% plot ternary axes
hold on;
HVS = {'handlevisibility','off'};
set(gca,'visible','off');

% prepare for plotting ternary plots
mm       = 5;
nn       = mm-1;
sin60    = sin(1/3*pi);
maxz     = 0;
plot3([0.5,1,0,0.5],[0,sin60,sin60,0],maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
plot3([0,1,0.5,0]-0.55,[0,0,sin60,0],maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
plot3([0,1,0.5,0]+0.55,[0,0,sin60,0],maxz.*ones(1,4),'color',[0,0,0],'LineWidth',2,HVS{:});
axis equal tight;

grids    = linspace(0,1,mm+1);
grids    = grids(1:end-1);
labels   = num2str(grids(2:end)');

[x3, y3] = terncoords(1-grids, grids, zeros(size(grids)));
[x2, y2] = terncoords(grids, zeros(size(grids)), 1-grids);
[x1, y1] = terncoords(zeros(size(grids)), 1-grids, grids);

% middle ternary
x1m = x1;    y1m = sin60 - y1;
x2m = x2;    y2m = sin60 - y2;
x3m = x3;    y3m = sin60 - y3;

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
text(x2m(2:end)-0.00,y2m(2:end)+0.07,labels,'FontSize',15,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');
% text(x1m(2:end)-0.03,y1m(2:end)-0.03,labels,'FontSize',15,TX{:},'HorizontalAlignment','right','VerticalAlignment','bottom');


% --- Left ternary grid and labels ---

% transform grid coordinates for left ternary
x1l = x1 - 0.55;    y1l = y1;
x2l = x2 - 0.55;    y2l = y2;
x3l = x3 - 0.55;    y3l = y3;

% plot grid lines for left ternary
for i = 1:nn
    plot3([x1l(i+1), x2l(nn-i+2)], [y1l(i+1), y2l(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x2l(i+1), x3l(nn-i+2)], [y2l(i+1), y3l(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x3l(i+1), x1l(nn-i+2)], [y3l(i+1), y1l(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x1l(i+1), x2l(nn-i+2)], [y1l(i+1), y2l(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x2l(i+1), x3l(nn-i+2)], [y2l(i+1), y3l(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x3l(i+1), x1l(nn-i+2)], [y3l(i+1), y1l(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
end

text(x2l(2:end)-0.00,y2l(2:end)-0.02,labels,'FontSize',16,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');
text(x1l(2:end)-0.02,y1l(2:end)-0.02,labels,'FontSize',16,TX{:},'HorizontalAlignment','right','VerticalAlignment','bottom');


% --- Right ternary grid and labels ---

% transform grid coordinates for right ternary
x1r = x1 + 0.55;    y1r = y1;
x2r = x2 + 0.55;    y2r = y2;
x3r = x3 + 0.55;    y3r = y3;

% plot grid lines for right ternary
for i = 1:nn
    plot3([x1r(i+1), x2r(nn-i+2)], [y1r(i+1), y2r(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x2r(i+1), x3r(nn-i+2)], [y2r(i+1), y3r(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x3r(i+1), x1r(nn-i+2)], [y3r(i+1), y1r(nn-i+2)], maxz*ones(2,1),'w',HVS{:});
    plot3([x1r(i+1), x2r(nn-i+2)], [y1r(i+1), y2r(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x2r(i+1), x3r(nn-i+2)], [y2r(i+1), y3r(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
    plot3([x3r(i+1), x1r(nn-i+2)], [y3r(i+1), y1r(nn-i+2)], maxz*ones(2,1),'k:',HVS{:});
end

text(x3r(2:end)+0.02,y3r(2:end)-0.02,labels,'FontSize',16,TX{:},'HorizontalAlignment','left','VerticalAlignment','bottom');
text(x2r(2:end)-0.00,y2r(2:end)-0.02,labels,'FontSize',16,TX{:},'HorizontalAlignment','center','VerticalAlignment','top');

% vertex labels
text(1.00,0.92,LBL1,'FontSize',18,TX{:},'HorizontalAlignment','left','VerticalAlignment','middle');
text(0.50,-0.09,LBL2,'FontSize',18,TX{:},'HorizontalAlignment','center','VerticalAlignment','bottom');
text(0.00,0.92,LBL3,'FontSize',18,TX{:},'HorizontalAlignment','right','VerticalAlignment','middle');
text(0.50-1.1,-0.09,LBL4,'FontSize',18,TX{:},'HorizontalAlignment','center','VerticalAlignment','bottom');
text(0.50+1.1,-0.09,LBL4,'FontSize',18,TX{:},'HorizontalAlignment','center','VerticalAlignment','bottom');
