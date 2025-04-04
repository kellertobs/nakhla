COL = {'color',[0.5,0.5,0.5]};

% Define points of composition delineators
p1=  [41,0];
p2=  [45,0];
p3=  [41,7];
p4=  [45,5];
p5=  [52,5];
p6=  [57,5.9];
p7=  [63,7];
p8=  [69,8];
p9=  [69,14];
p10= [62,14];
p11= [57.6,11.7];
p12= [53,9.3];
p13= [49.4,7.3];
p14= [45,9.4];
p15= [48.4,11.5];
p16= [52.5,14];
p17= [78,0];
p18= [41,3];
p19= [45,3];

% Define lines
l1 = [p1; p3; p14; p16];
l2 = [p2; p4; p13; p10];
l3 = [p4; p5; p8];
l4 = [p16; p11; p7; p7(1),0];
l5 = [p15; p12; p6; p6(1),0];
l6 = [p9; p8; p17];
l7 = [p14; p13; p5; p5(1),0];
l8 = [p18; p19];

% Draw lines
line( l1(:,1), l1(:,2), COL{:},LW{:});
line( l2(:,1), l2(:,2), COL{:},LW{:});
line( l3(:,1), l3(:,2), COL{:},LW{:});
line( l4(:,1), l4(:,2), COL{:},LW{:});
line( l5(:,1), l5(:,2), COL{:},LW{:});
line( l6(:,1), l6(:,2), COL{:},LW{:});
line( l7(:,1), l7(:,2), COL{:},LW{:});
line( l8(:,1), l8(:,2), COL{:},LW{:});

% Annotate the sections
cLabels = [0.2, 0.2, 0.2];
text(41.5, 2.0, {'picro-';'basalt'}, COL{:},FS{:},TX{:});
text(43.4, 7.1, {'tephrite/';'basanite'}, COL{:},FS{:},TX{:});
text(47.5, 9.2, {'phono-';'tephrite'}, COL{:},FS{:},TX{:});
text(51.0, 11.5, {'tephri-';'phonolite'}, COL{:},FS{:},TX{:});
text(45, 13, 'foidite', COL{:},FS{:},TX{:});
text(55, 13.4, 'phonolite', COL{:},FS{:},TX{:});
text(61.6, 10, {'trachyte/';'trachydacite'}, COL{:},FS{:},TX{:});
text(71.5, 10.5, 'rhyolite', COL{:},FS{:},TX{:});
text(65, 3.5, 'dacite', COL{:},FS{:},TX{:});
text(58, 3.0, 'andesite', COL{:},FS{:},TX{:});
text(52.5, 2.5, {'basaltic';'andesite'}, COL{:},FS{:},TX{:});
text(46.5, 2.5, 'basalt', COL{:},FS{:},TX{:});
text(47.5, 5.6, {'trachy-';'basalt'}, COL{:},FS{:},TX{:});
text(51.1, 6.8, {'basaltic';'trachy-';'andesite'}, COL{:},FS{:},TX{:});
text(55.6,8.6, {'trachy-';'andesite'}, COL{:},FS{:},TX{:});

text(77.5,14,'T [$^\circ$C]','FontSize',16,'Interpreter','latex','HorizontalAlignment','left','VerticalAlignment','bottom');
