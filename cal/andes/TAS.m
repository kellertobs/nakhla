CL = {'color',[0.5,0.5,0.5]};

% Define points of composition delineators
p1=  [41,1];
p2=  [45,1];
p3=  [41,7];
p4=  [45,5];
p5=  [52,5];
p6=  [57,5.9];
p7=  [63,7];
p8=  [69,8];
p9=  [69,13];
p10= [61,13.5];
p11= [57.6,11.7];
p12= [53,9.3];
p13= [49.4,7.3];
p14= [45,9.4];
p15= [48.4,11.5];
p16= [52.5,14];
p17= [78,1];
p18= [48.065, 16];

% Define lines
l1 = [p1; p3; p14; p16];
l2 = [p2; p4; p13; p10];
l3 = [p4; p5; p8];
l4 = [p18; p16; p11; p7; p7(1),1];
l5 = [p15; p12; p6; p6(1), 1];
l6 = [p9; p8; p17];
l7 = [p14; p13; p5; p5(1), 1];

% Draw lines
line( l1(:,1), l1(:,2), CL{:},LW{:});
line( l2(:,1), l2(:,2), CL{:},LW{:});
line( l3(:,1), l3(:,2), CL{:},LW{:});
line( l4(:,1), l4(:,2), CL{:},LW{:});
line( l5(:,1), l5(:,2), CL{:},LW{:});
line( l6(:,1), l6(:,2), CL{:},LW{:});
line( l7(:,1), l7(:,2), CL{:},LW{:});

% Annotate the sections
cLabels = [0.2, 0.2, 0.2];
text(42, 3, {'picro-';'basalt'}, CL{:},FS{:},TX{:});
text(44, 7.5, {'tephrite/';'basanite'}, CL{:},FS{:},TX{:});
text(48, 9.5, {'phono-';'tephrite'}, CL{:},FS{:},TX{:});
text(51.5, 12, {'tephri-';'phonolite'}, CL{:},FS{:},TX{:});
text(46, 13.5, 'foidite', CL{:},FS{:},TX{:});
text(56, 13.5, 'phonolite', CL{:},FS{:},TX{:});
text(62, 10, {'trachyte/';'trachydacite'}, CL{:},FS{:},TX{:});
text(72, 11, 'rhyolite', CL{:},FS{:},TX{:});
text(64.5, 3.0, 'dacite', CL{:},FS{:},TX{:});
text(58, 3.0, 'andesite', CL{:},FS{:},TX{:});
text(52.5, 2.5, {'basaltic';'andesite'}, CL{:},FS{:},TX{:});
text(46, 2.5, 'basalt', CL{:},FS{:},TX{:});
text(47.5, 5.75, {'trachy-';'basalt'}, CL{:},FS{:},TX{:});
text(51, 7, {'basaltic';'trachy-';'andesite'}, CL{:},FS{:},TX{:});
text(56,9, {'trachy-';'andesite'}, CL{:},FS{:},TX{:});
