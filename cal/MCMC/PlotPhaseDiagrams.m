
figure(200); clf;

clear cal;
cal_ASVZ;
var.c      = [linspace(0,1,100).',linspace(1,0,100).',zeros(1,100).',zeros(1,100).',zeros(1,100).',ones(1,100).'*0];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = 1200*ones(1,100).';          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,cala]  = meltmodel(var,cal,'T');

calh.Tliq = cala.Tliq;
for i=1:5
var.c      = [linspace(0,1,100).',linspace(1,0,100).',zeros(1,100).',zeros(1,100).',zeros(1,100).',ones(1,100).'*0.05];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = calh.Tliq;          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,calh]   = meltmodel(var,cal,'T');
end

subplot(2,2,1)
plot(var.c(:,2),cala.Tsol,'b-'); hold on; axis tight; box on;
plot(var.c(:,2),cala.Tliq,'r-');
plot(var.c(:,2),calh.Tsol,'b--');
plot(var.c(:,2),calh.Tliq,'r--');

%%
clear cal;
cal_ASVZ;
var.c      = [zeros(1,100).',linspace(0,1,100).',linspace(1,0,100).',zeros(1,100).',zeros(1,100).',ones(1,100).'*0];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = 1200*ones(1,100).';          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,cala]  = meltmodel(var,cal,'T');

calh.Tliq = cala.Tliq;
for i=1:5
var.c      = [zeros(1,100).',linspace(0,1,100).',linspace(1,0,100).',zeros(1,100).',zeros(1,100).',ones(1,100).'*0.05];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = calh.Tliq;          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,calh]   = meltmodel(var,cal,'T');
end

subplot(2,2,2)
plot(var.c(:,3),cala.Tsol,'b-'); hold on; axis tight; box on;
plot(var.c(:,3),cala.Tliq,'r-');
plot(var.c(:,3),calh.Tsol,'b--');
plot(var.c(:,3),calh.Tliq,'r--');

%%
clear cal;
cal_ASVZ;
var.c      = [zeros(1,100).',zeros(1,100).',linspace(0,1,100).',linspace(1,0,100).',zeros(1,100).',ones(1,100).'*0];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = 1200*ones(1,100).';          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,cala]  = meltmodel(var,cal,'T');

calh.Tliq = cala.Tliq;
for i=1:5
var.c      = [zeros(1,100).',zeros(1,100).',linspace(0,1,100).',linspace(1,0,100).',zeros(1,100).',ones(1,100).'*0.05];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = calh.Tliq;          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,calh]   = meltmodel(var,cal,'T');
end

subplot(2,2,3)
plot(var.c(:,4),cala.Tsol,'b-'); hold on; axis tight; box on;
plot(var.c(:,4),cala.Tliq,'r-');
plot(var.c(:,4),calh.Tsol,'b--');
plot(var.c(:,4),calh.Tliq,'r--');

%%
clear cal;
cal_ASVZ;
var.c      = [zeros(1,100).',zeros(1,100).',zeros(1,100).',linspace(0,1,100).',linspace(1,0,100).',ones(1,100).'*0];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = 1200*ones(1,100).';          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,cala]  = meltmodel(var,cal,'T');

calh.Tliq = cala.Tliq;
for i=1:5
var.c      = [zeros(1,100).',zeros(1,100).',zeros(1,100).',linspace(0,1,100).',linspace(1,0,100).',ones(1,100).'*0.05];          % component fractions [wt]
var.c      = var.c./sum(var.c,2);
var.T      = calh.Tliq;          % temperature [C]
var.P      = 1.5e8*ones(1,100).'/1e9;      % pressure [GPa]
var.H2O    = var.c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P,0,cal);
[~,calh]   = meltmodel(var,cal,'T');
end

subplot(2,2,4)
plot(var.c(:,5),cala.Tsol,'b-'); hold on; axis tight; box on;
plot(var.c(:,5),cala.Tliq,'r-');
plot(var.c(:,5),calh.Tsol,'b--');
plot(var.c(:,5),calh.Tliq,'r--');