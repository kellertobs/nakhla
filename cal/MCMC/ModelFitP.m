function [datafit,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,Tsolfit,Tliqfit,Tm] = ModelFitP(model,T,P,SYS,PHS,Psl,cal,reg)

% get number of data points
np = length(T);

% get model parameters from model vector
cal.T0   = model(               (1:cal.ncmp-1)).';
cal.A    = model(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.B    = model(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.r    = model(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.dTH2O= model(4*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem  = reshape(model(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd  = cmp_mem*cal.mem_oxd./100;

% find phase component compositions by non-negative least-squares
% MLT_cmp = zeros(np,cal.ncmp);
% for ip = 1:np
%     MLT_cmp(ip,:) = lsqnonneg(cmp_oxd.',MLT(ip,:).');
% end
% MLT_cmp(:,1:end-1) = MLT_cmp(:,1:end-1) + reg(1) ...
%                                         + reg(2).*(diff(MLT_cmp([1 1:end end],1:end-1),2,1)/8 + diff(MLT_cmp(:,[1 1:end-1 end-1]),2,2)/8);
% MLT_cmp = MLT_cmp./sum(MLT_cmp,2);
% 
% SOL_cmp = zeros(np,cal.ncmp);
% for ip = 1:np
%     SOL_cmp(ip,:) = lsqnonneg(cmp_mem.',SOL(ip,:).');
% end
% SOL_cmp(:,1:end-1) = SOL_cmp(:,1:end-1) + reg(1) ...
%                                         + reg(2).*(diff(SOL_cmp([1 1:end end],1:end-1),2,1)/8 + diff(SOL_cmp(:,[1 1:end-1 end-1]),2,2)/8);
% SOL_cmp = SOL_cmp./sum(SOL_cmp,2);
% 
% SYS_cmp = PHS(:,1)/100.*MLT_cmp + (1-PHS(:,1)/100).*SOL_cmp;

SYS_cmpfit = zeros(np,cal.ncmp);
for ip = 1:np
    SYS_cmpfit(ip,:) = lsqnonneg(cmp_oxd.',SYS(ip,:).');
end
SYS_cmpfit(:,1:end-1) = SYS_cmpfit(:,1:end-1) + reg(1) ...
                                        + reg(2).*(diff(SYS_cmpfit([1 1:end end],1:end-1),2,1)/8 + diff(SYS_cmpfit(:,[1 1:end-1 end-1]),2,2)/8);
SYS_cmpfit = SYS_cmpfit./sum(SYS_cmpfit,2);

% get fitted phase oxide compositions
SYS_oxdfit = SYS_cmpfit*cmp_oxd;

% get local phase equilibrium
var.m = ones(np,1)/2; var.x = var.m; var.f = 0*var.m;
var.c      = SYS_cmpfit;        % component fractions [wt]
var.X      = var.c*cal.cmp_oxd; % melt oxide fractions [wt %]
var.T      = T;             % temperature [C]
var.P      = P/10;         % pressure [GPa]
var.H2O    = SYS_cmpfit(:,end); % water concentration [wt]
cal.H2Osat = var.H2O./(PHS(:,1)/100);%fluidsat(var);
[var,cal]  = meltmodel(var,cal,'E');

% get fitted phase oxide compositions
MLT_cmpfit = var.cm;
SOL_cmpfit = var.cx;
MLT_oxdfit = var.cm*cmp_oxd;
SOL_oxdfit = var.cx*cmp_oxd;
SOL_memfit = var.cx*cmp_mem;

% update mineral systems oxide compositions for solid assemblage
PHS_oxdfit = zeros(np,cal.nmsy+1,cal.noxd);
PHS_oxdfit(:,1,:) = MLT_oxdfit;
for j = 1:cal.nmsy
    PHS_oxdfit(:,j+1,:) = SOL_memfit(:,cal.msy_mem(j,:)==1)*cal.mem_oxd(cal.msy_mem(j,:)==1,:)./sum(SOL_memfit(:,cal.msy_mem(j,:)==1)+1e-16,2);
end
PHS_oxdfit = PHS_oxdfit.*(PHS>0);

% get fitted phase fractions
PHS_frcfit = zeros(size(PHS));
PHS_frcfit(:,1) = var.m*100;
PHS_frcfit(:,2:end) = SOL_memfit*cal.msy_mem.';%.*(1-PHS_frcfit(:,1)/100);

% get solidus and liquidus at starting composition
cal = rmfield(cal,{'Tsol' 'Tliq'});
Psl = Psl(1:end/2);
var.m = ones(size(Psl))/2; var.x = var.m; var.f = 0*var.m;
var.c      = repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(Psl),1);
var.P      = Psl/10;                     % pressure [GPa]
var.T      = 1000+Psl(:)*5e-8;           % temperature [C]
var.H2O    = Psl(:)*0;                   % water concentration [wt]
cal.H2Osat = Psl(:)*0;
[~,cal]    = meltmodel(var,cal,'T');

Tsolfit = cal.Tsol;
Tliqfit = cal.Tliq;
Tm      = cal.Tm;

% get solidus and liquidus at end composition
cal = rmfield(cal,{'Tsol' 'Tliq'});
var.m = ones(size(Psl))/2; var.x = var.m; var.f = 0*var.m;
var.c      = repmat(mean(SYS_cmpfit(end-4:end,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(end-4:end,1:end-1),2),1),length(Psl),1);
var.P      = Psl(:)/10;                  % pressure [GPa]
var.T      = 1000+Psl(:)*5e-8;           % temperature [C]
var.H2O    = Psl(:)*0;                   % water concentration [wt]
cal.H2Osat = Psl(:)*0;
[~,cal]    = meltmodel(var,cal,'T');

Tsolfit = [Tsolfit;cal.Tsol];
Tliqfit = [Tliqfit;cal.Tliq];

datafit = [MLT_oxdfit(:);SOL_memfit(:);PHS_frcfit(:);Tsolfit(:);Tliqfit(:)];%repmat(PHSfit(:,1),6,1)];

end