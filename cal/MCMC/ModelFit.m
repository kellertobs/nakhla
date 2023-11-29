function [datafit,MLTfit,SOLfit,SYSfit,PHSfit,cmpSYS] = ModelFit(model,T,P,MLT,SOL,M,cal)

% get number of data points
np = length(T);

% get model parameters from model vector
cal.T0      = model(           (1:cal.ncmp-1)).';
cal.r       = model(cal.ncmp-1+(1:cal.ncmp-1)).';
cal.A       = (cal.T0+273.15)./350;
cmp_mem     = reshape(model(2*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cal.cmp_oxd = cmp_mem*cal.mem_oxd./100;


% find phase component compositions by non-negative least-squares
cmpMLT = zeros(np,cal.ncmp);
for ip = 1:np
    cmpMLT(ip,:) = lsqnonneg(cal.cmp_oxd.',MLT(ip,:).');
end
cmpMLT = cmpMLT./sum(cmpMLT,2);

cmpSOL = zeros(np,cal.ncmp);
for ip = 1:np
    cmpSOL(ip,:) = lsqnonneg(cal.cmp_oxd.',SOL(ip,:).');
end
cmpSOL = cmpSOL./sum(cmpSOL,2);

cmpSYS = M/100.*cmpMLT + (1-M/100).*cmpSOL;

% % get system composition in component fractions
% cmpSYS = zeros(np,cal.ncmp);
% for ip = 1:np
%     cmpSYS(ip,:) = lsqnonneg(cal.cmp_oxd.',SYS(ip,:).');
% end
% cmpSYS = cmpSYS./sum(cmpSYS,2);

% get fitted phase oxide compositions
SYSfit = cmpSYS*cal.cmp_oxd;
% MLTfit = cmpMLT*cal.cmp_oxd;
% SOLfit = cmpSOL*cal.cmp_oxd;

% get local phase equilibrium
var.m = ones(np,1); var.x = 0*var.m; var.f = 0*var.m;
var.c      = cmpSYS;        % component fractions [wt]
var.T      = T;             % temperature [C]
var.P      = P/1e9;         % pressure [GPa]
var.H2O    = cmpSYS(:,end); % water concentration [wt]
cal.H2Osat = cmpSYS(:,end)+0.01;
[var,cal]  = meltmodel(var,cal,'E');

% get fitted phase oxide compositions
MLTfit = var.cm*cal.cmp_oxd;
SOLfit = var.cx*cal.cmp_oxd;
% SYSfit = var.m.*MLTfit + (1-var.m).*SOLfit;

% get fitted phase fractions
PHSfit = zeros(np,cal.nmsy+1);
PHSfit(:,1) = var.m*100;
PHSfit(:,2:end) = (var.cx*cmp_mem)*cal.msy_mem.';%.*(1-PHSfit(:,1)/100);

datafit = [MLTfit(:);SOLfit(:);0.2*PHSfit(:)];%repmat(PHSfit(:,1),6,1)];

end