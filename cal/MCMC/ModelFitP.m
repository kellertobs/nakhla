function [datafit,MLTfit,SOLfit,SYSfit,PHSfit,cmpSYS,Tsolfit,Tliqfit,Tm] = ModelFitP(model,T,P,MLT,SOL,SYS,M,Psl,cal)

% get number of data points
np = length(T);

% get model parameters from model vector
cal.T0      = model(               (1:cal.ncmp-1)).';
% cal.A       = model(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.A       = (cal.T0+273.15)./350;
cal.B       = model(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cal.r       = model(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem     = reshape(model(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cal.cmp_oxd = cmp_mem*cal.mem_oxd./100;

% find phase component compositions by non-negative least-squares
cmpMLT = zeros(np,cal.ncmp);
for ip = 1:np
    cmpMLT(ip,:) = lsqnonneg(cal.cmp_oxd.',MLT(ip,:).');
end
cmpMLT(:,1:end-1) = cmpMLT(:,1:end-1) + diff(cmpMLT([1 1:end end],1:end-1),2,1)/8 + diff(cmpMLT(:,[1 1:end-1 end-1]),2,2)/8;
% for i=1:cal.ncmp-1; cmpMLT(i,max(1,i-2:end-1)) = max(0.05,cmpMLT(i,max(1,i-2:end-1))); end
cmpMLT = cmpMLT./sum(cmpMLT,2);

cmpSOL = zeros(np,cal.ncmp);
for ip = 1:np
    cmpSOL(ip,:) = lsqnonneg(cal.cmp_oxd.',SOL(ip,:).');
end
cmpSOL(:,1:end-1) = cmpSOL(:,1:end-1) + diff(cmpSOL([1 1:end end],1:end-1),2,1)/8 + diff(cmpSOL(:,[1 1:end-1 end-1]),2,2)/8;
cmpSOL = cmpSOL./sum(cmpSOL,2);

cmpSYS = M/100.*cmpMLT + (1-M/100).*cmpSOL;

% cmpSYS = zeros(np,cal.ncmp);
% for ip = 1:np
%     cmpSYS(ip,:) = lsqnonneg(cal.cmp_oxd.',SYS(ip,:).');
% end
% for i=1:cal.ncmp-1; cmpSYS(i,max(1,i-2:end-1)) = max(0.05,cmpSYS(i,max(1,i-2:end-1))); end
% cmpSYS = cmpSYS./sum(cmpSYS,2);


% get fitted phase oxide compositions
SYSfit = cmpSYS*cal.cmp_oxd;

% get local phase equilibrium
var.m = ones(np,1); var.x = 0*var.m; var.f = 0*var.m;
var.c      = cmpSYS;        % component fractions [wt]
var.T      = T;             % temperature [C]
var.P      = P/1e9;         % pressure [GPa]
var.H2O    = cmpSYS(:,end); % water concentration [wt]
cal.H2Osat = cmpMLT(:,end)+1e-6;
[var,cal]  = meltmodel(var,cal,'E');

% get fitted phase oxide compositions
MLTfit = var.cm*cal.cmp_oxd;
SOLfit = var.cx*cal.cmp_oxd;

% get fitted phase fractions
PHSfit = zeros(np,cal.nmsy+1);
PHSfit(:,1) = var.m*100;
PHSfit(:,2:end) = (var.cx*cmp_mem)*cal.msy_mem.';%.*(1-PHSfit(:,1)/100);

% get solidus and liquidus at starting composition
cal = rmfield(cal,{'Tsol' 'Tliq'});
var.m = ones(size(Psl)); var.x = 0*var.m; var.f = 0*var.m;
var.c      = repmat(cmpSYS(1,:),length(Psl),1);   % component fractions [wt]
var.P      = Psl;         % pressure [GPa]
var.T      = 1000+Psl*1e9*5e-8;             % temperature [C]
var.H2O    = var.c(1,end)+Psl*0; % water concentration [wt]
cal.H2Osat = cmpMLT(1,end)+1e-6+Psl*0;
[~,cal]    = meltmodel(var,cal,'T');

Tsolfit = cal.Tsol;
Tliqfit = cal.Tliq;
Tm      = cal.Tm;

datafit = [MLTfit(:);SOLfit(:);PHSfit(:);Tsolfit(:);Tliqfit(:)];%repmat(PHSfit(:,1),6,1)];

end