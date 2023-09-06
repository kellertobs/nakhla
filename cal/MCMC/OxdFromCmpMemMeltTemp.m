function [datafit] = OxdFromCmpMemMeltTemp(model,T,P,H2Osat,MLT,SOL,SYS,m,cal)

% get number of data points
np = length(SYS);

% reshape component mineral endmember compositions to ncmp x nmem
cmp_mem = reshape(model(1:cal.ncmp*cal.nmem),cal.ncmp,cal.nmem);
nn      = cal.ncmp*cal.nmem;

% get component oxide compositions
cal.cmp_oxd = cmp_mem*cal.mem_oxd./100;
cal.T0      = model(nn+           (1:cal.ncmp-1)).';
cal.r       = model(nn+cal.ncmp-1+(1:cal.ncmp-1)).';
cal.A       = (cal.T0+273.15)./300;

% unpack data
% SYS = reshape(data(1:np*cal.noxd),np,cal.noxd);
% MLT = reshape(data(1:np*cal.noxd),np,cal.noxd);
% SOL = reshape(data(  np*cal.noxd+1:2*np*cal.noxd),np,cal.noxd);
% PHS = reshape(data(2*np*cal.noxd+1:end),np,cal.nmsy+1);

% cmpSYS = zeros(np,cal.ncmp);
% for ip = 1:np
%     cmpSYS(ip,:) = lsqnonneg(cal.cmp_oxd.',SYS(ip,:).');
% end
% cmpSYS = cmpSYS./sum(cmpSYS,2);

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

cmpSYS = m/100.*cmpMLT + (1-m/100).*cmpSOL;

var.m = m/100; var.x = 0*var.m; var.f = 0*var.m;

% update local phase equilibrium
% var.c      = cmpSYS;        % component fractions [wt]
% var.T      = T;             % temperature [C]
% var.P      = P/1e9;         % pressure [GPa]
% var.H2O    = cmpSYS(:,end); % water concentration [wt]
% cal.H2Osat = H2Osat/100;
% [var,cal]  = meltmodel(var,cal,'E');

% get fitted phase oxide compositions
% MLTfit = var.cm*cal.cmp_oxd;
% SOLfit = var.cx*cal.cmp_oxd;

% get fitted phase oxide compositions
% SYSfit = cmpSYS*cmp_oxd;
MLTfit = cmpMLT*cal.cmp_oxd;
SOLfit = cmpSOL*cal.cmp_oxd;

% get fitted phase fractions
PHSfit = zeros(np,cal.nmsy+1);
% PHSfit(:,1) = var.m*100;
for ip = 1:np
    PHSfit(ip,1) = lsqnonneg((MLTfit(ip,:)-SOLfit(ip,:)).',(SYS(ip,:)-SOLfit(ip,:)).')*100;
end

% get mineral systems proportions in solid phase
memSOL          = cmpSOL*cmp_mem;
PHSfit(:,2:end) = memSOL*cal.msy_mem.'.*(1-PHSfit(:,1)/100);

datafit = [MLTfit(:);SOLfit(:);PHSfit(:)];

end