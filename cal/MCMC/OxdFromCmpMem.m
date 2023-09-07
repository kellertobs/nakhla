function [datafit] = OxdFromCmpMem(model,MLT,SOL,M,cal)

% get number of data points
np = length(MLT);
T  = (1:np).';

% reshape component mineral endmember compositions to ncmp x nmem
cmp_mem = reshape(model,cal.ncmp,cal.nmem);

% get component oxide compositions
cmp_oxd = cmp_mem*cal.mem_oxd./100;

% get phase component compositions by non-negative least-squares
cmpMLT = zeros(np,cal.ncmp);
for ip = 1:np
    cmpMLT(ip,:) = lsqnonneg(cmp_oxd.',MLT(ip,:).');
end
cmpMLT = cmpMLT./sum(cmpMLT,2);

cmpSOL = zeros(np,cal.ncmp);
for ip = 1:np
    cmpSOL(ip,:) = lsqnonneg(cmp_oxd.',SOL(ip,:).');
end
cmpSOL = cmpSOL./sum(cmpSOL,2);

% get fitted phase oxide compositions
% SYSfit = cmpSYS*cmp_oxd;
MLTfit = cmpMLT*cmp_oxd;
SOLfit = cmpSOL*cmp_oxd;

% get fitted melt fraction by non-negative least-squares
% PHSfit = zeros(np,cal.nmsy+1);
% for ip = 1:np
%     PHSfit(ip,1) = max(0,min(100,lsqnonneg((MLTfit(ip,:)-SOLfit(ip,:)).',(SYS(ip,:)-SOLfit(ip,:)).')*100));
% end
% 
% % get mineral systems proportions in solid phase
% memSOL          = cmpSOL*cmp_mem;
% PHSfit(:,2:end) = memSOL*cal.msy_mem.'.*(1-PHSfit(:,1)/100);

cmpSYS = cmpMLT.*M/100 + cmpSOL.*(1-M/100);
SYSfit = cmpSYS*cmp_oxd;

datafit = SYSfit;%[MLTfit(:);SOLfit(:)];%;PHSfit(:)];

end