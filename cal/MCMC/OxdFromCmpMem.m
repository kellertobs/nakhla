function [datafit] = OxdFromCmpMem(model,data,oxdSYS,cal)

% get number of data points
np = length(oxdSYS);

% reshape component mineral endmember compositions to ncmp x nmem
cmp_mem = reshape(model,cal.ncmp,cal.nmem);

% get component oxide compositions
cmp_oxd = cmp_mem*cal.mem_oxd./100;

% unpack data
oxdLIQ = reshape(data(1:np*cal.noxd),np,cal.noxd);
oxdSOL = reshape(data(  np*cal.noxd+1:2*np*cal.noxd),np,cal.noxd);
phs    = reshape(data(2*np*cal.noxd+1:end),np,cal.nmsy+1);

% find phase component compositions by non-negative least-squares
cmpLIQ = zeros(np,cal.ncmp);
for ip = 1:np
    cmpLIQ(ip,:) = lsqnonneg(cmp_oxd.',oxdLIQ(ip,:).');
end
cmpLIQ = cmpLIQ./sum(cmpLIQ,2);

cmpSOL = zeros(np,cal.ncmp);
for ip = 1:np
    cmpSOL(ip,:) = lsqnonneg(cmp_oxd.',oxdSOL(ip,:).');
end
cmpSOL = cmpSOL./sum(cmpSOL,2);

% get fitted phase oxide compositions
oxdLIQfit = cmpLIQ*cmp_oxd;
oxdSOLfit = cmpSOL*cmp_oxd;

% get fitted phase fractions
phsfit = zeros(np,cal.nmsy+1);
for ip = 1:np
    phsfit(ip,1) = lsqnonneg((oxdLIQfit(ip,:)-oxdSOLfit(ip,:)).',(oxdSYS(ip,:)-oxdSOLfit(ip,:)).')*100;
end

% get mineral systems proportions in solid phase
memSOL          = cmpSOL*cmp_mem;
phsfit(:,2:end) = memSOL*cal.msy_mem.'.*(1-phsfit(:,1)/100);

datafit = [oxdLIQfit(:);oxdSOLfit(:);phsfit(:)];

end