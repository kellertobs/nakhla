function [oxdfit] = OxdFromCmpMem(cmp_mem,oxd0,T,cal)

cmp_mem = reshape(cmp_mem,cal.ncmp,cal.nmem);

cmp_oxd = cmp_mem*cal.mem_oxd./100;

Xp = zeros(length(oxd0),cal.ncmp);
for ip = 1:length(oxd0)
    Xp(ip,:) = lsqnonneg(cmp_oxd.',oxd0(ip,:).');
end
Xp = Xp./sum(Xp,2);

oxdfit = Xp*cmp_oxd;
oxdfit = oxdfit(:);

end