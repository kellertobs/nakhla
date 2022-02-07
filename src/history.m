% record run history

% record model time
hist.time(step+1) = time;

% record total mass, heat, component mass in model
hist.sumM(step+1) = sum(sum(rho(2:end-1,2:end-1)));
hist.sumH(step+1) = sum(sum(H(2:end-1,2:end-1)));
hist.sumC(step+1) = sum(sum(C(2:end-1,2:end-1)));
hist.sumV(step+1) = sum(sum(V(2:end-1,2:end-1)));

if step==0
    hist.sumM0 = hist.sumM(1);
    hist.sumH0 = hist.sumH(1);
    hist.sumC0 = hist.sumC(1);
    hist.sumV0 = hist.sumV(1);
end

% record total mass change by mixture density change in fixed volume domain
dsumMdt = mean(mean(Div_V(2:end-1,2:end-1)));
if step>0; hist.DM(step+1) = hist.DM(step) + dsumMdt.*dt; else; hist.DM(step+1) = 0; end

% record heat/mass change by imposed boundary layers
dsumHdt = sum(sum(bndH(2:end-1,2:end-1)))./hist.sumM(step+1);
dsumCdt = sum(sum(bndC(2:end-1,2:end-1)))./hist.sumM(step+1);
dsumVdt = sum(sum(bndV(2:end-1,2:end-1)))./hist.sumM(step+1);

if step>0; hist.DH(step+1) = hist.DH(step) + dsumHdt.*dt; else; hist.DH(step+1) = 0; end
if step>0; hist.DC(step+1) = hist.DC(step) + dsumCdt.*dt; else; hist.DC(step+1) = 0; end
if step>0; hist.DV(step+1) = hist.DV(step) + dsumVdt.*dt; else; hist.DV(step+1) = 0; end

% record conservation error of mass M, heat H, major component C, volatile component V
hist.EM(step+1) =  hist.sumM(step+1)./hist.sumM0-1+hist.DM(step+1);
hist.EH(step+1) = (hist.sumH(step+1)./hist.sumM(step+1)-hist.DH(step+1))./(hist.sumH0./hist.sumM0)-1;
hist.EC(step+1) = (hist.sumC(step+1)./hist.sumM(step+1)-hist.DC(step+1))./(hist.sumC0./hist.sumM0)-1;
hist.EV(step+1) = (hist.sumV(step+1)./hist.sumM(step+1)-hist.DV(step+1))./(hist.sumV0./hist.sumM0)-1;

% record variable and coefficient diagnostics
hist.W(step+1,1) = min(min(W(:,2:end-1)));
hist.W(step+1,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(step+1,3) = max(max(W(:,2:end-1)));

hist.U(step+1,1) = min(min(U(:,2:end-1)));
hist.U(step+1,2) = mean(mean(abs(U(:,2:end-1))));
hist.U(step+1,3) = max(max(U(:,2:end-1)));

hist.P(step+1,1) = min(min(P(:,2:end-1)));
hist.P(step+1,2) = mean(mean(abs(P(:,2:end-1))));
hist.P(step+1,3) = max(max(P(:,2:end-1)));

hist.x(step+1,1) = min(min(x(:,2:end-1)));
hist.x(step+1,2) = mean(mean(x(:,2:end-1)));
hist.x(step+1,3) = max(max(x(:,2:end-1)));

hist.f(step+1,1) = min(min(f(:,2:end-1)));
hist.f(step+1,2) = mean(mean(f(:,2:end-1)));
hist.f(step+1,3) = max(max(f(:,2:end-1)));

hist.m(step+1,1) = min(min(m(:,2:end-1)));
hist.m(step+1,2) = mean(mean(m(:,2:end-1)));
hist.m(step+1,3) = max(max(m(:,2:end-1)));

hist.chi(step+1,1) = min(min(chi(:,2:end-1)));
hist.chi(step+1,2) = mean(mean(chi(:,2:end-1)));
hist.chi(step+1,3) = max(max(chi(:,2:end-1)));

hist.phi(step+1,1) = min(min(phi(:,2:end-1)));
hist.phi(step+1,2) = mean(mean(phi(:,2:end-1)));
hist.phi(step+1,3) = max(max(phi(:,2:end-1)));

hist.mu(step+1,1) = min(min(mu(:,2:end-1)));
hist.mu(step+1,2) = mean(mean(mu(:,2:end-1)));
hist.mu(step+1,3) = max(max(mu(:,2:end-1)));

hist.T(step+1,1) = min(min(T(:,2:end-1)));
hist.T(step+1,2) = mean(mean(T(:,2:end-1)));
hist.T(step+1,3) = max(max(T(:,2:end-1)));

hist.c(step+1,1) = min(min(c(:,2:end-1)));
hist.c(step+1,2) = mean(mean(c(:,2:end-1)));
hist.c(step+1,3) = max(max(c(:,2:end-1)));

hist.v(step+1,1) = min(min(v(:,2:end-1)));
hist.v(step+1,2) = mean(mean(v(:,2:end-1)));
hist.v(step+1,3) = max(max(v(:,2:end-1)));

hist.rho(step+1,1) = min(min(rho(:,2:end-1)));
hist.rho(step+1,2) = mean(mean(rho(:,2:end-1)));
hist.rho(step+1,3) = max(max(rho(:,2:end-1)));

hist.eta(step+1,1) = min(min(eta(:,2:end-1)));
hist.eta(step+1,2) = geomean(geomean(eta(:,2:end-1)));
hist.eta(step+1,3) = max(max(eta(:,2:end-1)));
