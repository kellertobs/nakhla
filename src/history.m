% record run history

% record model time
hist.time(step+1) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(step+1) = sum(sum(rho(2:end-1,2:end-1)*h*h*1));  % [kg]
hist.sumH(step+1) = sum(sum(  H(2:end-1,2:end-1)*h*h*1));  % [J]
hist.sumC(step+1) = sum(sum(  C(2:end-1,2:end-1)*h*h*1));  % [kg]
hist.sumV(step+1) = sum(sum(  V(2:end-1,2:end-1)*h*h*1));  % [kg]

if step==0
    hist.sumM0 = hist.sumM(1);
    hist.sumH0 = hist.sumH(1);
    hist.sumC0 = hist.sumC(1);
    hist.sumV0 = hist.sumV(1);
end

% record expected rates of change by volume change and imposed boundaries layers
dsumMdt = sum(rho(2,2:end-1).*W(1,2:end-1)*h*1) - sum(rho(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(rho(2:end-1,2).*U(2:end-1,1)*h*1) - sum(rho(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [kg/s]
dsumHdt = sum(sum(bndH(2:end-1,2:end-1)*h*h*1)) ...
        + sum(  H(2,2:end-1).*W(1,2:end-1)*h*1) - sum(  H(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(  H(2:end-1,2).*U(2:end-1,1)*h*1) - sum(  H(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [J /s]
dsumCdt = sum(sum(bndC(2:end-1,2:end-1)*h*h*1)) ...
        + sum(  C(2,2:end-1).*W(1,2:end-1)*h*1) - sum(  C(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(  C(2:end-1,2).*U(2:end-1,1)*h*1) - sum(  C(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [kg/s]
dsumVdt = sum(sum(bndV(2:end-1,2:end-1)*h*h*1)) ...
        + sum(  V(2,2:end-1).*W(1,2:end-1)*h*1) - sum(  V(end-1,2:end-1).*W(end,2:end-1)*h*1) ...
        + sum(  V(2:end-1,2).*U(2:end-1,1)*h*1) - sum(  V(2:end-1,end-1).*U(2:end-1,end)*h*1);  % [kg/s]

if step>0; hist.DM(step+1) = hist.DM(step) + dsumMdt.*dt; else; hist.DM(step+1) = 0; end  % [kg]
if step>0; hist.DH(step+1) = hist.DH(step) + dsumHdt.*dt; else; hist.DH(step+1) = 0; end  % [J ]
if step>0; hist.DC(step+1) = hist.DC(step) + dsumCdt.*dt; else; hist.DC(step+1) = 0; end  % [kg]
if step>0; hist.DV(step+1) = hist.DV(step) + dsumVdt.*dt; else; hist.DV(step+1) = 0; end  % [kg]

% record conservation error of mass M, heat H, major component C, volatile component V
hist.EM(step+1) = (hist.sumM(step+1) - hist.DM(step+1))./hist.sumM0 - 1;  % [kg/kg]
hist.EH(step+1) = (hist.sumH(step+1) - hist.DH(step+1))./hist.sumH0 - 1;  % [J /J ]
hist.EC(step+1) = (hist.sumC(step+1) - hist.DC(step+1))./hist.sumC0 - 1;  % [kg/kg]
hist.EV(step+1) = (hist.sumV(step+1) - hist.DV(step+1))./hist.sumV0 - 1;  % [kg/kg]

% record variable and coefficient diagnostics
hist.W(step+1,1) = min(min(-W(:,2:end-1)));
hist.W(step+1,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(step+1,3) = max(max(-W(:,2:end-1)));

hist.U(step+1,1) = min(min(U(2:end-1,:)));
hist.U(step+1,2) = mean(mean(abs(U(2:end-1,:))));
hist.U(step+1,3) = max(max(U(2:end-1,:)));

hist.P(step+1,1) = min(min(P(2:end-1,2:end-1)));
hist.P(step+1,2) = mean(mean(abs(P(2:end-1,2:end-1))));
hist.P(step+1,3) = max(max(P(2:end-1,2:end-1)));

hist.x(step+1,1) = min(min(x(2:end-1,2:end-1)));
hist.x(step+1,2) = mean(mean(x(2:end-1,2:end-1)));
hist.x(step+1,3) = max(max(x(2:end-1,2:end-1)));

hist.f(step+1,1) = min(min(f(2:end-1,2:end-1)));
hist.f(step+1,2) = mean(mean(f(2:end-1,2:end-1)));
hist.f(step+1,3) = max(max(f(2:end-1,2:end-1)));

hist.m(step+1,1) = min(min(m(2:end-1,2:end-1)));
hist.m(step+1,2) = mean(mean(m(2:end-1,2:end-1)));
hist.m(step+1,3) = max(max(m(2:end-1,2:end-1)));

hist.chi(step+1,1) = min(min(chi(2:end-1,2:end-1)));
hist.chi(step+1,2) = mean(mean(chi(2:end-1,2:end-1)));
hist.chi(step+1,3) = max(max(chi(2:end-1,2:end-1)));

hist.phi(step+1,1) = min(min(phi(2:end-1,2:end-1)));
hist.phi(step+1,2) = mean(mean(phi(2:end-1,2:end-1)));
hist.phi(step+1,3) = max(max(phi(2:end-1,2:end-1)));

hist.mu(step+1,1) = min(min(mu(2:end-1,2:end-1)));
hist.mu(step+1,2) = mean(mean(mu(2:end-1,2:end-1)));
hist.mu(step+1,3) = max(max(mu(2:end-1,2:end-1)));

hist.T(step+1,1) = min(min(T(2:end-1,2:end-1)));
hist.T(step+1,2) = mean(mean(T(2:end-1,2:end-1)));
hist.T(step+1,3) = max(max(T(2:end-1,2:end-1)));

hist.c(step+1,1) = min(min(c(2:end-1,2:end-1)));
hist.c(step+1,2) = mean(mean(c(2:end-1,2:end-1)));
hist.c(step+1,3) = max(max(c(2:end-1,2:end-1)));

hist.v(step+1,1) = min(min(v(2:end-1,2:end-1)));
hist.v(step+1,2) = mean(mean(v(2:end-1,2:end-1)));
hist.v(step+1,3) = max(max(v(2:end-1,2:end-1)));

hist.rho(step+1,1) = min(min(rho(2:end-1,2:end-1)));
hist.rho(step+1,2) = mean(mean(rho(2:end-1,2:end-1)));
hist.rho(step+1,3) = max(max(rho(2:end-1,2:end-1)));

hist.eta(step+1,1) = min(min(eta(2:end-1,2:end-1)));
hist.eta(step+1,2) = geomean(geomean(eta(2:end-1,2:end-1)));
hist.eta(step+1,3) = max(max(eta(2:end-1,2:end-1)));

hist.wx(step+1,1) = min(min(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)));
hist.wx(step+1,2) = mean(mean(abs((chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1))));
hist.wx(step+1,3) = max(max(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)));

hist.wf(step+1,1) = min(min(-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1)));
hist.wf(step+1,2) = mean(mean(abs((phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1))));
hist.wf(step+1,3) = max(max(-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1)));
