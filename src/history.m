% record run history

dsumMdto = dsumMdt;
dsumHdto = dsumHdt;
dsumCdto = dsumCdt;
dsumVdto = dsumVdt;

stp = max(1,step);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(stp) = sum(sum(rho(2:end-1,2:end-1)*h*h*1))+TINY;  % [kg]
hist.sumH(stp) = sum(sum(  H(2:end-1,2:end-1)*h*h*1))+TINY;  % [J]
hist.sumC(stp) = sum(sum(  C(2:end-1,2:end-1)*h*h*1))+TINY;  % [kg]
hist.sumV(stp) = sum(sum(  V(2:end-1,2:end-1)*h*h*1))+TINY;  % [kg]

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

% if step>1; hist.DM(stp) = hist.DM(stp-1) + (THETA*dsumMdt + (1-THETA)*dsumMdto).*dt; else; hist.DM(stp) = 0; end  % [kg]
if step>1; hist.DM(stp) = hist.DM(stp-1) + dsumMdt.*dt; else; hist.DM(stp) = 0; end  % [kg]
if step>1; hist.DH(stp) = hist.DH(stp-1) + (THETA*dsumHdt + (1-THETA)*dsumHdto).*dt; else; hist.DH(stp) = 0; end  % [J ]
if step>1; hist.DC(stp) = hist.DC(stp-1) + (THETA*dsumCdt + (1-THETA)*dsumCdto).*dt; else; hist.DC(stp) = 0; end  % [kg]
if step>1; hist.DV(stp) = hist.DV(stp-1) + (THETA*dsumVdt + (1-THETA)*dsumVdto).*dt; else; hist.DV(stp) = 0; end  % [kg]

% record conservation error of mass M, heat H, major component C, volatile component V
hist.EM(stp) = (hist.sumM(stp) - hist.DM(stp))./hist.sumM(1) - 1;  % [kg/kg]
hist.EH(stp) = (hist.sumH(stp) - hist.DH(stp))./hist.sumH(1) - 1;  % [J /J ]
hist.EC(stp) = (hist.sumC(stp) - hist.DC(stp))./hist.sumC(1) - 1;  % [kg/kg]
hist.EV(stp) = (hist.sumV(stp) - hist.DV(stp))./hist.sumV(1) - 1;  % [kg/kg]

% record variable and coefficient diagnostics
hist.W(stp,1) = min(min(-W(:,2:end-1)));
hist.W(stp,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(stp,3) = max(max(-W(:,2:end-1)));

hist.U(stp,1) = min(min(U(2:end-1,:)));
hist.U(stp,2) = mean(mean(abs(U(2:end-1,:))));
hist.U(stp,3) = max(max(U(2:end-1,:)));

hist.P(stp,1) = min(min(P(2:end-1,2:end-1)));
hist.P(stp,2) = mean(mean(abs(P(2:end-1,2:end-1))));
hist.P(stp,3) = max(max(P(2:end-1,2:end-1)));

hist.x(stp,1) = min(min(x(2:end-1,2:end-1)));
hist.x(stp,2) = mean(mean(x(2:end-1,2:end-1)));
hist.x(stp,3) = max(max(x(2:end-1,2:end-1)));

hist.f(stp,1) = min(min(f(2:end-1,2:end-1)));
hist.f(stp,2) = mean(mean(f(2:end-1,2:end-1)));
hist.f(stp,3) = max(max(f(2:end-1,2:end-1)));

hist.m(stp,1) = min(min(m(2:end-1,2:end-1)));
hist.m(stp,2) = mean(mean(m(2:end-1,2:end-1)));
hist.m(stp,3) = max(max(m(2:end-1,2:end-1)));

hist.chi(stp,1) = min(min(chi(2:end-1,2:end-1)));
hist.chi(stp,2) = mean(mean(chi(2:end-1,2:end-1)));
hist.chi(stp,3) = max(max(chi(2:end-1,2:end-1)));

hist.phi(stp,1) = min(min(phi(2:end-1,2:end-1)));
hist.phi(stp,2) = mean(mean(phi(2:end-1,2:end-1)));
hist.phi(stp,3) = max(max(phi(2:end-1,2:end-1)));

hist.mu(stp,1) = min(min(mu(2:end-1,2:end-1)));
hist.mu(stp,2) = mean(mean(mu(2:end-1,2:end-1)));
hist.mu(stp,3) = max(max(mu(2:end-1,2:end-1)));

hist.T(stp,1) = min(min(T(2:end-1,2:end-1)));
hist.T(stp,2) = mean(mean(T(2:end-1,2:end-1)));
hist.T(stp,3) = max(max(T(2:end-1,2:end-1)));

hist.c(stp,1) = min(min(c(2:end-1,2:end-1)));
hist.c(stp,2) = mean(mean(c(2:end-1,2:end-1)));
hist.c(stp,3) = max(max(c(2:end-1,2:end-1)));

hist.v(stp,1) = min(min(v(2:end-1,2:end-1)));
hist.v(stp,2) = mean(mean(v(2:end-1,2:end-1)));
hist.v(stp,3) = max(max(v(2:end-1,2:end-1)));

hist.cx(stp,1) = min(min(cx(2:end-1,2:end-1)));
hist.cx(stp,2) = mean(mean(cx(2:end-1,2:end-1)));
hist.cx(stp,3) = max(max(cx(2:end-1,2:end-1)));

hist.cm(stp,1) = min(min(cm(2:end-1,2:end-1)));
hist.cm(stp,2) = mean(mean(cm(2:end-1,2:end-1)));
hist.cm(stp,3) = max(max(cm(2:end-1,2:end-1)));

hist.vf(stp,1) = min(min(vf(2:end-1,2:end-1)));
hist.vf(stp,2) = mean(mean(vf(2:end-1,2:end-1)));
hist.vf(stp,3) = max(max(vf(2:end-1,2:end-1)));

hist.vm(stp,1) = min(min(vm(2:end-1,2:end-1)));
hist.vm(stp,2) = mean(mean(vm(2:end-1,2:end-1)));
hist.vm(stp,3) = max(max(vm(2:end-1,2:end-1)));

hist.Gx(stp,1) = min(min(Gx(2:end-1,2:end-1)));
hist.Gx(stp,2) = mean(mean(Gx(2:end-1,2:end-1)));
hist.Gx(stp,3) = max(max(Gx(2:end-1,2:end-1)));

hist.Gf(stp,1) = min(min(Gf(2:end-1,2:end-1)));
hist.Gf(stp,2) = mean(mean(Gf(2:end-1,2:end-1)));
hist.Gf(stp,3) = max(max(Gf(2:end-1,2:end-1)));

hist.dV(stp,1) = min(min(VolSrc(2:end-1,2:end-1)));
hist.dV(stp,2) = mean(mean(VolSrc(2:end-1,2:end-1)));
hist.dV(stp,3) = max(max(VolSrc(2:end-1,2:end-1)));

hist.rho(stp,1) = min(min(rho(2:end-1,2:end-1)));
hist.rho(stp,2) = mean(mean(rho(2:end-1,2:end-1)));
hist.rho(stp,3) = max(max(rho(2:end-1,2:end-1)));

hist.eta(stp,1) = min(min(eta(2:end-1,2:end-1)));
hist.eta(stp,2) = geomean(geomean(eta(2:end-1,2:end-1)));
hist.eta(stp,3) = max(max(eta(2:end-1,2:end-1)));

hist.wx(stp,1) = min(min(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)));
hist.wx(stp,2) = mean(mean(abs((chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1))));
hist.wx(stp,3) = max(max(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1)));

hist.wf(stp,1) = min(min(-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1)));
hist.wf(stp,2) = mean(mean(abs((phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1))));
hist.wf(stp,3) = max(max(-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1)));

hist.it(stp,1) = min(min(it(2:end-1,2:end-1)));
hist.it(stp,2) = mean(mean(it(2:end-1,2:end-1)));
hist.it(stp,3) = max(max(it(2:end-1,2:end-1)));

hist.ct(stp,1) = min(min(ct(2:end-1,2:end-1)));
hist.ct(stp,2) = mean(mean(ct(2:end-1,2:end-1)));
hist.ct(stp,3) = max(max(ct(2:end-1,2:end-1)));

hist.si(stp,1) = min(min(si(2:end-1,2:end-1)));
hist.si(stp,2) = mean(mean(si(2:end-1,2:end-1)));
hist.si(stp,3) = max(max(si(2:end-1,2:end-1)));

hist.rip(stp,1) = min(min(rip(2:end-1,2:end-1)));
hist.rip(stp,2) = mean(mean(rip(2:end-1,2:end-1)));
hist.rip(stp,3) = max(max(rip(2:end-1,2:end-1)));

hist.rid(stp,1) = min(min(rid(2:end-1,2:end-1)));
hist.rid(stp,2) = mean(mean(rid(2:end-1,2:end-1)));
hist.rid(stp,3) = max(max(rid(2:end-1,2:end-1)));