% record run history

dsumMdtoo = dsumMdto; dsumMdto = dsumMdt;
dsumSdtoo = dsumSdto; dsumSdto = dsumSdt;
dsumCdtoo = dsumCdto; dsumCdto = dsumCdt;
dsumVdtoo = dsumVdto; dsumVdto = dsumVdt;

stp = max(1,step+1);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumM(stp) = sum(sum(rho(2:end-1,2:end-1)*h*h*1));  % [kg]
hist.sumS(stp) = sum(sum(  S(2:end-1,2:end-1)*h*h*1));  % [J]
hist.sumC(stp) = sum(sum(  C(2:end-1,2:end-1)*h*h*1));  % [kg]
hist.sumV(stp) = sum(sum(  V(2:end-1,2:end-1)*h*h*1));  % [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumMdt = sum(  X(2,2:end-1).*Wx(1,2:end-1)*h*1) - sum(X(end-1,2:end-1).*Wx(end,2:end-1)*h*1) ...
        + sum(  X(2:end-1,2).*Ux(2:end-1,1)*h*1) - sum(X(2:end-1,end-1).*Ux(2:end-1,end)*h*1) ...
        + sum(  F(2,2:end-1).*Wf(1,2:end-1)*h*1) - sum(F(end-1,2:end-1).*Wf(end,2:end-1)*h*1) ...
        + sum(  F(2:end-1,2).*Uf(2:end-1,1)*h*1) - sum(F(2:end-1,end-1).*Uf(2:end-1,end)*h*1) ...
        + sum(  M(2,2:end-1).*Wm(1,2:end-1)*h*1) - sum(M(end-1,2:end-1).*Wm(end,2:end-1)*h*1) ...
        + sum(  M(2:end-1,2).*Um(2:end-1,1)*h*1) - sum(M(2:end-1,end-1).*Um(2:end-1,end)*h*1);  % [kg/s]
dsumSdt = sum(sum(bnd_S*h*h*1)) + sum(sum(diss_h*h*h*1)) ...
        + sum(  X(2,2:end-1).*sx(2,2:end-1).*Wx(1,2:end-1)*h*1) - sum(X(end-1,2:end-1).*sx(end-1,2:end-1).*Wx(end,2:end-1)*h*1) ...
        + sum(  X(2:end-1,2).*sx(2:end-1,2).*Ux(2:end-1,1)*h*1) - sum(X(2:end-1,end-1).*sx(2:end-1,end-1).*Ux(2:end-1,end)*h*1) ...
        + sum(  F(2,2:end-1).*sf(2,2:end-1).*Wf(1,2:end-1)*h*1) - sum(F(end-1,2:end-1).*sf(end-1,2:end-1).*Wf(end,2:end-1)*h*1) ...
        + sum(  F(2:end-1,2).*sf(2:end-1,2).*Uf(2:end-1,1)*h*1) - sum(F(2:end-1,end-1).*sf(2:end-1,end-1).*Uf(2:end-1,end)*h*1) ...
        + sum(  M(2,2:end-1).*sm(2,2:end-1).*Wm(1,2:end-1)*h*1) - sum(M(end-1,2:end-1).*sm(end-1,2:end-1).*Wm(end,2:end-1)*h*1) ...
        + sum(  M(2:end-1,2).*sm(2:end-1,2).*Um(2:end-1,1)*h*1) - sum(M(2:end-1,end-1).*sm(2:end-1,end-1).*Um(2:end-1,end)*h*1);  % [J /s]
dsumCdt = sum(sum(bnd_C*h*h*1)) ...
        + sum(  X(2,2:end-1).*cx(2,2:end-1).*Wx(1,2:end-1)*h*1) - sum(X(end-1,2:end-1).*cx(end-1,2:end-1).*Wx(end,2:end-1)*h*1) ...
        + sum(  X(2:end-1,2).*cx(2:end-1,2).*Ux(2:end-1,1)*h*1) - sum(X(2:end-1,end-1).*cx(2:end-1,end-1).*Ux(2:end-1,end)*h*1) ...
        + sum(  M(2,2:end-1).*cm(2,2:end-1).*Wm(1,2:end-1)*h*1) - sum(M(end-1,2:end-1).*cm(end-1,2:end-1).*Wm(end,2:end-1)*h*1) ...
        + sum(  M(2:end-1,2).*cm(2:end-1,2).*Um(2:end-1,1)*h*1) - sum(M(2:end-1,end-1).*cm(2:end-1,end-1).*Um(2:end-1,end)*h*1);  % [kg/s]
dsumVdt = sum(sum(bnd_V*h*h*1)) ...
        + sum(  F(2,2:end-1).*vf(2,2:end-1).*Wf(1,2:end-1)*h*1) - sum(F(end-1,2:end-1).*vf(end-1,2:end-1).*Wf(end,2:end-1)*h*1) ...
        + sum(  F(2:end-1,2).*vf(2:end-1,2).*Uf(2:end-1,1)*h*1) - sum(F(2:end-1,end-1).*vf(2:end-1,end-1).*Uf(2:end-1,end)*h*1) ...
        + sum(  M(2,2:end-1).*vm(2,2:end-1).*Wm(1,2:end-1)*h*1) - sum(M(end-1,2:end-1).*vm(end-1,2:end-1).*Wm(end,2:end-1)*h*1) ...
        + sum(  M(2:end-1,2).*vm(2:end-1,2).*Um(2:end-1,1)*h*1) - sum(M(2:end-1,end-1).*vm(2:end-1,end-1).*Um(2:end-1,end)*h*1);  % [kg/s]

if step>=1; hist.DM(stp) = (alpha2*hist.DM(stp-1) + alpha3*hist.DM(max(1,stp-2)) + (beta1*dsumMdt + beta2*dsumMdto + beta3*dsumMdtoo)*dt)/alpha1; else; hist.DM(stp) = 0; end  % [kg]
if step>=1; hist.DS(stp) = (alpha2*hist.DS(stp-1) + alpha3*hist.DS(max(1,stp-2)) + (beta1*dsumSdt + beta2*dsumSdto + beta3*dsumSdtoo)*dt)/alpha1; else; hist.DS(stp) = 0; end  % [kg]
if step>=1; hist.DC(stp) = (alpha2*hist.DC(stp-1) + alpha3*hist.DC(max(1,stp-2)) + (beta1*dsumCdt + beta2*dsumCdto + beta3*dsumCdtoo)*dt)/alpha1; else; hist.DC(stp) = 0; end  % [kg]
if step>=1; hist.DV(stp) = (alpha2*hist.DV(stp-1) + alpha3*hist.DV(max(1,stp-2)) + (beta1*dsumVdt + beta2*dsumVdto + beta3*dsumVdtoo)*dt)/alpha1 .* any(v(:)>1e-6); else; hist.DV(stp) = 0; end  % [kg]

% record conservation error of mass M, heat H, major component C, volatile component V
hist.EM(stp) = (hist.sumM(stp) - hist.DM(stp))./hist.sumM(1) - 1;  % [kg/kg]
hist.ES(stp) = (hist.sumS(stp) - hist.DS(stp))./hist.sumS(1) - 1;  % [J /J ]
hist.EC(stp) = (hist.sumC(stp) - hist.DC(stp))./hist.sumC(1) - 1;  % [kg/kg]
hist.EV(stp) =((hist.sumV(stp) - hist.DV(stp))./(hist.sumV(1)+TINY) - 1);  % [kg/kg]

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

for i=1:cal.nc
    hist.c_cmp(stp,1,i) = min(min(c_cmp(2:end-1,2:end-1,i)));
    hist.c_cmp(stp,2,i) = mean(mean(c_cmp(2:end-1,2:end-1,i)));
    hist.c_cmp(stp,3,i) = max(max(c_cmp(2:end-1,2:end-1,i)));

    hist.c_oxd(stp,1,i) = min(min(c_oxd(2:end-1,2:end-1,i)));
    hist.c_oxd(stp,2,i) = mean(mean(c_oxd(2:end-1,2:end-1,i)));
    hist.c_oxd(stp,3,i) = max(max(c_oxd(2:end-1,2:end-1,i)));
end

if any(v(:)>1e-6)
    hist.v(stp,1) = min(min(v(2:end-1,2:end-1)));
    hist.v(stp,2) = mean(mean(v(2:end-1,2:end-1)));
    hist.v(stp,3) = max(max(v(2:end-1,2:end-1)));
else
    hist.v(stp,1:3) = NaN;
end

indx = repmat(x>1e-6,1,1,cal.nc);
if any(indx(:)>0)
    hist.cx(stp,1) = min(min(cx(indx(2:end-1,2:end-1,1))));
    hist.cx(stp,2) = sum(sum(cx(2:end-1,2:end-1).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
    hist.cx(stp,3) = max(max(cx(indx(2:end-1,2:end-1,1))));

    for i=1:cal.nc
        hist.cx_cmp(stp,1,i) = min(min(cx_cmp(indx(2:end-1,2:end-1,i))));
        hist.cx_cmp(stp,2,i) = sum(sum(cx_cmp(2:end-1,2:end-1,i).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cx_cmp(stp,3,i) = max(max(cx_cmp(indx(2:end-1,2:end-1,i))));

        hist.cx_oxd(stp,1,i) = min(min(cx_oxd(indx(2:end-1,2:end-1,i))));
        hist.cx_oxd(stp,2,i) = sum(sum(cx_oxd(2:end-1,2:end-1,i).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cx_oxd(stp,3,i) = max(max(cx_oxd(indx(2:end-1,2:end-1,i))));
    end

    hist.rhox(stp,1) = min(min(rhox(indx(2:end-1,2:end-1,1))));
    hist.rhox(stp,2) = sum(sum(rhox(2:end-1,2:end-1).*x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(x(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
    hist.rhox(stp,3) = max(max(rhox(indx(2:end-1,2:end-1,1))));
else
    hist.cx(stp,1:3) = NaN;
    hist.cx_cmp(stp,1:3,1:cal.nc) = NaN;
    hist.cx_oxd(stp,1:3,1:cal.nc) = NaN;
    hist.rhox(stp,1:3) = NaN;
end

indm = repmat(m>1e-6,1,1,cal.nc);
if any(indm(:)>0)
    hist.cm(stp,1) = min(min(cm(indm(2:end-1,2:end-1,1))));
    hist.cm(stp,2) = sum(sum(cm(2:end-1,2:end-1).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
    hist.cm(stp,3) = max(max(cm(indm(2:end-1,2:end-1,1))));
    
    for i=1:cal.nc
        hist.cm_cmp(stp,1,i) = min(min(cm_cmp(indm(2:end-1,2:end-1,i))));
        hist.cm_cmp(stp,2,i) = sum(sum(cm_cmp(2:end-1,2:end-1,i).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cm_cmp(stp,3,i) = max(max(cm_cmp(indm(2:end-1,2:end-1,i))));

        hist.cm_oxd(stp,1,i) = min(min(cm_oxd(indm(2:end-1,2:end-1,i))));
        hist.cm_oxd(stp,2,i) = sum(sum(cm_oxd(2:end-1,2:end-1,i).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.cm_oxd(stp,3,i) = max(max(cm_oxd(indm(2:end-1,2:end-1,i))));
    end

    if any(v(:)>1e-6)
        hist.vm(stp,1) = min(min(vm(indm(2:end-1,2:end-1,1))));
        hist.vm(stp,2) = sum(sum(vm(2:end-1,2:end-1).*m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
        hist.vm(stp,3) = max(max(vm(indm(2:end-1,2:end-1,1))));
    else
        hist.vm(stp,1:3) = NaN;
    end
    
    hist.rhom(stp,1) = min(min(rhom(indm(2:end-1,2:end-1,1))));
    hist.rhom(stp,2) = sum(sum(rhom(2:end-1,2:end-1).*m(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1)));
    hist.rhom(stp,3) = max(max(rhom(indm(2:end-1,2:end-1,1))));

    hist.etam(stp,1) = min(min(etam(indm(2:end-1,2:end-1,1))));
    hist.etam(stp,2) = sum(sum(etam(2:end-1,2:end-1).*m(2:end-1,2:end-1)))./sum(sum(m(2:end-1,2:end-1)));
    hist.etam(stp,3) = max(max(etam(indm(2:end-1,2:end-1,1))));
else
    hist.cm(stp,1:3) = NaN;
    hist.cm_cmp(stp,1:3,1:cal.nc) = NaN;
    hist.cm_oxd(stp,1:3,1:cal.nc) = NaN;
    hist.vm(stp,1:3) = NaN;
    hist.rhom(stp,1:3) = NaN;
    hist.etam(stp,1:3) = NaN;
end

indf = f>1e-6;
if any(indf(:)>0)
    hist.vf(stp,1) = min(min(vf(indf(2:end-1,2:end-1))));
    hist.vf(stp,2) = sum(sum(vf(2:end-1,2:end-1).*f(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(f(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
    hist.vf(stp,3) = max(max(vf(indf(2:end-1,2:end-1))));
    
    hist.rhof(stp,1) = min(min(rhof(indf(2:end-1,2:end-1))));
    hist.rhof(stp,2) = sum(sum(rhof(2:end-1,2:end-1).*f(2:end-1,2:end-1).*rho(2:end-1,2:end-1)))./sum(sum(f(2:end-1,2:end-1).*rho(2:end-1,2:end-1)));
    hist.rhof(stp,3) = max(max(rhof(indf(2:end-1,2:end-1))));
else
    hist.vf(stp,1:3) = NaN;
    hist.rhof(stp,1:3) = NaN;
end

hist.Gx(stp,1) = min(min(Gx));
hist.Gx(stp,2) = mean(mean(Gx));
hist.Gx(stp,3) = max(max(Gx));

hist.Gf(stp,1) = min(min(Gf));
hist.Gf(stp,2) = mean(mean(Gf));
hist.Gf(stp,3) = max(max(Gf));

hist.dV(stp,1) = min(min(Div_V));
hist.dV(stp,2) = mean(mean(Div_V));
hist.dV(stp,3) = max(max(Div_V));

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

for k = 1:length(te0)
    hist.te(stp,k,1) = min(min(te(2:end-1,2:end-1,k)));
    hist.te(stp,k,2) = mean(mean(te(2:end-1,2:end-1,k)));
    hist.te(stp,k,3) = max(max(te(2:end-1,2:end-1,k)));
end
for k = 1:length(ir0)
    hist.ir(stp,k,1) = min(min(ir(2:end-1,2:end-1,k)));
    hist.ir(stp,k,2) = mean(mean(ir(2:end-1,2:end-1,k)));
    hist.ir(stp,k,3) = max(max(ir(2:end-1,2:end-1,k)));
end

% fraction, composition, and temperature of eruptible magma suspension (mu>0.55)
indmagma = max(0,min(1,(1+erf((mu-0.55)./0.05))/2));
hist.Fmagma(stp) = sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Cmagma(stp) = sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*c(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*h^2));
hist.Tmagma(stp) = sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmagma(2:end-1,2:end-1).*h^2));

% fraction, composition, and temperature of plutonic rock (mu<0.15)
indpluton = max(0,min(1,(1+erf((chi-0.85)./0.05))/2));
hist.Fpluton(stp) = sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Cpluton(stp) = sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*c(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*h^2));
hist.Tpluton(stp) = sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indpluton(2:end-1,2:end-1).*h^2));

% fraction, composition, and temperature of magma mush (0.15<mu<0.55)
indmush = max(0,min(1,1-indmagma-indpluton));
hist.Fmush(stp) = sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Cmush(stp) = sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*c(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*h^2));
hist.Tmush(stp) = sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmush(2:end-1,2:end-1).*h^2));

% fraction, crystallinity, and temperature of felsic materials (c > (perCm_cphs1)/2)
indfelsic = max(0,min(1,(1+erf((c-(cal.perCm+cal.cphs1)/2)./0.005))/2));
hist.Ffelsic(stp) = sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Xfelsic(stp) = sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*h^2));
hist.Tfelsic(stp) = sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indfelsic(2:end-1,2:end-1).*h^2));

% fraction, crystallinity, and temperature of intermediate materials (perCm < c < (perCm_cphs1)/2)
indinterm = max(0,min(1,(1+erf((c-cal.perCm)./0.005))/2 .* (1-indfelsic)));
hist.Finterm(stp) = sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Xinterm(stp) = sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*h^2));
hist.Tinterm(stp) = sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indinterm(2:end-1,2:end-1).*h^2));

% fraction, crystallinity, and temperature of mafic materials (perCx < c < perCm)
indmafic = max(0,min(1,(1+erf((c-cal.perCx)./0.005))/2 .* (1-indinterm-indfelsic)));
hist.Fmafic(stp) = sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Xmafic(stp) = sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*h^2));
hist.Tmafic(stp) = sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indmafic(2:end-1,2:end-1).*h^2));

% fraction, crystallinity, and temperature of ultramafic materials (c < perCx)
indultram = max(0,min(1,1-indmafic-indinterm-indfelsic));
hist.Fultram(stp) = sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2));
hist.Xultram(stp) = sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*x(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*h^2));
hist.Tultram(stp) = sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*T(2:end-1,2:end-1).*h^2))./sum(sum(rho(2:end-1,2:end-1).*indultram(2:end-1,2:end-1).*h^2));

% differentiation index
nobnd = false(size(c));
nobnd(inz,inx) = bndshape<1e-2;
if any(nobnd(:))
    hist.Rdiff(stp) = (max(max(c(nobnd)))-min(min(c(nobnd))))./(cal.cphs1-cal.cphs0);
end

% index of assimilation
if step>1
    hist.Ra (stp) = hist.Ra (stp-1) + sum(sum(bndshape.*rho(2:end-1,2:end-1)./tau_a.*h^2))./sum(sum(rho(2:end-1,2:end-1).*h^2)).*dt;
else
    hist.Ra(stp,1) = 0;
end