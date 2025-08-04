% record run history

dsumSdtoo = dsumSdto; dsumSdto = dsumSdt;
dsumBdtoo = dsumBdto; dsumBdto = dsumBdt;
dsumMdtoo = dsumMdto; dsumMdto = dsumMdt;
dsumXdtoo = dsumXdto; dsumXdto = dsumXdt;
dsumFdtoo = dsumFdto; dsumFdto = dsumFdt;
dsumCdtoo = dsumCdto; dsumCdto = dsumCdt;
dsumTdtoo = dsumTdto; dsumTdto = dsumTdt;

stp = max(1,step/nrh);

% record model time
hist.time(stp) = time;

% record total mass, heat, component mass in model (assume hy = 1, unit length in third dimension)
hist.sumS(stp  ) = sum(  S(:)*h*h*1)+eps;  % [J/K]
hist.sumB(stp  ) = sum(rho(:)*h*h*1)+eps;  % [kg]
hist.sumM(stp  ) = sum(  M(:)*h*h*1)+eps;  % [kg]
hist.sumX(stp  ) = sum(  X(:)*h*h*1)+eps;  % [kg]
hist.sumF(stp  ) = sum(  F(:)*h*h*1)+eps;  % [kg]
hist.sumC(stp,:) = squeeze(sum(sum(C  *h*h*1,1),2))+eps; % [kg]
hist.sumT(stp,:) = squeeze(sum(sum(TRC*h*h*1,1),2))+eps; % [kg]

% record expected rates of change by volume change and imposed boundaries layers
dsumSdt =  sum(sum(bnd_S*h*h*1)) + sum(sum(diss_h*h*h*1)) ...
        +  sum(X(1,:).*sx(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*sx(end,:).*Wx(end,2:end-1)*h*1) ...
        + (sum(X(:,1).*sx(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*sx(:,end).*Ux(2:end-1,end)*h*1))*(1-periodic) ...
        +  sum(F(1,:).*sf(1,:).*Wf(1,2:end-1)*h*1) - sum(F(end,:).*sf(end,:).*Wf(end,2:end-1)*h*1) ...
        + (sum(F(:,1).*sf(:,1).*Uf(2:end-1,1)*h*1) - sum(F(:,end).*sf(:,end).*Uf(2:end-1,end)*h*1))*(1-periodic) ...
        +  sum(M(1,:).*sm(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*sm(end,:).*Wm(end,2:end-1)*h*1) ...
        + (sum(M(:,1).*sm(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*sm(:,end).*Um(2:end-1,end)*h*1))*(1-periodic);  % [J/K/s]
dsumBdt =  sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1) ...
        + (sum(X(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*Ux(2:end-1,end)*h*1))*(1-periodic) ...
        +  sum(F(1,:).*Wf(1,2:end-1)*h*1) - sum(F(end,:).*Wf(end,2:end-1)*h*1) ...
        + (sum(F(:,1).*Uf(2:end-1,1)*h*1) - sum(F(:,end).*Uf(2:end-1,end)*h*1))*(1-periodic) ...
        +  sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1) ...
        + (sum(M(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*Um(2:end-1,end)*h*1))*(1-periodic);  % [kg/s]
dsumMdt =  sum(sum(Gm*h*h*1)) ...
        +  sum(M(1,:).*Wm(1,2:end-1)*h*1) - sum(M(end,:).*Wm(end,2:end-1)*h*1) ...
        + (sum(M(:,1).*Um(2:end-1,1)*h*1) - sum(M(:,end).*Um(2:end-1,end)*h*1))*(1-periodic);  % [kg/s]
dsumXdt =  sum(sum(Gx*h*h*1)) ...
        +  sum(X(1,:).*Wx(1,2:end-1)*h*1) - sum(X(end,:).*Wx(end,2:end-1)*h*1) ...
        + (sum(X(:,1).*Ux(2:end-1,1)*h*1) - sum(X(:,end).*Ux(2:end-1,end)*h*1))*(1-periodic);  % [kg/s]
dsumFdt =  sum(sum(Gf*h*h*1)) ...
        +  sum(F(1,:).*Wf(1,2:end-1)*h*1) - sum(F(end,:).*Wf(end,2:end-1)*h*1) ...
        + (sum(F(:,1).*Uf(2:end-1,1)*h*1) - sum(F(:,end).*Uf(2:end-1,end)*h*1))*(1-periodic);  % [kg/s]
dsumCdt =  squeeze(sum(sum(bnd_C*h*h*1,1),2) ...
        +  sum(X(1,:).*cx(1,:,:).*Wx(1,2:end-1)*h*1,2) - sum(X(end,:).*cx(end,:,:).*Wx(end,2:end-1)*h*1,2) ...
        + (sum(X(:,1).*cx(:,1,:).*Ux(2:end-1,1)*h*1,1) - sum(X(:,end).*cx(:,end,:).*Ux(2:end-1,end)*h*1,1))*(1-periodic) ...
        +  sum(F(1,:).*cf(1,:,:).*Wf(1,2:end-1)*h*1,2) - sum(F(end,:).*cf(end,:,:).*Wf(end,2:end-1)*h*1,2) ...
        + (sum(F(:,1).*cf(:,1,:).*Uf(2:end-1,1)*h*1,1) - sum(F(:,end).*cf(:,end,:).*Uf(2:end-1,end)*h*1,1))*(1-periodic) ...
        +  sum(M(1,:).*cm(1,:,:).*Wm(1,2:end-1)*h*1,2) - sum(M(end,:).*cm(end,:,:).*Wm(end,2:end-1)*h*1,2) ...
        + (sum(M(:,1).*cm(:,1,:).*Um(2:end-1,1)*h*1,1) - sum(M(:,end).*cm(:,end,:).*Um(2:end-1,end)*h*1,1))*(1-periodic)).';  % [kg/s]
dsumTdt =  squeeze(sum(sum(bnd_TRC*h*h*1,1),2) ...
        +  sum(X(1,:).*trcx(1,:,:).*Wx(1,2:end-1)*h*1,2) - sum(X(end,:).*trcx(end,:,:).*Wx(end,2:end-1)*h*1,2) ...
        + (sum(X(:,1).*trcx(:,1,:).*Ux(2:end-1,1)*h*1,1) - sum(X(:,end).*trcx(:,end,:).*Ux(2:end-1,end)*h*1,1))*(1-periodic) ...
        +  sum(M(1,:).*trcm(1,:,:).*Wm(1,2:end-1)*h*1,2) - sum(M(end,:).*trcm(end,:,:).*Wm(end,2:end-1)*h*1,2) ...
        + (sum(M(:,1).*trcm(:,1,:).*Um(2:end-1,1)*h*1,1) - sum(M(:,end).*trcm(:,end,:).*Um(2:end-1,end)*h*1,1))*(1-periodic)).';  % [kg/s]

if step>=2; hist.DS(stp  ) = (a2*hist.DS(max(1,stp-1)  ) + a3*hist.DS(max(1,stp-2)  ) + (b1*dsumSdt + b2*dsumSdto + b3*dsumSdtoo)*dt)/a1; else; hist.DS(stp  ) = 0; end  % [kg]
if step>=2; hist.DB(stp  ) = (a2*hist.DB(max(1,stp-1)  ) + a3*hist.DB(max(1,stp-2)  ) + (b1*dsumBdt + b2*dsumBdto + b3*dsumBdtoo)*dt)/a1; else; hist.DB(stp  ) = 0; end  % [kg]
if step>=2; hist.DM(stp  ) = (a2*hist.DM(max(1,stp-1)  ) + a3*hist.DM(max(1,stp-2)  ) + (b1*dsumMdt + b2*dsumMdto + b3*dsumMdtoo)*dt)/a1; else; hist.DM(stp  ) = 0; end  % [kg]
if step>=2; hist.DX(stp  ) = (a2*hist.DX(max(1,stp-1)  ) + a3*hist.DX(max(1,stp-2)  ) + (b1*dsumXdt + b2*dsumXdto + b3*dsumXdtoo)*dt)/a1; else; hist.DX(stp  ) = 0; end  % [kg]
if step>=2; hist.DF(stp  ) = (a2*hist.DF(max(1,stp-1)  ) + a3*hist.DF(max(1,stp-2)  ) + (b1*dsumFdt + b2*dsumFdto + b3*dsumFdtoo)*dt)/a1; else; hist.DF(stp  ) = 0; end  % [kg]
if step>=2; hist.DC(stp,:) = (a2*hist.DC(max(1,stp-1),:) + a3*hist.DC(max(1,stp-2),:) + (b1*dsumCdt + b2*dsumCdto + b3*dsumCdtoo)*dt)/a1; else; hist.DC(stp,:) = zeros(1,cal.ncmp); end  % [kg]
if step>=2; hist.DT(stp,:) = (a2*hist.DT(max(1,stp-1),:) + a3*hist.DT(max(1,stp-2),:) + (b1*dsumTdt + b2*dsumTdto + b3*dsumTdtoo)*dt)/a1; else; hist.DT(stp,:) = zeros(1,cal.ntrc); end  % [kg]

% record conservation error of mass M, heat S, components C
hist.ES(stp  ) = (hist.sumS(stp  ) - hist.DS(stp  ))./hist.sumS(1  ) - 1;  % [JK/JK]
hist.EB(stp  ) = (hist.sumB(stp  ) - hist.DB(stp  ))./hist.sumB(1  ) - 1;  % [kg/kg]
hist.EM(stp  ) = (hist.sumM(stp  ) - hist.DM(stp  ))./hist.sumM(1  ) - 1;  % [kg/kg]
hist.EX(stp  ) = (hist.sumX(stp  ) - hist.DX(stp  ))./hist.sumX(1  ) - 1;  % [kg/kg]
hist.EF(stp  ) = (hist.sumF(stp  ) - hist.DF(stp  ))./hist.sumF(1  ) - 1;  % [kg/kg]
hist.EC(stp,:) = (hist.sumC(stp,:) - hist.DC(stp,:))./hist.sumC(1,:) - 1;  % [kg/kg]
hist.ET(stp,:) = (hist.sumT(stp,:) - hist.DT(stp,:))./hist.sumT(1,:) - 1;  % [kg/kg]

% % alternative way to calculate i/o rate, tested but accumulates more error
% if stp>=2; hist.DS2(stp  ) = hist.DS2(max(1,stp-1)  ) + dsumSdt*dt; else; hist.DS2(stp  ) = 0; end  % [J/K]
% if stp>=2; hist.DB2(stp  ) = hist.DB2(max(1,stp-1)  ) + dsumBdt*dt; else; hist.DB2(stp  ) = 0; end  % [kg]
% if stp>=2; hist.DM2(stp  ) = hist.DM2(max(1,stp-1)  ) + dsumMdt*dt; else; hist.DM2(stp  ) = 0; end  % [kg]
% if stp>=2; hist.DX2(stp  ) = hist.DX2(max(1,stp-1)  ) + dsumXdt*dt; else; hist.DX2(stp  ) = 0; end  % [kg]
% if stp>=2; hist.DF2(stp  ) = hist.DF2(max(1,stp-1)  ) + dsumFdt*dt; else; hist.DF2(stp  ) = 0; end  % [kg]
% if stp>=2; hist.DC2(stp,:) = hist.DC2(max(1,stp-1),:) + dsumCdt*dt; else; hist.DC2(stp,:) = zeros(1,cal.ncmp); end  % [kg]
% if stp>=2; hist.DT2(stp,:) = hist.DT2(max(1,stp-1),:) + dsumTdt*dt; else; hist.DT2(stp,:) = zeros(1,cal.ntrc); end  % [kg]
% 
% hist.ES2(stp  ) = (hist.sumS(stp  ) - hist.DS2(stp  ))./hist.sumS(1  ) - 1;  % [JK/JK]
% hist.EB2(stp  ) = (hist.sumB(stp  ) - hist.DB2(stp  ))./hist.sumB(1  ) - 1;  % [kg/kg]
% hist.EM2(stp  ) = (hist.sumM(stp  ) - hist.DM2(stp  ))./hist.sumM(1  ) - 1;  % [kg/kg]
% hist.EX2(stp  ) = (hist.sumX(stp  ) - hist.DX2(stp  ))./hist.sumX(1  ) - 1;  % [kg/kg]
% hist.EF2(stp  ) = (hist.sumF(stp  ) - hist.DF2(stp  ))./hist.sumF(1  ) - 1;  % [kg/kg]
% hist.EC2(stp,:) = (hist.sumC(stp,:) - hist.DC2(stp,:))./hist.sumC(1,:) - 1;  % [kg/kg]
% hist.ET2(stp,:) = (hist.sumT(stp,:) - hist.DT2(stp,:))./hist.sumT(1,:) - 1;  % [kg/kg]

% record variable and coefficient diagnostics
hist.W(stp,1) = min(min(-W(:,2:end-1)));
hist.W(stp,2) = mean(mean(abs(W(:,2:end-1))));
hist.W(stp,3) = max(max(-W(:,2:end-1)));

hist.U(stp,1) = min(min(U(2:end-1,:)));
hist.U(stp,2) = mean(mean(abs(U(2:end-1,:))));
hist.U(stp,3) = max(max(U(2:end-1,:)));

hist.P(stp,1) = min(min(P(2:end-1,2:end-1)));
hist.P(stp,2) = mean(mean(P(2:end-1,2:end-1)));
hist.P(stp,3) = max(max(P(2:end-1,2:end-1)));

hist.Pt(stp,1) = min(min(Pt));
hist.Pt(stp,2) = mean(mean(abs(Pt)));
hist.Pt(stp,3) = max(max(Pt));

hist.Pchmb(stp,1) = Pchmb;

hist.x(stp,1) = min(min(x));
hist.x(stp,2) = mean(mean(x));
hist.x(stp,3) = max(max(x));

hist.f(stp,1) = min(min(f));
hist.f(stp,2) = mean(mean(f));
hist.f(stp,3) = max(max(f));

hist.m(stp,1) = min(min(m));
hist.m(stp,2) = mean(mean(m));
hist.m(stp,3) = max(max(m));

hist.chi(stp,1) = min(min(chi));
hist.chi(stp,2) = mean(mean(chi));
hist.chi(stp,3) = max(max(chi));

hist.phi(stp,1) = min(min(phi));
hist.phi(stp,2) = mean(mean(phi));
hist.phi(stp,3) = max(max(phi));

hist.mu(stp,1) = min(min(mu));
hist.mu(stp,2) = mean(mean(mu));
hist.mu(stp,3) = max(max(mu));

hist.T(stp,1) = min(min(T));
hist.T(stp,2) = mean(mean(T));
hist.T(stp,3) = max(max(T));

hist.Tsol(stp,1) = min(min(cal.Tsol));
hist.Tsol(stp,2) = mean(mean(cal.Tsol));
hist.Tsol(stp,3) = max(max(cal.Tsol));

hist.Tliq(stp,1) = min(min(cal.Tliq));
hist.Tliq(stp,2) = mean(mean(cal.Tliq));
hist.Tliq(stp,3) = max(max(cal.Tliq));

if Nz==1 && dPdT
    hist.Tm(stp,:) = cal.Tm;
end

for i=1:cal.ncmp
    hist.c(stp,1,i) = min(min(c(:,:,i)));
    hist.c(stp,2,i) = mean(mean(c(:,:,i)));
    hist.c(stp,3,i) = max(max(c(:,:,i)));
end
for i=1:cal.noxd
    hist.c_oxd(stp,1,i) = min(min(c_oxd(:,:,i)));
    hist.c_oxd(stp,2,i) = mean(mean(c_oxd(:,:,i)));
    hist.c_oxd(stp,3,i) = max(max(c_oxd(:,:,i)));
end
for i=1:cal.nmem
    hist.c_mem(stp,1,i) = min(min(c_mem(:,:,i)));
    hist.c_mem(stp,2,i) = mean(mean(c_mem(:,:,i)));
    hist.c_mem(stp,3,i) = max(max(c_mem(:,:,i)));
end
for i=1:cal.nmsy
    hist.c_msy(stp,1,i) = min(min(c_msy(:,:,i)));
    hist.c_msy(stp,2,i) = mean(mean(c_msy(:,:,i)));
    hist.c_msy(stp,3,i) = max(max(c_msy(:,:,i)));
end

for i=1:cal.ncmp
    hist.cx(stp,1,i) = min(min(cx(:,:,i)));
    hist.cx(stp,2,i) = sum(sum(cx(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx(stp,3,i) = max(max(cx(:,:,i)));
end
for i=1:cal.noxd
    hist.cx_oxd(stp,1,i) = min(min(cx_oxd(:,:,i)));
    hist.cx_oxd(stp,2,i) = sum(sum(cx_oxd(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx_oxd(stp,3,i) = max(max(cx_oxd(:,:,i)));
end
for i=1:cal.nmem
    hist.cx_mem(stp,1,i) = min(min(cx_mem(:,:,i)));
    hist.cx_mem(stp,2,i) = sum(sum(cx_mem(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx_mem(stp,3,i) = max(max(cx_mem(:,:,i)));
end
for i=1:cal.nmsy
    hist.cx_msy(stp,1,i) = min(min(cx_msy(:,:,i)));
    hist.cx_msy(stp,2,i) = sum(sum(cx_msy(:,:,i).*x.*rho))./sum(sum(x.*rho));
    hist.cx_msy(stp,3,i) = max(max(cx_msy(:,:,i)));
end

hist.rhox(stp,1) = min(min(rhox));
hist.rhox(stp,2) = sum(sum(rhox.*x.*rho))./sum(sum(x.*rho));
hist.rhox(stp,3) = max(max(rhox));

for i=1:cal.ncmp
    hist.cm(stp,1,i) = min(min(cm(:,:,i)));
    hist.cm(stp,2,i) = sum(sum(cm(:,:,i).*m.*rho))./sum(sum(m.*rho));
    hist.cm(stp,3,i) = max(max(cm(:,:,i)));
end
for i=1:cal.noxd
    hist.cm_oxd(stp,1,i) = min(min(cm_oxd(:,:,i)));
    hist.cm_oxd(stp,2,i) = sum(sum(cm_oxd(:,:,i).*m.*rho))./sum(sum(m.*rho));
    hist.cm_oxd(stp,3,i) = max(max(cm_oxd(:,:,i)));
end
for i=1:cal.nmem
    hist.cm_mem(stp,1,i) = min(min(cm_mem(:,:,i)));
    hist.cm_mem(stp,2,i) = sum(sum(cm_mem(:,:,i).*m.*rho))./sum(sum(m.*rho));
    hist.cm_mem(stp,3,i) = max(max(cm_mem(:,:,i)));
end
for i=1:cal.nmsy
    hist.cm_msy(stp,1,i) = min(min(cm_msy(:,:,i)));
    hist.cm_msy(stp,2,i) = sum(sum(cm_msy(:,:,i).*m.*rho))./sum(sum(m.*rho));
    hist.cm_msy(stp,3,i) = max(max(cm_msy(:,:,i)));
end

hist.rhom(stp,1) = min(min(rhom));
hist.rhom(stp,2) = sum(sum(rhom.*m))./sum(sum(m));
hist.rhom(stp,3) = max(max(rhom));

hist.etam(stp,1) = min(min(etam));
hist.etam(stp,2) = sum(sum(etam.*m))./sum(sum(m));
hist.etam(stp,3) = max(max(etam));

hist.Gm(stp,1) = min(min(Gm));
hist.Gm(stp,2) = mean(mean(Gm));
hist.Gm(stp,3) = max(max(Gm));

hist.Gx(stp,1) = min(min(Gx));
hist.Gx(stp,2) = mean(mean(Gx));
hist.Gx(stp,3) = max(max(Gx));

hist.Gf(stp,1) = min(min(Gf));
hist.Gf(stp,2) = mean(mean(Gf));
hist.Gf(stp,3) = max(max(Gf));

hist.dV(stp,1) = min(min(Div_V));
hist.dV(stp,2) = mean(mean(Div_V));
hist.dV(stp,3) = max(max(Div_V));

hist.rho(stp,1) = min(min(rho));
hist.rho(stp,2) = mean(mean(rho));
hist.rho(stp,3) = max(max(rho));

hist.eta(stp,1) = min(min(eta));
hist.eta(stp,2) = geomean(geomean(eta));
hist.eta(stp,3) = max(max(eta));

hist.wx(stp,1) = min(min(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)));
hist.wx(stp,2) = mean(mean(abs((chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1))));
hist.wx(stp,3) = max(max(-(chi([1,1:end],:)+chi([1:end,end],:))/2.*wx(:,2:end-1)));

hist.wf(stp,1) = min(min(-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1)));
hist.wf(stp,2) = mean(mean(abs((phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1))));
hist.wf(stp,3) = max(max(-(phi([1,1:end],:)+phi([1:end,end],:))/2.*wf(:,2:end-1)));

hist.dH(stp,1) = -sum(bnd_S(1:round(Nz/2)  ,:).*T(1:round(Nz/2)  ,:).*h^2,'all')/L;
hist.dH(stp,2) = -sum(bnd_S(round(Nz/2):end,:).*T(round(Nz/2):end,:).*h^2,'all')/L;
hist.dH(stp,3) = -sum(bnd_S.*T.*h^2,'all')/L;


% % fraction, composition, and temperature of eruptible magma suspension (mu>0.55)
% indmagma = max(0,min(1,(1+erf((mu-0.55)./0.05))/2));
% hist.Fmagma(stp) = sum(sum(rho.*indmagma.*h^2))./sum(sum(rho.*h^2));
% hist.Cmagma(stp) = sum(sum(rho.*indmagma.*c.*h^2))./sum(sum(rho.*indmagma.*h^2));
% hist.Tmagma(stp) = sum(sum(rho.*indmagma.*T.*h^2))./sum(sum(rho.*indmagma.*h^2));
% 
% % fraction, composition, and temperature of plutonic rock (mu<0.15)
% indpluton = max(0,min(1,(1+erf((chi-0.85)./0.05))/2));
% hist.Fpluton(stp) = sum(sum(rho.*indpluton.*h^2))./sum(sum(rho.*h^2));
% hist.Cpluton(stp) = sum(sum(rho.*indpluton.*c.*h^2))./sum(sum(rho.*indpluton.*h^2));
% hist.Tpluton(stp) = sum(sum(rho.*indpluton.*T.*h^2))./sum(sum(rho.*indpluton.*h^2));
% 
% % fraction, composition, and temperature of magma mush (0.15<mu<0.55)
% indmush = max(0,min(1,1-indmagma-indpluton));
% hist.Fmush(stp) = sum(sum(rho.*indmush.*h^2))./sum(sum(rho.*h^2));
% hist.Cmush(stp) = sum(sum(rho.*indmush.*c.*h^2))./sum(sum(rho.*indmush.*h^2));
% hist.Tmush(stp) = sum(sum(rho.*indmush.*T.*h^2))./sum(sum(rho.*indmush.*h^2));

% % fraction, crystallinity, and temperature of felsic materials (c > (perCm_cphs1)/2)
% indfelsic = max(0,min(1,(1+erf((c(:,:,1)-(cal.perCm+cal.cphs1)/2)./0.005))/2));
% hist.Ffelsic(stp) = sum(sum(rho.*indfelsic.*h^2))./sum(sum(rho.*h^2));
% hist.Xfelsic(stp) = sum(sum(rho.*indfelsic.*x.*h^2))./sum(sum(rho.*indfelsic.*h^2));
% hist.Tfelsic(stp) = sum(sum(rho.*indfelsic.*T.*h^2))./sum(sum(rho.*indfelsic.*h^2));
% 
% % fraction, crystallinity, and temperature of intermediate materials (perCm < c < (perCm_cphs1)/2)
% indinterm = max(0,min(1,(1+erf((c-cal.perCm)./0.005))/2 .* (1-indfelsic)));
% hist.Finterm(stp) = sum(sum(rho.*indinterm.*h^2))./sum(sum(rho.*h^2));
% hist.Xinterm(stp) = sum(sum(rho.*indinterm.*x.*h^2))./sum(sum(rho.*indinterm.*h^2));
% hist.Tinterm(stp) = sum(sum(rho.*indinterm.*T.*h^2))./sum(sum(rho.*indinterm.*h^2));
% 
% % fraction, crystallinity, and temperature of mafic materials (perCx < c < perCm)
% indmafic = max(0,min(1,(1+erf((c-cal.perCx)./0.005))/2 .* (1-indinterm-indfelsic)));
% hist.Fmafic(stp) = sum(sum(rho.*indmafic.*h^2))./sum(sum(rho.*h^2));
% hist.Xmafic(stp) = sum(sum(rho.*indmafic.*x.*h^2))./sum(sum(rho.*indmafic.*h^2));
% hist.Tmafic(stp) = sum(sum(rho.*indmafic.*T.*h^2))./sum(sum(rho.*indmafic.*h^2));
% 
% % fraction, crystallinity, and temperature of ultramafic materials (c < perCx)
% indultram = max(0,min(1,1-indmafic-indinterm-indfelsic));
% hist.Fultram(stp) = sum(sum(rho.*indultram.*h^2))./sum(sum(rho.*h^2));
% hist.Xultram(stp) = sum(sum(rho.*indultram.*x.*h^2))./sum(sum(rho.*indultram.*h^2));
% hist.Tultram(stp) = sum(sum(rho.*indultram.*T.*h^2))./sum(sum(rho.*indultram.*h^2));

% % differentiation index
% nobnd = topshape<1e-3 & botshape<1e-3 & sdsshape<1e-3;
% if any(nobnd(:))
%     hist.Rdiff(stp) = (max(max(c(nobnd)))-min(min(c(nobnd))))./(cal.cphs1-cal.cphs0);
% end

% % index of assimilation
% if step>1
%     hist.Ra (stp) = hist.Ra (stp-1) + sum(sum((topshape+botshape+sdsshape).*rho./tau_a.*h^2))./sum(sum(rho.*h^2)).*dt;
% else
%     hist.Ra(stp,1) = 0;
% end