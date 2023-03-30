% calculate phase equilibrium

function [xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T0,c,v,P,SiO2m,cal,TINY)

Tphs0d = cal.Tphs0;
Tphs1d = cal.Tphs1;
cphs0  = cal.cphs0;
cphs1  = cal.cphs1;
perTdx = cal.perTx;
perTdm = cal.perTm;
perCx  = cal.perCx;
perCm  = cal.perCm;
clap   = cal.clap;
dTH2O  = cal.dTH2O;
PhDg   = cal.PhDg;

perCm = (perCm-cphs0)./(cphs1-cphs0);
perCx = (perCx-cphs0)./(cphs1-cphs0);
SiO2m = (SiO2m-cphs0)./(cphs1-cphs0);

iter    = 1;
maxit   = 75;
resnorm = 1;
tol     = 1e-10;
alpha   = 0.50;

vmq_c0 = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T0+273.15)-1./(1200+273.15))); % Katz et al., 2003; Moore et al., 1998
vmq_c1 = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T0+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015
vmq0   = (1-SiO2m).*vmq_c0 + SiO2m.*vmq_c1;

while resnorm > tol && iter < maxit
    xqi = xq;  fqi = fq;

    [resx,resf,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTdx,perTdm,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);
    
    vmq0 = (1-cmq).*vmq_c0 + cmq.*vmq_c1;

    dresx_dx = cmq-cxq;

    if any(v(:)>10*TINY)
        dresf_df = vmq-vfq;
    else
        resf     = fq;
        dresf_df = ones(size(fq));
    end

    xq = xq - alpha.*resx./dresx_dx;
    fq = fq - alpha.*resf./dresf_df;

    resnorm = norm(xq-xqi)./sqrt(length(xq(:))) ...
            + norm(fq-fqi)./sqrt(length(fq(:)));

    iter    = iter+1;
end

[~,~,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTdx,perTdm,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);

xq = max(0,min(1-fq, xq ));
fq = max(0,min(1   , fq ));

end


function [resx,resf,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTdx,perTdm,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY)

if any(v(:)>10*TINY)
    vfq = ones(size(v));
    vmq = max(0,min(v./(1-xq),vmq0));
    
    Tphs0 = Tphs0d  - dTH2O(1).*vmq.^0.75;
    Tphs1 = Tphs1d  - dTH2O(3).*vmq.^0.75;
    perTx = (perTdx - dTH2O(2).*vmq.^0.75 -Tphs0)./(Tphs1-Tphs0);
    perTm = (perTdm - dTH2O(2).*vmq.^0.75 -Tphs0)./(Tphs1-Tphs0);
else
    vfq =  ones(size(v));
    vmq = zeros(size(v));
    
    Tphs0 = Tphs0d;
    Tphs1 = Tphs1d;
    perTx = (perTdx-Tphs0)./(Tphs1-Tphs0).*ones(size(v));
    perTm = (perTdm-Tphs0)./(Tphs1-Tphs0).*ones(size(v));
end

T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0)));

a = 15;
b = 0.05;

cx1 = max(0,min(1,          perCx .*erfc(PhDg(1).*(T-perTx)./(1-perTx))));
cx2 = max(0,min(1, perCx+(1-perCx).*erfc(PhDg(2).*(T-    0)./   perTx) ));

dcdT = -2*perCx*PhDg(1)/sqrt(pi)./(1-perTx);
cx1(T<perTx) = perCx + dcdT(T<perTx).*(T(T<perTx)-perTx(T<perTx));

ind1 = T>=perTx+b;
ind2 = T< perTx-b;
ind3 = T>=perTx-b & T<perTx+b;
cxq = zeros(size(T));
cxq(ind1) =  cx1(ind1);
cxq(ind2) =  cx2(ind2);
cxq(ind3) = (cx1(ind3).^-a+cx2(ind3).^-a).^-(1/a);

cm1 = max(0,min(1,          perCm .*erf(PhDg(3).*(1    -T)./(1-perTm))./erf(PhDg(3))));
cm2 = max(0,min(1, perCm+(1-perCm).*erf(PhDg(4).*(perTm-T)./(  perTm))./erf(PhDg(4))));

ind1 = T>=perTm+b;
ind2 = T< perTm-b;
ind3 = T>=perTm-b & T<perTm+b;
cmq = zeros(size(T));
cmq(ind1) =  cm1(ind1);
cmq(ind2) =  cm2(ind2);
cmq(ind3) = (cm1(ind3).^a+cm2(ind3).^a).^(1/a);

cxq = max(cphs0,min(c./(1-fq),cphs0 + cxq.*(cphs1-cphs0)));
cmq = min(cphs1,max(c./(1-fq),cphs0 + cmq.*(cphs1-cphs0)));

resx = c - (xq.*cxq + (1-xq-fq).*cmq);
resf = v - (fq.*vfq + (1-xq-fq).*vmq);

end