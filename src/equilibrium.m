% calculate phase equilibrium

function [xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T0,c,v,P,cal,TINY)

Tphs0d = cal.Tphs0;
Tphs1d = cal.Tphs1;
cphs0  = cal.cphs0;
cphs1  = cal.cphs1;
perTd  = cal.perT;
perCx  = cal.perCx;
perCm  = cal.perCm;
clap   = cal.clap;
dTH2O  = cal.dTH2O;
PhDg   = cal.PhDg;

perCm = (perCm-cphs0)./(cphs1-cphs0);
perCx = (perCx-cphs0)./(cphs1-cphs0);

iter    = 1;
maxit   = 50;
resnorm = 1;
tol     = 1e-9;
eps     = 1e-6;

vmq_c0 = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T0+273.15)-1./(perTd+273.15))); % Katz et al., 2003; Moore et al., 1998
vmq_c1 = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T0+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015
vmq0   = (1-c).*vmq_c0 + c.*vmq_c1;

while resnorm > tol && iter < maxit
   
    beta = 1-exp(-iter/10);

    [resx,resf,~,cmq,~,~] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);
    
    [resx_xp,~,~,~,~,~] = res_xf(xq+eps,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);
    [resx_xm,~,~,~,~,~] = res_xf(xq-eps,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);
    
    dresx_dx = (resx_xp-resx_xm)/2/eps;

    if any(v(:)>10*TINY)
        vmq0   = (1-cmq).*vmq_c0 + cmq.*vmq_c1;

        [~,resf_fp,~,~,~,~] = res_xf(xq,fq+eps,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);
        [~,resf_fm,~,~,~,~] = res_xf(xq,fq-eps,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);
        
        dresf_df = (resf_fp-resf_fm)/2/eps;
    else
        resf     = fq - TINY;
        dresf_df = ones(size(fq));
    end
    
    xq = xq - beta*(resx./dresx_dx);
    fq = fq - beta*(resf./dresf_df);
    
    resnorm = (norm(resx(:),2) + norm(resf(:),2))/sqrt(length(xq(:)));
    
    iter    = iter+1;
end

xq = max(TINY,min(1-fq-TINY,xq));
fq = max(TINY,min(1-xq-TINY,fq));

[~,~,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY);

if size(xq,1)>1
    xq([1 end],:) = xq([2 end-1],:);
    xq(:,[1 end]) = xq(:,[2 end-1]);
    fq([1 end],:) = fq([2 end-1],:);
    fq(:,[1 end]) = fq(:,[2 end-1]);
    
    cxq([1 end],:) = cxq([2 end-1],:);
    cxq(:,[1 end]) = cxq(:,[2 end-1]);
    cmq([1 end],:) = cmq([2 end-1],:);
    cmq(:,[1 end]) = cmq(:,[2 end-1]);
    
    vfq([1 end],:) = vfq([2 end-1],:);
    vfq(:,[1 end]) = vfq(:,[2 end-1]);
    vmq([1 end],:) = vmq([2 end-1],:);
    vmq(:,[1 end]) = vmq(:,[2 end-1]);
end

end

function [resx,resf,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg,TINY)

if any(v(:)>10*TINY)
    vfq = ones(size(v));
    vmq = max(TINY,min(v./(1-xq),vmq0));
    
    Tphs0 = Tphs0d - dTH2O(1).*vmq.^0.75;
    Tphs1 = Tphs1d - dTH2O(3).*vmq.^0.75;
    perT  = (perTd - dTH2O(2).*vmq.^0.75 -Tphs0)./(Tphs1-Tphs0);
else
    vfq =  ones(size(v));
    vmq = zeros(size(v))+TINY;
    
    Tphs0 = Tphs0d;
    Tphs1 = Tphs1d;
    perT  = (perTd-Tphs0)./(Tphs1-Tphs0).*ones(size(v));
end

T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0)));

a = 20;
ind1 = T>=perT+0.02;
ind2 = T< perT-0.02;
ind3 = T>=perT-0.02 & T<perT+0.02;

cx1 = max(TINY,min(1-TINY,          perCx .*erfc(PhDg(1).*(T-perT)./(1-perT))));
cx2 = max(TINY,min(1-TINY, perCx+(1-perCx).*erfc(PhDg(2).*(T-   0)./   perT) ));

dcdT = -2*perCx*PhDg(1)/sqrt(pi)./(1-perT);
cx1(T<perT) = perCx + dcdT(T<perT).*(T(T<perT)-perT(T<perT));

cxq = zeros(size(T));
cxq(ind1) =  cx1(ind1);
cxq(ind2) =  cx2(ind2);
cxq(ind3) = (cx1(ind3).^-a+cx2(ind3).^-a).^-(1/a);

cm1 = max(TINY,min(1-TINY,          perCm .*erf(PhDg(3).*(1   -T)./(1-perT))./erf(PhDg(3))));
cm2 = max(TINY,min(1-TINY, perCm+(1-perCm).*erf(PhDg(4).*(perT-T)./(  perT))./erf(PhDg(4))));

cmq = zeros(size(T));
cmq(ind1) =  cm1(ind1);
cmq(ind2) =  cm2(ind2);
cmq(ind3) = (cm1(ind3).^a+cm2(ind3).^a).^(1/a);

cxq = min(min(cphs1,c./(1-fq)),cphs0 + cxq.*(cphs1-cphs0));
cmq = max(min(cphs1,c./(1-fq)),cphs0 + cmq.*(cphs1-cphs0));

resx = c - (xq.*cxq + (1-xq-fq).*cmq);
resf = v - (fq.*vfq + (1-xq-fq).*vmq);

end