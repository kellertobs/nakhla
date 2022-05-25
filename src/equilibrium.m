% calculate phase equilibrium

function [xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,PhDg,beta)

TINY  = 1e-16;

perCm = (perCm-cphs0 )./(cphs1-cphs0  );
perCx = (perCx-cphs0 )./(cphs1-cphs0  );

iter    = 0;
maxit   = 1e2;
resnorm = 1;
tol     = 1e-14;
eps     = 1e-6;

vmq_c0 = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T0+273.15)-1./(perTd+273.15))); % Katz et al., 2003; Moore et al., 1998
vmq_c1 = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T0+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015
vmq0   = (1-c).*vmq_c0 + c.*vmq_c1;

while resnorm > tol && iter < maxit
    xi  = xq;  fi = fq;
    
    [resx ,resf ] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg);
    
    [resx_xp,~,~,~,~,~] = res_xf(xq+eps,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg);
    [resx_xm,~,~,~,~,~] = res_xf(xq-eps,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg);
    
    dresx_dx = (resx_xp-resx_xm)/2/eps;
    
    if any(v(:)>0)
        [~,resf_fp,~,~,~,~] = res_xf(xq,fq+eps,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg);
        [~,resf_fm,~,~,~,~] = res_xf(xq,fq-eps,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg);
        
        dresf_df = (resf_fp-resf_fm)/2/eps;
    else
        resf     = 0*fq;
        dresf_df = TINY;
    end
    
    fq = max(TINY,min(1   -TINY, fi - (1-beta)*resf./dresf_df ));
    xq = max(TINY,min(1-fq-TINY, xi - (1-beta)*resx./dresx_dx ));
    
    resnorm = (norm(resx(:),2) + norm(resf(:),2))./sqrt(length(xq(:)));
    
    iter    = iter+1;
end

[~,~,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg);

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

function [resx,resf,cxq,cmq,vfq,vmq] = res_xf(xq,fq,T0,c,v,P,Tphs0d,Tphs1d,cphs0,cphs1,perTd,perCx,perCm,clap,dTH2O,vmq0,PhDg)

TINY = 1e-16;

if any(v(:)>0)
    vfq = ones(size(v));
    vmq = max(0,min(1,min(v./(1-xq+TINY),vmq0)));
    
    Tphs0 = Tphs0d - dTH2O(1).*vmq.^0.75;
    Tphs1 = Tphs1d - dTH2O(3).*vmq.^0.75;
    perT  = (perTd - dTH2O(2).*vmq.^0.75 -Tphs0)./(Tphs1-Tphs0);
else
    vfq = ones(size(v));
    vmq = zeros(size(v));
    
    Tphs0 = Tphs0d;
    Tphs1 = Tphs1d;
    perT  = (perTd-Tphs0)./(Tphs1-Tphs0);
end

T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0)));

cx1 = max(TINY,min(1-TINY,          perCx .*erfc(PhDg(1).*(T-perT)./(1-perT))));
cx2 = max(TINY,min(1-TINY, perCx+(1-perCx).*erfc(PhDg(2).*(T     )./   perT) ));

cxq = zeros(size(T));
cxq(T>=perT) = cx1(T>=perT);
cxq(T< perT) = cx2(T< perT);

cm1 = max(TINY,min(1-TINY,          perCm .*erf(PhDg(3).*(1   -T)./(1-perT))./erf(PhDg(3))));
cm2 = max(TINY,min(1-TINY, perCm+(1-perCm).*erf(PhDg(4).*(perT-T)./(  perT))./erf(PhDg(4))));

cmq = zeros(size(T));
cmq(T>=perT) = cm1(T>=perT);
cmq(T< perT) = cm2(T< perT);

cxq = min(c./(1-fq),cphs0 + cxq.*(cphs1-cphs0));
cmq = max(c./(1-fq),cphs0 + cmq.*(cphs1-cphs0));

resf = fq - max(TINY,min(1   -TINY, (v-(1-xq).*vmq)./(vfq-vmq) ));
resx = xq - max(TINY,min(1-fq-TINY, (c-(1-fq).*cmq)./(cxq-cmq) ));

end