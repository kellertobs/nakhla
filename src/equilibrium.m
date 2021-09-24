% calculate phase equilibrium

function [xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T0,c,v,P,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTv,PhDg)

TINY  = 0e-16;

perCm = (perCm-cphs0)./(cphs1-cphs0);
perCx = (perCx-cphs0)./(cphs1-cphs0);
perT  = (perT -Tphs0)./(Tphs1-Tphs0);

vmq   = (4.8e-5.*P.^0.6 + 1e-9.*P)./100;
vfq   = ones(size(v));
Kf    = vfq./vmq;

it    = 0;
maxit = 500;
res   = 1;
tol   = 1e-12;
alpha = 0.00;

if any(v>1e-6)
    
    while res > tol && it < maxit
        xi  = xq;  fi = fq;
               
        vm  = min(vmq,(v - fq.*vfq)./(1-fq-xq));
        
        T   = max(0,min(1,(T0 - P*clap + dTv.*vm.^0.75 -Tphs0)./(Tphs1-Tphs0)));
        
        cx1 = max(TINY,min(1-TINY,          perCx .*erfc((2+PhDg).*(T-perT)./(1-perT))));
        cx2 = max(TINY,min(1-TINY, perCx+(1-perCx).*erfc((1+PhDg).*(T     )./   perT) ));
        
        cxq = zeros(size(T));
        cxq(T>=perT) = cx1(T>=perT);
        cxq(T< perT) = cx2(T< perT);
        
        cm1 = max(TINY,min(1-TINY,          perCm .*erf((1.0+PhDg/10).*(1-T)         ./(1-perT))./erf((1.0+PhDg/10))));
        cm2 = max(TINY,min(1-TINY, perCm+(1-perCm).*erf((0.9+PhDg/10).*(1-T-(1-perT))./(  perT))./erf((0.9+PhDg/10))));
        
        cmq = zeros(size(T));
        cmq(T>=perT) = cm1(T>=perT);
        cmq(T< perT) = cm2(T< perT);
        
        cxq = min(c,cphs0 + cxq.*(cphs1-cphs0));
        cmq = max(c,cphs0 + cmq.*(cphs1-cphs0));
        
        xq  = max(TINY,min(1-TINY, alpha.*xi + (1-alpha) .* (c-(1-fq).*cmq)./(cxq-cmq) ));
        fq  = max(TINY,min(1-TINY, alpha.*fi + (1-alpha) .* (v-(1-xq).*vmq)./(vfq-vmq) ));
        
        res = (norm(xq(:)-xi(:),2) + norm(fq(:)-fi(:),2))./sqrt(2*length(xq(:)));
        it = it+1;
    end
    
    xq(T==1) = 0;
    xq(T==0) = 1;
    
    vmq = v./(fq.*Kf + (1-fq-xq));

else
    
    T   = max(0,min(1,(T0 - P*clap -Tphs0)./(Tphs1-Tphs0))) ;
    
    cx1 = max(TINY,min(1-TINY,          perCx .*erfc((2+PhDg).*(T-perT)./(1-perT))));
    cx2 = max(TINY,min(1-TINY, perCx+(1-perCx).*erfc((1+PhDg).*(T     )./   perT) ));
    
    cxq = zeros(size(T));
    cxq(T>=perT) = cx1(T>=perT);
    cxq(T< perT) = cx2(T< perT);
    
    cm1 = max(TINY,min(1-TINY,          perCm .*erf((1.0+PhDg/10).*(1-T)         ./(1-perT))./erf((1.0+PhDg/10))));
    cm2 = max(TINY,min(1-TINY, perCm+(1-perCm).*erf((0.9+PhDg/10).*(1-T-(1-perT))./(  perT))./erf((0.9+PhDg/10))));
    
    cmq = zeros(size(T));
    cmq(T>=perT) = cm1(T>=perT);
    cmq(T< perT) = cm2(T< perT);
    
    cxq = min(c,cphs0 + cxq.*(cphs1-cphs0));
    cmq = max(c,cphs0 + cmq.*(cphs1-cphs0));
    
    xq  = max(TINY,min(1-TINY, (c-cmq)./(cxq-cmq) ));
    fq  = 0.*xq;
    vmq = v./(fq.*Kf + (1-fq-xq));
    
    xq(T==1) = 0;
    xq(T==0) = 1;

end


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