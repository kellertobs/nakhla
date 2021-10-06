% calculate phase equilibrium

function [xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T0,c,v,P,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTv,PhDg)

TINY  = 1e-16;

perCm = (perCm-cphs0)./(cphs1-cphs0);
perCx = (perCx-cphs0)./(cphs1-cphs0);
perT  = (perT -Tphs0)./(Tphs1-Tphs0);

iter  = 0;
maxit = 500;
res   = 1;
tol   = 1e-15;
alpha = 0.50;

if any(v>1e-6)
    
    while res > tol && iter < maxit
        xi  = xq;  fi = fq;
               
        vmq = (4.8e-5.*P.^0.6 + 1e-9.*P)./100;
        vfq = ones(size(v));

%         vmq = min((v - fq.*vfq)./(1-fq-xq),vmq);
        vmq = min(v./(1-xq),vmq);

        T   = max(0,min(1,(T0 - P*clap + dTv.*vmq.^0.75 -Tphs0)./(Tphs1-Tphs0)));
        
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
        
        cmq = max(c./(1-fq),cphs0 + cmq.*(cphs1-cphs0));
        cxq = min(c./(1-fq),cphs0 + cxq.*(cphs1-cphs0));

        xq  = alpha.*xi + (1-alpha) .* max(TINY,min(1-fq-TINY, (c-(1-fq).*cmq)./(cxq-cmq) ));
        fq  = alpha.*fi + (1-alpha) .* max(TINY,min(1-xq-TINY, (v-(1-xq).*vmq)./(vfq-vmq) ));

%         xq(T<=0) = alpha.*xi(T<=0) + (1-alpha) .* 1-fq(T<=0);
%         xq(T>=1) = alpha.*xi(T>=1) + (1-alpha) .* 0;
        
%         fq  = max(TINY,min(1-xq-TINY, alpha.*fi + (1-alpha) .* (1-xq).*(v./(1-xq)-vmq)./(vfq-vmq) ));
%         xq  = max(TINY,min(1-fq-TINY, alpha.*xi + (1-alpha) .* (1-fq).*(c./(1-fq)-cmq)./(min(c./(1-fq),cxq)-cmq) ));
        
        res = (norm(xq(:)-xi(:),2) + norm(fq(:)-fi(:),2))./sqrt(2*length(xq(:)));
        iter = iter+1;
    end
    
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
    vfq = ones(size(v));
    vmq = zeros(size(v));

end

% cxq = min(c./(1-fq),cxq);
% cmq = max(c./(1-fq),cmq);

% vmq = vm;

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