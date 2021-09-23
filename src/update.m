% update solution-dependent parameter and auxiliary variable fields
eII(2:end-1,2:end-1) = 1e-16 + (0.5.*(exx(2:end-1,2:end-1).^2 + ezz(2:end-1,2:end-1).^2 ...  % get strain rate magnitude
    + 2.*(exz(1:end-1,1:end-1).^2.*exz(2:end,1:end-1).^2.*exz(1:end-1,2:end).^2.*exz(2:end,2:end).^2).^0.25)).^0.5;
eII(:,[1 end]) = eII(:,[2 end-1]);
eII([1 end],:) = eII([2 end-1],:);

tII(2:end-1,2:end-1) = 1e-16 + (0.5.*(txx(2:end-1,2:end-1).^2 + tzz(2:end-1,2:end-1).^2 ...  % get stress magnitude
    + 2.*(txz(1:end-1,1:end-1).^2.*txz(2:end,1:end-1).^2.*txz(1:end-1,2:end).^2.*txz(2:end,2:end).^2).^0.25)).^0.5;  
tII(:,[1 end]) = tII(:,[2 end-1]);
tII([1 end],:) = tII([2 end-1],:);

eta   = eta0 .* max(1e-6,1-phi/0.5).^-A .* max(1e-6,1-chi/0.5).^-B;            % bubble-crystal-dep. magma viscosity
eta   = 1./(1./(eta + eta0./sqrt(etactr)) + 1./(eta0.*sqrt(etactr)));
for d = 1:ceil(delta)
    dd  = delta/ceil(delta);
    eta = eta + dd.*(diff(eta([1,1:end,end],:),2,1)./8 + diff(eta(:,[1,1:end,end]),2,2)./8);
end
etac  = (eta(1:end-1,1:end-1)+eta(2:end,1:end-1) ...                       % viscosity in cell corners
       + eta(1:end-1,2:end  )+eta(2:end,2:end  ))./4;
    
rhoref = mean(mean(rho(2:end-1,2:end-1)));
rhoBF  = (rho(2:end-2,2:end-1)+rho(3:end-1,2:end-1))./2 - rhoref;          % relative density for bouancy force term

Cp    = mu.*Cm + chi.*Cx + phi.*Cf;                                        % bubble-crystal-dep. magma heat capacity

kT    = mu.*kTm + chi.*kTx + phi.*kTf;                                    % bubble-crystal-dep. magma thermal conductivity

wf = 2/9 .* (rhof-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*df^2./((eta(1:end-1,:)+eta(2:end,:))/2); % bubble flotation speed
wf([1 end],:) = 0;
for d = 1:1
    wf(2:end-1,:) = wf(2:end-1,:) + diff(wf,2,1)./8;
    wf([1 end],:) = 0;
end

wx = 2/9 .* (rhox-(rho(1:end-1,:)+rho(2:end,:))/2)*g0*dx^2./((eta(1:end-1,:)+eta(2:end,:))/2); % crystal foundering speed
wx([1 end],:) = 0;
for d = 1:1
    wx(2:end-1,:) = wx(2:end-1,:) + diff(wx,2,1)./8;
    wx([1 end],:) = 0;
end
