% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TRC = zeros(Nz,Nx,cal.ntrc);
Ktrc    = zeros(Nz,Nx,cal.ntrc);
for i = 1:cal.ntrc
    
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Ktrc(:,:,i) = Ktrc(:,:,i) + cal.Ktrc_mem(i,j) .* cx_mem(:,:,j)./100; end

    % update trace element phase compositions
    trcm(:,:,i) = trc(:,:,i)./(m + x.*Ktrc(:,:,i));
    trcx(:,:,i) = trc(:,:,i)./(m./Ktrc(:,:,i) + x);

    % get trace element advection
    adv_TRC(:,:,i) = - advect(M.*trcm(:,:,i),Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
                     - advect(X.*trcx(:,:,i),Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

    % get trace element diffusion (regularisation)
    dff_TRC(:,:,i) = diffus(trcm(:,:,i),M.*kc,h,[1,2],BCD) + diffus(trcx(:,:,i),X.*kc,h,[1,2],BCD);

    % get trace element assimilation
    if ~isnan(trcwall(1,i)); bnd_TRC(:,:,i) = bnd_TRC(:,:,i) + (RHO.*trcwall(1,i)-TRC(:,:,i)).*mu./tau_a .* topshape; end
    if ~isnan(trcwall(2,i)); bnd_TRC(:,:,i) = bnd_TRC(:,:,i) + (RHO.*trcwall(2,i)-TRC(:,:,i)).*mu./tau_a .* botshape; end
    if ~isnan(trcwall(3,i)); bnd_TRC(:,:,i) = bnd_TRC(:,:,i) + (RHO.*trcwall(3,i)-TRC(:,:,i)).*mu./tau_a .* sdsshape; end
end

% get total rate of change
dTRCdt = adv_TRC + dff_TRC + bnd_TRC;

% residual of trace element evolution
res_TRC = (a1*TRC-a2*TRCo-a3*TRCoo)/dt - (b1*dTRCdt + b2*dTRCdto + b3*dTRCdtoo);

% semi-implicit update of trace element density
upd_TRC = - alpha*res_TRC*dt/a1 + beta*upd_TRC;
TRC     = TRC + upd_TRC;

% convert from densites to concentrations
for i = 1:cal.ntrc; trc(:,:,i) = TRC(:,:,i)./RHO; end
