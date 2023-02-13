% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

TEii = TEi;  TEi = TE;
IRii = IRi;  IRi = IR;

% *****  Trace Elements  **************************************************

bnd_TE = zeros(Nz,Nx,cal.nte);
adv_TE = zeros(Nz,Nx,cal.nte);
Kte    = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    
    % update bulk partitioning coefficients
    for j=1:cal.nc; Kte(:,:,i) = Kte(:,:,i) + cal.Kte_cmp(i,j) .* c_cmp(:,:,j)./100; end

    % update trace element phase compositions
    tem(:,:,i) = te(:,:,i)./(m + x.*Kte(:,:,i));
    tex(:,:,i) = te(:,:,i)./(m./Kte(:,:,i) + x);

    % get trace element advection
    adv_TE(:,:,i) = - advect(M.*tem(:,:,i),Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X.*tex(:,:,i),Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

    % get trace element assimilation
    if ~isnan(tewall(i)); bnd_TE(:,:,i) = bnd_TE(:,:,i) + (rho.*tewall(i)-TE(:,:,i))./tau_a .* bndshape; end
end

% get total rate of change
dTEdt = adv_TE + bnd_TE;

% residual of trace element evolution
res_TE = (a1*TE-a2*TEo-a3*TEoo)/dt - (b1*dTEdt + b2*dTEdto + b3*dTEdtoo);

% update trace element concentrations
TE = TE - alpha*res_TE*dt + beta*(TEii-TEi);
TE = max(0, TE );                                                          % enforce min bound


% *****  Isotope Ratios  **************************************************

bnd_IR = zeros(Nz,Nx,cal.nir);
adv_IR = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir

    % update trace element phase compositions
    irm(:,:,i) = ir(:,:,i)./(1-f);
    irx(:,:,i) = ir(:,:,i)./(1-f);

    % get isotope ratio advection
    adv_IR(:,:,i) = - advect(M.*irm(:,:,i),Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X.*irx(:,:,i),Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

    % get isotope ratio assimilation
    if ~isnan(irwall(i)); bnd_IR(:,:,i) = bnd_IR(:,:,i) + (rho.*irwall(i)-IR(:,:,i))./tau_a .* bndshape; end
end

% get total rate of change
dIRdt = adv_IR + bnd_IR;

% residual of isotope ratio evolution
res_IR = (a1*IR-a2*IRo-a3*IRoo)/dt - (b1*dIRdt + b2*dIRdto + b3*dIRdtoo);

% update isotope ratio concentrations
IR = IR - alpha*res_IR*dt + beta*(IRii-IRi);

% convert from mixture density to concentration
for i = 1:cal.nte; te(:,:,i) = TE(:,:,i)./rho; end
for i = 1:cal.nir; ir(:,:,i) = IR(:,:,i)./rho; end