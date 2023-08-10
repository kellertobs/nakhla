% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TE = zeros(Nz,Nx,cal.nte);
adv_TE = zeros(Nz,Nx,cal.nte);
dff_TE = zeros(Nz,Nx,cal.nte);
Kte    = zeros(Nz,Nx,cal.nte);
for i = 1:cal.nte
    
    % update bulk partitioning coefficients
    for j=1:cal.nmem; Kte(:,:,i) = Kte(:,:,i) + cal.Kte_mem(i,j) .* c_mem(:,:,j)./100; end

    % update trace element phase compositions
    tem(:,:,i) = te(:,:,i)./(m + x.*Kte(:,:,i));
    tex(:,:,i) = te(:,:,i)./(m./Kte(:,:,i) + x);

    % get trace element advection
    adv_TE(:,:,i) = - advect(M.*tem(:,:,i),Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X.*tex(:,:,i),Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

    % get trace element diffusion (regularisation)
    dff_TE(:,:,i) = diffus(te(:,:,i),kc,h,[1,2],BCD);

    % get trace element assimilation
    if ~isnan(tewall(1,i)); bnd_TE(:,:,i) = bnd_TE(:,:,i) + (rho.*tewall(1,i).*(1-f)-TE(:,:,i)).*mu./tau_a .* topshape; end
    if ~isnan(tewall(2,i)); bnd_TE(:,:,i) = bnd_TE(:,:,i) + (rho.*tewall(2,i).*(1-f)-TE(:,:,i)).*mu./tau_a .* botshape; end
    if ~isnan(tewall(3,i)); bnd_TE(:,:,i) = bnd_TE(:,:,i) + (rho.*tewall(3,i).*(1-f)-TE(:,:,i)).*mu./tau_a .* sdsshape; end
end

% get total rate of change
dTEdt = adv_TE + dff_TE + bnd_TE;

% residual of trace element evolution
res_TE = (a1*TE-a2*TEo-a3*TEoo)/dt - (b1*dTEdt + b2*dTEdto + b3*dTEdtoo);

% update trace element concentrations
TE = (a2*TEo+a3*TEoo + (b1*dTEdt + b2*dTEdto + b3*dTEdtoo)*dt)/a1;

% enforce min bound
TE = max(0, TE );                                                          


% *****  Isotope Ratios  **************************************************

bnd_IR = zeros(Nz,Nx,cal.nir);
adv_IR = zeros(Nz,Nx,cal.nir);
dff_IR = zeros(Nz,Nx,cal.nir);
for i = 1:cal.nir

    % update trace element phase compositions
    irm(:,:,i) = ir(:,:,i)./(1-f);
    irx(:,:,i) = ir(:,:,i)./(1-f);

    % get isotope ratio advection
    adv_IR(:,:,i) = - advect(M.*irm(:,:,i),Um(2:end-1,:),Wm(:,2:end-1),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X.*irx(:,:,i),Ux(2:end-1,:),Wx(:,2:end-1),h,{ADVN,''},[1,2],BCA);

    % get trace element diffusion (regularisation)
    dff_IR(:,:,i) = diffus(ir(:,:,i),kc,h,[1,2],BCD);

    % get isotope ratio assimilation
    if ~isnan(irwall(1,i)); bnd_IR(:,:,i) = bnd_IR(:,:,i) + (rho.*irwall(1,i).*(1-f)-IR(:,:,i)).*mu./tau_a .* topshape; end
    if ~isnan(irwall(2,i)); bnd_IR(:,:,i) = bnd_IR(:,:,i) + (rho.*irwall(2,i).*(1-f)-IR(:,:,i)).*mu./tau_a .* botshape; end
    if ~isnan(irwall(3,i)); bnd_IR(:,:,i) = bnd_IR(:,:,i) + (rho.*irwall(3,i).*(1-f)-IR(:,:,i)).*mu./tau_a .* sdsshape; end
end

% get total rate of change
dIRdt = adv_IR + dff_IR + bnd_IR;

% residual of isotope evolution
res_IR = (a1*IR-a2*IRo-a3*IRoo)/dt - (b1*dIRdt + b2*dIRdto + b3*dIRdtoo);

% update isotope ratio concentrations
IR = (a2*IRo+a3*IRoo + (b1*dIRdt + b2*dIRdto + b3*dIRdtoo)*dt)/a1;


% convert from densites to concentrations
for i = 1:cal.nte; te(:,:,i) = TE(:,:,i)./RHO; end
for i = 1:cal.nir; ir(:,:,i) = IR(:,:,i)./RHO; end