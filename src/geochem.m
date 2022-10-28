% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TE = zeros(size(TE(inz,inx,:)));
adv_TE = zeros(size(TE(inz,inx,:)));
for k = 1:length(te0)
    
    % update trace element phase compositions
    tem(:,:,k)  = te(:,:,k)./(m + x.*Kte(k) );
    tex(:,:,k)  = te(:,:,k)./(m./Kte(k)  + x);

    % get trace element advection
    adv_TE(:,:,k) = - advect(rho(inz,inx).*m(inz,inx).*tem(inz,inx,k),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(rho(inz,inx).*x(inz,inx).*tex(inz,inx,k),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    % get trace element assimilation
    if ~isnan(tewall(k)); bnd_TE(:,:,k) = bnd_TE(:,:,k) + (tewall(k)-te(inz,inx,k)).*rho(inz,inx)./tau_a .* bndshape; end
end

% get total rate of change
dTEdt = adv_TE + bnd_TE;

% update trace element concentrations
TE(inz,inx,:) = TEo(inz,inx,:) + (theta.*dTEdt + (1-theta).*dTEdto).*dt;   % explicit update
TE = max(0, TE );                                                          % enforce min bound
TE([1 end],:,:) = TE([2 end-1],:,:);                                       % boundary conditions
TE(:,[1 end],:) = TE(:,[2 end-1],:);


% *****  Isotope Ratios  **************************************************

bnd_IR = zeros(size(IR(inz,inx,:)));
adv_IR = zeros(size(IR(inz,inx,:)));
for k = 1:length(ir0)
    % get isotope ratio advection
    adv_IR(:,:,k) = - advect(rho(inz,inx).*m(inz,inx).*ir(inz,inx,k),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(rho(inz,inx).*x(inz,inx).*ir(inz,inx,k),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    % get isotope ratio assimilation
    if ~isnan(irwall(k)); bnd_IR(:,:,k) = bnd_IR(:,:,k) + (irwall(k)-ir(inz,inx,k)).*rho(inz,inx)./tau_a .* bndshape; end
end

% get total rate of change
dIRdt = adv_IR + bnd_IR;

% update isotope ratio concentrations
IR(inz,inx,:) = IRo(inz,inx,:) + (theta.*dIRdt + (1-theta).*dIRdto).*dt;   % explicit update
IR([1 end],:,:) = IR([2 end-1],:,:);                                       % boundary conditions
IR(:,[1 end],:) = IR(:,[2 end-1],:);

for k = 1:length(te0); te(:,:,k) = TE(:,:,k)./rho; end
for k = 1:length(ir0); ir(:,:,k) = IR(:,:,k)./rho; end