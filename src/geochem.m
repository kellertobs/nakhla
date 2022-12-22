% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

bnd_TE = zeros(size(TE(inz,inx,:)));
adv_TE = zeros(size(TE(inz,inx,:)));
for i = 1:length(te0)
    
    % update trace element phase compositions
    tem(:,:,i)  = te(:,:,i)./(m + x.*Kte(i) );
    tex(:,:,i)  = te(:,:,i)./(m./Kte(i)  + x);

    % get trace element advection
    adv_TE(:,:,i) = - advect(M(inz,inx).*tem(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X(inz,inx).*tex(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    % get trace element assimilation
    if ~isnan(tewall(i)); bnd_TE(:,:,i) = bnd_TE(:,:,i) + (tewall(i)-te(inz,inx,i)).*rho(inz,inx)./tau_a .* bndshape; end
end

% get total rate of change
dTEdt = adv_TE + bnd_TE;

% update trace element concentrations
TE(inz,inx,:)   = TEo(inz,inx,:) + (theta.*dTEdt + (1-theta).*dTEdto).*dt; % explicit update
TE = max(0, TE );                                                          % enforce min bound
TE([1 end],:,:) = TE([2 end-1],:,:);                                       % boundary conditions
TE(:,[1 end],:) = TE(:,[2 end-1],:);


% *****  Isotope Ratios  **************************************************

bnd_IR = zeros(size(IR(inz,inx,:)));
adv_IR = zeros(size(IR(inz,inx,:)));
for i = 1:length(ir0)

    % update trace element phase compositions
    irm(:,:,i)  = ir(:,:,i)./(m + x*1);
    irx(:,:,i)  = ir(:,:,i)./(m/1 + x);

    % get isotope ratio advection
    adv_IR(:,:,i) = - advect(M(inz,inx).*irm(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                    - advect(X(inz,inx).*irx(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    % get isotope ratio assimilation
    if ~isnan(irwall(i)); bnd_IR(:,:,i) = bnd_IR(:,:,i) + (irwall(i)-ir(inz,inx,i)).*rho(inz,inx)./tau_a .* bndshape; end
end

% get total rate of change
dIRdt = adv_IR + bnd_IR;

% update isotope ratio concentrations
IR(inz,inx,:)   = IRo(inz,inx,:) + (theta.*dIRdt + (1-theta).*dIRdto).*dt; % explicit update
IR([1 end],:,:) = IR([2 end-1],:,:);                                       % boundary conditions
IR(:,[1 end],:) = IR(:,[2 end-1],:);

% convert from mixture density to concentration
for i = 1:length(te0); te(:,:,i) = TE(:,:,i)./rho; end
for i = 1:length(ir0); ir(:,:,i) = IR(:,:,i)./rho; end