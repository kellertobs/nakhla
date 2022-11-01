% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Trace Elements  **************************************************

advn_TE = zeros(size(TE(inz,inx,:)));
diff_TE = zeros(size(TE(inz,inx,:)));
for i = 1:length(te0)
    
    % update trace element phase compositions
    tem(:,:,i)  = te(:,:,i)./(m + x.*Kte(i) );
    tex(:,:,i)  = te(:,:,i)./(m./Kte(i)  + x);

    % get trace element advection
    advn_TE(:,:,i) = - advect(rho(inz,inx).*m(inz,inx).*tem(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                     - advect(rho(inz,inx).*x(inz,inx).*tex(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    % get isotope ratio diffusion
    qz = - (kc(1:end-1,:)+kc(2:end,:))./2 .* ddz(te(:,:,i),h);
    qx = - (kc(:,1:end-1)+kc(:,2:end))./2 .* ddx(te(:,:,i),h);
    diff_TE(:,:,i) = (- ddz(qz(:,inx),h)  ...
                      - ddx(qx(inz,:),h));
end

% get total rate of change
dTEdt = advn_TE + diff_TE;

% update trace element concentrations
TE(inz,inx,:) = TEo(inz,inx,:) + (theta.*dTEdt + (1-theta).*dTEdto).*dt;   % explicit update
TE = max(0, TE );                                                          % enforce min bound
TE([1 end],:,:) = TE([2 end-1],:,:);                                       % boundary conditions
TE(:,[1 end],:) = TE(:,[2 end-1],:);


% *****  Isotope Ratios  **************************************************

advn_IR = zeros(size(IR(inz,inx,:)));
diff_IR = zeros(size(IR(inz,inx,:)));
for i = 1:length(ir0)
    % get isotope ratio advection
    advn_IR(:,:,i) = - advect(rho(inz,inx).*m(inz,inx).*ir(inz,inx,i),Um(inz,:),Wm(:,inx),h,{ADVN,''},[1,2],BCA) ...
                     - advect(rho(inz,inx).*x(inz,inx).*ir(inz,inx,i),Ux(inz,:),Wx(:,inx),h,{ADVN,''},[1,2],BCA);

    % get isotope ratio diffusion
    qz = - (kc(1:end-1,:)+kc(2:end,:))./2 .* ddz(ir(:,:,i),h);
    qx = - (kc(:,1:end-1)+kc(:,2:end))./2 .* ddx(ir(:,:,i),h);
    diff_IR(:,:,i) = (- ddz(qz(:,inx),h)  ...
                      - ddx(qx(inz,:),h));
end

% get total rate of change
dIRdt = advn_IR + diff_IR;

% update isotope ratio concentrations
IR(inz,inx,:) = IRo(inz,inx,:) + (theta.*dIRdt + (1-theta).*dIRdto).*dt;   % explicit update
IR([1 end],:,:) = IR([2 end-1],:,:);                                       % boundary conditions
IR(:,[1 end],:) = IR(:,[2 end-1],:);

for i = 1:length(te0); te(:,:,i) = TE(:,:,i)./rho; end
for i = 1:length(ir0); ir(:,:,i) = IR(:,:,i)./rho; end