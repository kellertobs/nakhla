% ***** TRACE & ISOTOPE GEOCHEMISTRY  *************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
itm = it./(m + x.*KIT);
itx = it./(m./KIT + x);

% update incompatible trace element composition
advn_IT = - advect(rho(inz,inx).*m(inz,inx).*itm(inz,inx),Um(inz,:),Wm(:,inx),h,SCHM,[1,2],BCA) ...  % advection
          - advect(rho(inz,inx).*x(inz,inx).*itx(inz,inx),Ux(inz,:),Wx(:,inx),h,SCHM,[1,2],BCA);

qz   = - (kc(1:end-1,:)+kc(2:end,:))/2 .* ddz(it,h);
qx   = - (kc(:,1:end-1)+kc(:,2:end))/2 .* ddx(it,h);
diff_IT = - ddz(qz(:,2:end-1),h) ...                                       % diffusion in melt
          - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(itwall); assim = assim + (itwall-it).*rho./tau_a .* bndshape; end% impose wall assimilation
    
dITdt = advn_IT + diff_IT + assim(inz,inx);                                % total rate of change

IT(inz,inx) = ITo(inz,inx) + (theta.*dITdt + (1-theta).*dITdto).*dt;       % explicit update
IT = max(0+TINY, IT );
IT([1 end],:) = IT([2 end-1],:);                                           % boundary conditions
IT(:,[1 end]) = IT(:,[2 end-1]);


% *****  COMPATIBLE TRACE ELEMENT  ****************************************

% update compatible trace element phase compositions
ctm = ct./(m + x.*KCT);
ctx = ct./(m./KCT + x);

% update compatible trace element composition
advn_CT = - advect(rho(inz,inx).*m(inz,inx).*ctm(inz,inx),Um(inz,:),Wm(:,inx),h,SCHM,[1,2],BCA) ...  % advection
          - advect(rho(inz,inx).*x(inz,inx).*ctx(inz,inx),Ux(inz,:),Wx(:,inx),h,SCHM,[1,2],BCA);

qz   = - (kc(1:end-1,:)+kc(2:end,:))/2 .* ddz(ct,h);
qx   = - (kc(:,1:end-1)+kc(:,2:end))/2 .* ddx(ct,h);
diff_CT = - ddz(qz(:,2:end-1),h) ...                                       % diffusion in melt
          - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(itwall); assim = assim + (itwall-it).*rho./tau_a .* bndshape; end% impose wall assimilation
    
dCTdt = advn_CT + diff_CT + assim(inz,inx);                                % total rate of change

CT(inz,inx) = CTo(inz,inx) + (theta.*dCTdt + (1-theta).*dCTdto).*dt;       % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);


% *****  STABLE ISOTOPE RATIO  ********************************************

% update stable isotope in melt
advn_SI = - advect(rho(inz,inx).*m(inz,inx).*si(inz,inx),Um(inz,:),Wm(:,inx),h,SCHM,[1,2],BCA) ...  % advection
          - advect(rho(inz,inx).*x(inz,inx).*si(inz,inx),Ux(inz,:),Wx(:,inx),h,SCHM,[1,2],BCA) ...
          - advect(rho(inz,inx).*f(inz,inx).*si(inz,inx),Uf(inz,:),Wf(:,inx),h,SCHM,[1,2],BCA);

qz   = - (kc(1:end-1,:)+kc(2:end,:))/2 .* ddz(si,h);
qx   = - (kc(:,1:end-1)+kc(:,2:end))/2 .* ddx(si,h);
diff_SI = - ddz(qz(:,2:end-1),h) ...                                       % diffusion in melt
          - ddx(qx(2:end-1,:),h);
                       
bndSI = zeros(size(si));
if ~isnan(siwall); bndSI = bndSI + (siwall-si).*rho.*m./tau_a .* bndshape; end % impose wall assimilation

dSIdt = advn_SI + diff_SI + assim(inz,inx);                                % total rate of change

SI(inz,inx) = SIo(inz,inx) + (theta.*dSIdt + (1-theta).*dSIdto).*dt;       % explicit update
SI([1 end],:) = SI([2 end-1],:);                                           % boundary conditions
SI(:,[1 end]) = SI(:,[2 end-1]);


% *****  RADIOGENIC ISOTOPES  *********************************************

% decay rate of radiogenic isotope
dcy_rip = rho.*rip./HLRIP.*log(2);
dcy_rid = rho.*rid./HLRID.*log(2);

% update radiogenic parent isotope phase compositions
ripm = rip./(m + x.*KRIP);
ripx = rip./(m./KRIP + x);

% update radiogenic parent isotope composition
advn_RIP = - advect(rho(inz,inx).*m(inz,inx).*ripm(inz,inx),Um(inz,:),Wm(:,inx),h,SCHM,[1,2],BCA) ...  % advection
           - advect(rho(inz,inx).*x(inz,inx).*ripx(inz,inx),Ux(inz,:),Wx(:,inx),h,SCHM,[1,2],BCA);

qz   = - (kc(1:end-1,:)+kc(2:end,:))/2 .* ddz(rip,h);
qx   = - (kc(:,1:end-1)+kc(:,2:end))/2 .* ddx(rip,h);
diff_RIP = - ddz(qz(:,2:end-1),h) ...                                      % diffusion in melt
           - ddx(qx(2:end-1,:),h);

assim = zeros(size(rip));
if ~isnan(riwall); assim = assim + (riwall-rip).*rho./tau_a .* bndshape; end % impose wall assimilation
                                                % secular equilibrium!
dRIPdt = advn_RIP + diff_RIP + assim(inz,inx) - dcy_rip(inz,inx) + dcy_rip(inz,inx);  % total rate of change
                                       
RIP(inz,inx) = RIPo(inz,inx) + (theta.*dRIPdt + (1-theta).*dRIPdto).*dt;   % explicit update
RIP = max(0+TINY, RIP );
RIP([1 end],:) = RIP([2 end-1],:);                                         % boundary conditions
RIP(:,[1 end]) = RIP(:,[2 end-1]);


% update radiogenic daughter isotope phase compositions
ridm = rid./(m + x.*KRID);
ridx = rid./(m./KRID + x);

% update radiogenic daughter isotope composition
advn_RID = - advect(rho(inz,inx).*m(inz,inx).*ridm(inz,inx),Um(inz,:),Wm(:,inx),h,SCHM,[1,2],BCA) ...  % advection
           - advect(rho(inz,inx).*x(inz,inx).*ridx(inz,inx),Ux(inz,:),Wx(:,inx),h,SCHM,[1,2],BCA);

qz   = - (kc(1:end-1,:)+kc(2:end,:))/2 .* ddz(rid,h);
qx   = - (kc(:,1:end-1)+kc(:,2:end))/2 .* ddx(rid,h);
diff_RID = - ddz(qz(:,2:end-1),h) ...                                      % diffusion in melt
           - ddx(qx(2:end-1,:),h);

assim = zeros(size(rid));
if ~isnan(riwall); assim = assim + (riwall.*HLRID./HLRIP-rid).*rho./tau_a .* bndshape; end % impose wall assimilation
    
dRIDdt = advn_RID + diff_RID + assim(inz,inx) - dcy_rid(inz,inx) + dcy_rip(inz,inx);  % total rate of change

RID(inz,inx) = RIDo(inz,inx) + (theta.*dRIDdt + (1-theta).*dRIDdto).*dt;   % explicit update
RID = max(0+TINY, RID );
RID([1 end],:) = RID([2 end-1],:);                                         % boundary conditions
RID(:,[1 end]) = RID(:,[2 end-1]);

it  = IT./rho;
ct  = CT./rho;
si  = SI./rho;
rip = RIP./rho;
rid = RID./rho;