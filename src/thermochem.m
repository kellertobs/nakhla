% thermochemical solver

% update time step
dtk = min((h/2)^2./max([kT(:)./rhoCp(:);kc;kv]));                          % diffusive time step size
dth = 1e32;%1./max(abs((dHdt(:)./rhoCp(:)))+1e-16);                        % compaction time step size
% dth = 0.01./max(abs(x(:)-xo(:))+1e-16);                                  % compaction time step size
dta = min(min(h/2/max(abs([UBG(:);WBG(:);Um(:);Wm(:);Uf(:);Wf(:);Ux(:);Wx(:)]+1e-16)))); % advective time step size
dt  = min(2*dto,CFL*min([dtk,dth,dta]));                                   % physical time step size


%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% update temperature
if ~isotherm

    advn_H = advection(rhom.*mu .*(Cpm.*T + 0  ),Um,Wm,h,ADVN,'flx') ...
           + advection(rhox.*chi.*(Cpx.*T + DLx),Ux,Wx,h,ADVN,'flx') ...
           + advection(rhof.*phi.*(Cpf.*T + DLf),Uf,Wf,h,ADVN,'flx');
                           
    qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                 % heat diffusion z-flux
    qTx    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddx(T,h);                 % heat diffusion x-flux
    diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...                 % heat diffusion
                               - ddx(qTx(2:end-1,:),h));
    
    switch coolmode
        case 0
            coolshape = zeros(size(T));
        case 1
            coolshape = exp( ( -ZZ)/dw);
        case 2
            coolshape = exp( ( -ZZ)/dw) ...
                + exp(-(D-ZZ)/dw);
        case 3
            coolshape = exp( ( -ZZ)/dw) ...
                + exp(-(D-ZZ)/dw) ...
                + exp( ( -XX)/dw) ...
                + exp(-(L-XX)/dw);
    end
    coolshape = min(1,coolshape);
    cool = rhoCp .* (Twall-T)./tau_T .* coolshape;
        
    dHdt   = - advn_H + diff_T + cool;                                     % total rate of change
    
    if step>0; H = Ho + (theta.*dHdt + (1-theta).*dHdto).*dt; end          % explicit update of enthalpy
    H([1 end],:) = H([2 end-1],:);                                         % apply boundary conditions
    H(:,[1 end]) = H(:,[2 end-1]);    
    
end


% update composition
if ~isochem
  
    % update major component
    advn_C = advection(rhom.*mu .*cm,Um,Wm,h,ADVN,'flx') ...
           + advection(rhox.*chi.*cx,Ux,Wx,h,ADVN,'flx');
                          
    qcz   = - kc.*rhom.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(cm,h);        % major component diffusion z-flux
    qcx   = - kc.*rhom.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(cm,h);        % major component diffusion x-flux
    diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                  % major component diffusion
                              - ddx(qcx(2:end-1,:),h);
    
    assim = zeros(size(c));
    if ~isnan(ctop); assim = assim + rhob .* (ctop-c)./tau_c.* exp((-ZZ)/dw); end % impose top assimilation rate
    
    dCdt = - advn_C + diff_c + assim;                                      % total rate of change
    
    if step>0; C = Co + (theta.*dCdt + (1-theta).*dCdto).*dt; end          % explicit update of major component density
    C = min(rho.*cphs1-TINY,max(rho.*cphs0+TINY, C ));
    C([1 end],:) = C([2 end-1],:);                                         % apply boundary conditions
    C(:,[1 end]) = C(:,[2 end-1]);  
    
    % update volatile component
    advn_V = advection(rhom.*mu .*vm,Um,Wm,h,ADVN,'flx') ...
           + advection(rhof.*phi.*vf,Uf,Wf,h,ADVN,'flx');
        
    qvz   = - kc.*rhom.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(vm,h);        % volatile component diffusion z-flux
    qvx   = - kc.*rhom.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(vm,h);        % volatile component diffusion x-flux
    diff_v(2:end-1,2:end-1) = - ddz(qvz(:,2:end-1),h) ...                  % volatile component diffusion
                              - ddx(qvx(2:end-1,:),h);
    
    degas = zeros(size(v));
    if ~isnan(ftop); degas = degas + min(0,ftop-f).*rho./tau_f.* exp((-ZZ)/dw); end % impose top degassing rate
    
    dVdt = - advn_V + diff_v + degas;                                      % total rate of change
    
    if step>0; V = Vo + (theta.*dVdt + (1-theta).*dVdto).*dt; end          % explicit update of volatile component density
    V = min(rho.*1-TINY,max(rho.*0+TINY, V ));
    V([1 end],:) = V([2 end-1],:);                                         % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);  
    
end


% convert enthalpy and component densities to temperature and concentrations
T = (H - chi.*rhox.*DLx - phi.*rhof.*DLf)./rhoCp;
c = C./rho;
v = V./rho;


%% *****  update local phase equilibrium  *********************************

[xqi,cxq,cmq,fqi,vfq,vmq] = equilibrium(x,f,(T+To)./2,(c+co)./2,(v+vo)./2,(Pt+Pto)./2,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);
% [xqi,cxqi,cmqi,fqi,vfqi,vmqi] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);

xq  = beta.*xq + (1-beta).*xqi;
fq  = beta.*fq + (1-beta).*fqi;

    
% update crystal fraction
if diseq
    
    Gxi = (xq-x)./max(4.*dt,tau_r);
    Gx  = beta.*Gx + (1-beta).*Gxi;                                        % crystallisation rate
    
    advn_x = advection(x,Ux,Wx,h,ADVN,'flx');                              % get advection term
    
    dxdt   = - advn_x + Gx;                                                % total rate of change
    
    if step>0; x = x + (theta.*dxdt + (1-theta).*dxdto).*dt; end           % explicit update of crystal fraction
    x = min(1-f-TINY,max(TINY,x));                                         % enforce [0,1] limit
    x([1 end],:) = x([2 end-1],:);                                         % apply boundary conditions
    x(:,[1 end]) = x(:,[2 end-1]);
    
    cx  = cxq;
    cm  = max(c+TINY,min(cphs1-TINY,(c - x.*cx)./(1-x-f)));
    cm(m < 1e-3) = cmq(m < 1e-3);
    
else
    
    x  = xq;                                                               % enforce local phase equilibrium
    cx = cxq;
    cm = cmq;
    Gx = (x-xo)./dt + advection(x,Ux,Wx,h,ADVN,'flx');                     % reconstruct crystallisation rate
    
end


% update bubble fraction
if diseq
    
    Gfi = (fq-f)./max(4.*dt,tau_r);
    Gf  = beta.*Gf + (1-beta).*Gfi;
    
    advn_f = advection(f,Uf,Wf,h,ADVN,'flx');                              % get advection term
    
    dfdt   = - advn_f + Gf;                                                % total rate of change
    
    if step>0; f = fo + (theta.*dfdt + (1-theta).*dfdto).*dt; end          % explicit update of bubble fraction
    f = min(1-x-TINY,max(TINY,f));                                         % enforce [0,1-x] limit
    f([1 end],:) = f([2 end-1],:);                                         % apply boundary conditions
    f(:,[1 end]) = f(:,[2 end-1]);
    
    vf  = vfq;
    vm  = max(v+TINY,min(1 -TINY,(v - f.*vf)./(1-f-x)));
    vm(m < 1e-3) = vmq(m < 1e-3);
    
else
    
    f  = fq;
    vf = vfq;
    vm = vmq;
    Gf = (f-fo)./dt + advection(f,Uf,Wf,h,ADVN,'flx');                     % reconstruct exsolution rate
    
end


% update melt fraction
m = 1-x-f;


%% ***** Trace & isotope geochemistry  ************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
itm = it./(m + x.*KIT);
itx = it./(m./KIT + x);

% update incompatible trace element composition
advn_IT = advection(rhom.*mu .*itm,Um,Wm,h,ADVN,'flx') ...
        + advection(rhox.*chi.*itx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*rhom.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(itm,h);
qx   = - kc.*rhom.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(itm,h);
diff_it(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

dITdt = - advn_IT + diff_it;                                               % total rate of change

if step>0; IT = ITo + (theta.*dITdt + (1-theta).*dITdto).*dt; end          % explicit update
IT = max(0+TINY, IT );
IT([1 end],:) = IT([2 end-1],:);                                           % boundary conditions
IT(:,[1 end]) = IT(:,[2 end-1]);

it = IT./rho;


% *****  COMPATIBLE TRACE ELEMENT  ****************************************

% update compatible trace element phase compositions
ctm = ct./(m + x.*KCT);
ctx = ct./(m./KCT + x);

% update compatible trace element composition
advn_CT = advection(rhom.*mu .*ctm,Um,Wm,h,ADVN,'flx') ...
        + advection(rhox.*chi.*ctx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*rhom.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(ctm,h);
qx   = - kc.*rhom.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(ctm,h);
diff_ct(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

dCTdt = - advn_CT + diff_ct;                                               % total rate of change

if step>0; CT = CTo + (theta.*dCTdt + (1-theta).*dCTdto).*dt; end          % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);

ct = CT./rho;


% *****  STABLE ISOTOPE  **************************************************

% update stable isotope composition
advn_si = advection(mu .*si,Um,Wm,h,ADVN,'flx') ...
        + advection(chi.*si,Ux,Wx,h,ADVN,'flx') ...
        + advection(phi.*si,Uf,Wf,h,ADVN,'flx');

qz   = - kc.*(mu(1:end-1,:)+mu(2:end,:))/2 .* ddz(si,h);
qx   = - kc.*(mu(:,1:end-1)+mu(:,2:end))/2 .* ddx(si,h);
diff_si(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

dsidt = - advn_si + diff_si;                                               % total rate of change

if step>0; si = sio + (theta.*dsidt + (1-theta).*dsidto).*dt; end          % explicit update
si([1 end],:) = si([2 end-1],:);                                           % boundary conditions
si(:,[1 end]) = si(:,[2 end-1]);


% if step == 0; dHdto = dHdt;  dCdto = dCdt;  dVdto = dVdt;  dxdto = dxdt;  dfdto = dfdt;  dITdto = dITdt;  dCTdto = dCTdt;  dSIdto = dSIdt; end


%     % update time step
%     Vel    = [U(:)+u(:);W(:)+w(:);fUf(:);fWf(:);fUs(:);fWs(:)];            % combine all velocity components
%     DfDt   = Div_fV+RctR_f(ic,ic);
% %     Div_fV = diff(fUs(ic,:),1,2)./h + diff(fWs(:,ic),1,1)./h;            % phase advection/compaction
%     Div_fVi = diff(fUs(ic,:),1,2)./h + diff(fWs(:,ic),1,1)./h/2 ...
%             - diff(fUf(ic,:),1,2)./h - diff(fWf(:,ic),1,1)./h/2;           % phase advection/compaction
%     Div_fVi(f(ic,ic)<=flim | (1-f(ic,ic))<=flim) = 0;
%     Div_fV = (1-alpha).*Div_fVi + alpha.*Div_fV;
%
%     PeTeff = PeT .* (1 - 0.75./(1+exp(-(f-0.4).*16)));
%     PeCeff = PeC .* (1 - 0.75./(1+exp(-(f-0.4).*16)));
%
%     dt     = min(2*dto,CFL*min([ (h/2)^2*min(PeTeff(:)) , h/2/max(abs(Vel)) , 0.005./max(abs(DfDt(:)))]));   % ,  physical time step
%
%     % update temperature
%     advn_T = flxdiv(T,fUs,fWs,h,ADVN,'adv') + flxdiv(T,fUf,fWf,h,ADVN,'adv'); % advection
%
%     cmpt_T = -1/St.*Div_fV;                                                % advection/compaction
%
%     diff_T = diff( (1./PeTeff(im,ic)+1./PeTeff(ip,ic))/2 .* diff(T(:,ic),1,1) ,1,1)./h^2 ...
%            + diff( (1./PeTeff(ic,im)+1./PeTeff(ic,ip))/2 .* diff(T(ic,:),1,2) ,1,2)./h^2;
%
%     RTin   = (Tin-T)./(5*dt) .* exp(-(L-Z)./(2*h));                        % base injection rate
%
%     RT(ic,ic) = - advn_T - cmpt_T + diff_T + RTin(ic,ic);                  % total rate of change
%
%     res_T = ((T + f/St)-(To + fo/St))./dt - (theta.*RT + (1-theta).*RTo);  % residual temperature evolution equation
%
%     T     = T - res_T.*dt/4;                                               % update temperature solution
%
%     if isotherm_topbot                                                     % apply boundary conditions
%         T([1 end],:) = To([1 end],:);
%     else
%         T([1 end],:) = T(ibz,:);
%     end
%     if isotherm_sides
%         T(:,[1 end]) = To(:,[1 end]);
%     else
%         T(:,[1 end]) = T(:,ibx);
%     end
%     res_T([1 end],:) = 0;
%     res_T(:,[1 end]) = 0;
%
%
%     % update major element composition
%     advn_MAJ = flxdiv(MAJs,fUs,fWs,h,ADVN,'adv') + flxdiv(MAJf,fUf,fWf,h,ADVN,'adv'); % advection
%
%     cmpt_MAJ = (MAJs(ic,ic)-MAJf(ic,ic)).*Div_fV;                          % advection/compaction
%
%     diff_MAJ = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(MAJ(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(MAJ(ic,:),1,2) ,1,2)./h^2;
%
%     RMAJin   = (MAJin-MAJ)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate
%
%     RMAJ(ic,ic) = - advn_MAJ - cmpt_MAJ + diff_MAJ + RMAJin(ic,ic);        % total rate of change
%
%     res_MAJ = (MAJ-MAJo)./dt - (theta.*RMAJ + (1-theta).*RMAJo);           % residual composition evolution equation
%
%     MAJ     = MAJ - res_MAJ.*dt/4;                                         % update composition solution
%     MAJ     = max(1e-16,min(MAJf,MAJ));                                    % enforce min/max bounds
%
%     MAJ([1 end],:) = MAJ(ibz,:);                                           % apply boundary conditions
%     MAJ(:,[1 end]) = MAJ(:,ibx);
%     res_MAJ([1 end],:) = 0;
%     res_MAJ(:,[1 end]) = 0;
%
%
%     %*****  PHASE EQUILIBRIUM  ************************************************
%
%     % update equilibrium
%     [fq,MAJsq,MAJfq]  =  equilibrium(T,MAJ,Pt,perT,perCs,perCf,clap,PhDg);
%
%     if diseq
%
%         % update reaction rate
%         RctR_fi = min(Da,0.5/dt).*(fq-f);
%         for k  = 1:ceil(kappa)                                             % regularisation
%             kk = kappa/ceil(kappa);
%             RctR_fi(ic,ic) = RctR_fi(ic,ic) + kk.*(diff(RctR_fi(ic,:),2,2)+diff(RctR_fi(:,ic),2,1))./8;
%             RctR_fi([1 end],:) = RctR_fi(ibz,:);
%             RctR_fi(:,[1 end]) = RctR_fi(:,ibx);
%         end
%         RctR_f = (1-alpha)*RctR_fi + alpha*RctR_f;
%
%         % update disequilibrium melt fraction
%         Rfin   = (fin-f)./(5*dt) .* exp(-(L-Z)./(2*h));                     % base injection rate
%
%         Rf(ic,ic) = + Div_fV + RctR_f(ic,ic) + Rfin(ic,ic);                % total rate of change
%
%         res_f = (f-fo)./dt - (theta.*Rf + (1-theta).*Rfo);                 % residual composition evolution equation
%
%         f     = f - res_f.*dt/4;                                           % update composition solution
%         f     = max(1e-16,min(1-1e-16,f));                                 % enforce min/max bounds
%
%         qind = f<=flim | fq<=flim | (1-f)<=flim | (1-fq)<=flim;
%         f(qind)      = fq(qind);                                           % apply equilibrium outside cutoff
%
%         f([1 end],:) = f(ibz,:);                                           % apply boundary conditions
%         f(:,[1 end]) = f(:,ibx);
%         res_f([1 end],:) = 0;
%         res_f(:,[1 end]) = 0;
%
%     else
%
%         % update equilibrium melt fraction
%         f              =  (1-alpha)*fq + alpha*f;
%
%         % update reaction rate
%         RctR_fi        =  0.*RctR_f;
%         RctR_fi(ic,ic) = (f(ic,ic)-fo(ic,ic))./dt - Div_fV;
%         RctR_fi(f<=flim | (1-f)<=flim) = 0;
%         for k  = 1:ceil(kappa)                                             % regularisation
%             kk = kappa/ceil(kappa);
%             RctR_fi(ic,ic) = RctR_fi(ic,ic) + kk.*(diff(RctR_fi(ic,:),2,2)+diff(RctR_fi(:,ic),2,1))./8;
%             RctR_fi([1 end],:) = RctR_fi(ibz,:);
%             RctR_fi(:,[1 end]) = RctR_fi(:,ibx);
%         end
%         RctR_f = RctR_fi;
%
%     end
%
%     % update phase major element composition
%     KMAJ = MAJsq./MAJfq;
%     MAJf = MAJ./(f + (1-f).*KMAJ);
%     MAJs = MAJ./(f./KMAJ + (1-f));
%     MAJs(f<=  flim) = MAJ(f<=  flim);  MAJs(f>=1-flim) = MAJsq(f>=1-flim);
%     MAJf(f>=1-flim) = MAJ(f>=1-flim);  MAJf(f<=  flim) = MAJfq(f<=  flim);
%
%     % update incompatible trace element phase compositions
%     TRIf = TRI./(f + (1-f).*KTRI);
%     TRIs = TRI./(f./KTRI + (1-f));
%     TRIs(f<=  flim) = TRI(f<=  flim);
%     TRIf(f>=1-flim) = TRI(f>=1-flim);
%
%     % update compatible trace element phase compositions
%     TRCf = TRC./(f + (1-f).*KTRC);
%     TRCs = TRC./(f./KTRC + (1-f));
%     TRCs(f<=  flim) = TRC(f<=  flim);
%     TRCf(f>=1-flim) = TRC(f>=1-flim);
%
%     % update radiogenic parent isotope phase compositions
%     IRPf = IRP./(f + (1-f).*KIRP);
%     IRPs = IRP./(f./KIRP + (1-f));
%     IRPs(f<=  flim) = IRP(f<=  flim);
%     IRPf(f>=1-flim) = IRP(f>=1-flim);
%
%     % update radiogenic daughter isotope phase compositions
%     IRDf = IRD./(f + (1-f).*KIRD);
%     IRDs = IRD./(f./KIRD + (1-f));
%     IRDs(f<=  flim) = IRD(f<=  flim);
%     IRDf(f>=1-flim) = IRD(f>=1-flim);
%
%     % radiogenic isotope decay rate
%     DcyR_IRP = IRP./DIRP.*log(2);
%     DcyR_IRD = IRD./DIRD.*log(2);
%
%     % stable isotope transfer composition
%     ISR = double(RctR_f>0).*ISS + double(RctR_f<=0).*ISF;
%
%
%     % *****  TRACE ELEMENTS  **********************************************
%
%     % update incompatible trace element composition
%     advn_TRI = flxdiv(TRIs,fUs,fWs,h,ADVN,'adv') + flxdiv(TRIf,fUf,fWf,h,ADVN,'adv'); % advection
%
%     cmpt_TRI = (TRIs(ic,ic)-TRIf(ic,ic)).*Div_fV;                          % advection/compaction
%
%     diff_TRI = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(TRI(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(TRI(ic,:),1,2) ,1,2)./h^2;
%
%     RTRIin   = (TRIin-TRI)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate
%
%     RTRI(ic,ic) = - advn_TRI - cmpt_TRI + diff_TRI + RTRIin(ic,ic);    % total rate of change
%
%     res_TRI = (TRI-TRIo)./dt - (theta.*RTRI + (1-theta).*RTRIo);           % residual composition evolution equation
%
%     TRI     = TRI - res_TRI.*dt/4;                                         % update composition solution
%     TRI     = max(1e-16,min(TRIf,TRI));                                    % enforce min/max bounds
%
%     TRI([1 end],:) = TRI(ibz,:);                                           % apply boundary conditions
%     TRI(:,[1 end]) = TRI(:,ibx);
%     res_TRI([1 end],:) = 0;
%     res_TRI(:,[1 end]) = 0;
%
%
%     % update compatible trace element composition
%     advn_TRC = flxdiv(TRCs,fUs,fWs,h,ADVN,'adv') + flxdiv(TRCf,fUf,fWf,h,ADVN,'adv'); % advection
%
%     cmpt_TRC = (TRCs(ic,ic)-TRCf(ic,ic)).*Div_fV;                          % advection/compaction
%
%     diff_TRC = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(TRC(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(TRC(ic,:),1,2) ,1,2)./h^2;
%
%     RTRCin   = (TRCin-TRC)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate
%
%     RTRC(ic,ic) = - advn_TRC - cmpt_TRC + diff_TRC + RTRCin(ic,ic);    % total rate of change
%
%     res_TRC = (TRC-TRCo)./dt - (theta.*RTRC + (1-theta).*RTRCo);           % residual composition evolution equation
%
%     TRC     = TRC - res_TRC.*dt/4;                                         % update composition solution
%     TRC     = max(1e-16,min(TRCs,TRC));                                    % enforce min/max bounds
%
%     TRC([1 end],:) = TRC(ibz,:);                                           % apply boundary conditions
%     TRC(:,[1 end]) = TRC(:,ibx);
%     res_TRC([1 end],:) = 0;
%     res_TRC(:,[1 end]) = 0;
%
%
%     % *****  RADIOGENIC ISOTOPES  *****************************************
%
%     % update parent isotope composition
%     advn_IRP = flxdiv(IRPs,fUs,fWs,h,ADVN,'adv') + flxdiv(IRPf,fUf,fWf,h,ADVN,'adv'); % advection
%
%     cmpt_IRP = (IRPs(ic,ic)-IRPf(ic,ic)).*Div_fV;                          % advection/compaction
%
%     diff_IRP = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(IRP(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(IRP(ic,:),1,2) ,1,2)./h^2;
%
%     RIRPin   = (IRPin.*exp(-time./DIRP.*log(2))-IRP)./(5*dt) .* exp(-(L-Z)./(2*h)); % base injection rate
%
%     RIRP(ic,ic) = - advn_IRP - cmpt_IRP + diff_IRP - DcyR_IRP(ic,ic) + RIRPin(ic,ic); % total rate of change;
%
%     res_IRP = (IRP-IRPo)./dt - (theta.*RIRP + (1-theta).*RIRPo);           % residual composition evolution equation
%
%     IRP     = IRP - res_IRP.*dt/4;                                         % update composition solution
%     IRP     = max(1e-16,min(1e3,IRP));                                     % enforce min/max bounds
%
%     IRP([1 end],:) = IRP(ibz,:);                                           % apply boundary conditions
%     IRP(:,[1 end]) = IRP(:,ibx);
%     res_IRP([1 end],:) = 0;
%     res_IRP(:,[1 end]) = 0;
%
%
%     % update daughter isotope composition
%     advn_IRD = flxdiv(IRDs,fUs,fWs,h,ADVN,'adv') + flxdiv(IRDf,fUf,fWf,h,ADVN,'adv'); % advection
%
%     cmpt_IRD = (IRDs(ic,ic)-IRDf(ic,ic)).*Div_fV;                          % advection/compaction
%
%     diff_IRD = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(IRD(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(IRD(ic,:),1,2) ,1,2)./h^2;
%
%     RIRDin   = (IRDin+IRPin.*(1-exp(-time./DIRP.*log(2)))-IRD)./(5*dt) .* exp(-(L-Z)./(2*h)); % base injection rate
%
%     RIRD(ic,ic) = - advn_IRD - cmpt_IRD + diff_IRD + DcyR_IRP(ic,ic) - DcyR_IRD(ic,ic) + RIRDin(ic,ic); % total rate of change;
%
%     res_IRD = (IRD-IRDo)./dt - (theta.*RIRD + (1-theta).*RIRDo);           % residual composition evolution equation
%
%     IRD     = IRD - res_IRD.*dt/4;                                         % update composition solution
%     IRD     = max(1e-16,min(1e3,IRD));                                     % enforce min/max bounds
%
%     IRD([1 end],:) = IRD(ibz,:);                                           % apply boundary conditions
%     IRD(:,[1 end]) = IRD(:,ibx);
%     res_IRD([1 end],:) = 0;
%     res_IRD(:,[1 end]) = 0;
%
%
%     % *****  STABLE ISOTOPES  *********************************************
%
%     % update solid stable isotope composition
%     advn_ISS = flxdiv(ISS,U,W,h,ADVN,'adv');                               % advection/compaction
%
%     diff_ISS = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(ISS(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(ISS(ic,:),1,2) ,1,2)./h^2;
%
%     RctR_ISS = -(ISR-ISS).*RctR_f./max(flim,1-f);                          % reactive transfer rate
%
%     RISSin   = (ISSin-ISS)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate
%
%     RISS(ic,ic) = - advn_ISS + diff_ISS + RctR_ISS(ic,ic) + RISSin(ic,ic);  % total rate of change
%
%     res_ISS = (ISS-ISSo)./dt - (theta.*RISS + (1-theta).*RISSo);           % residual composition evolution equation
%
%     ISS     = ISS - res_ISS.*dt/4;                                         % update composition solution
%     ISS     = max(min(ISSin(:)),min(max(ISSin(:)),ISS));                   % enforce min/max bounds
%
%     ISS([1 end],:) = ISS(ibz,:);                                           % apply boundary conditions
%     ISS(:,[1 end]) = ISS(:,ibx);
%     res_ISS([1 end],:) = 0;
%     res_ISS(:,[1 end]) = 0;
%
%
%     % update fluid stable isotope composition
%     advn_ISF = flxdiv(ISF,Uf,Wf,h,ADVN,'adv');                             % advection/compaction
%
%     diff_ISF = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(ISF(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(ISF(ic,:),1,2) ,1,2)./h^2;
%
%     RctR_ISF = (ISR-ISF).*RctR_f./max(flim,f);                             % reactive transfer rate
%
%     RISFin   = (ISFin-ISF)./(5*dt) .* exp(-(L-Z)./(2*h));                  % base injection rate
%
%     RISF(ic,ic) = - advn_ISF + diff_ISF + RctR_ISF(ic,ic) + RISFin(ic,ic); % total rate of change
%
%     res_ISF = (ISF-ISFo)./dt - (theta.*RISF + (1-theta).*RISFo);           % residual composition evolution equation
%
%     ISF     = ISF - res_ISF.*dt/4;                                         % update composition solution
%     ISF     = max(min(ISFin(:)),min(max(ISFin(:)),ISF));                   % enforce min/max bounds
%
%     ISF([1 end],:) = ISF(ibz,:);                                           % apply boundary conditions
%     ISF(:,[1 end]) = ISF(:,ibx);
%     res_ISF([1 end],:) = 0;
%     res_ISF(:,[1 end]) = 0;
%
%
%     %*****  DAMAGE EVOLUTION  *************************************************
%
%     advn_DMG = flxdiv(DMG,U,W,h,ADVN,'adv');                               % damage advection
%
%     diff_DMG = diff( (1./PeCeff(im,ic)+1./PeCeff(ip,ic))/2 .* diff(DMG(:,ic),1,1) ,1,1)./h^2 ...
%              + diff( (1./PeCeff(ic,im)+1./PeCeff(ic,ip))/2 .* diff(DMG(ic,:),1,2) ,1,2)./h^2;
%
%     RDMG(ic,ic) = -advn_DMG + diff_DMG + epsDMG(ic,ic) - RHEAL.*DMG(ic,ic);% total rate of change
%
%     res_DMG = (DMG-DMGo)./dt - (theta.*RDMG + RDMGo);                      % residual damage evolution equation
%
%     DMG = DMG - res_DMG.*dt/4;                                               % update composition solution
%     DMG = max(1e-16,DMG);                                                  % enforce min bound
%
%     DMG([1 end],:) = DMG(ibz,:);                                           % apply boundary conditions
%     DMG(:,[1 end]) = DMG(:,ibx);
%     res_DMG([1 end],:) = 0;
%     res_DMG(:,[1 end]) = 0;
%
% end
