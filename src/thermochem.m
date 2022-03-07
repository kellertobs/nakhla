%% update time step
dtk = min((h/2)^2./max([kT(:)./rhoCp(:);kc]));                             % diffusive time step size
dta = min(min(h/2/max(abs([UBG(:);WBG(:);Um(:);Wm(:);Uf(:);Wf(:);Ux(:);Wx(:)]+1e-16)))); % advective time step size
dt  = min([2*dto,dtmax,CFL*min(dtk,dta)]);                                 % physical time step size


%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; xi = x;

% update temperature
advn_H = advection(rho.*m.*(Cpm + 0  ).*T,Um,Wm,h,ADVN,'flx') ...
       + advection(rho.*x.*(Cpx + Dsx).*T,Ux,Wx,h,ADVN,'flx') ...
       + advection(rho.*f.*(Cpf + Dsf).*T,Uf,Wf,h,ADVN,'flx');
                           
qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                     % heat diffusion z-flux
qTx    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddx(T,h);                     % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...                     % heat diffusion
                           - ddx(qTx(2:end-1,:),h));
    
bndH = zeros(size(T));
if ~isnan(Twall); bndH = bndH + rhoCp.*(Twall-T)./tau_T .* bndshape; end   % impose boundary layer

dHdt = - advn_H + diff_T + bndH;                                           % total rate of change
    
if step>0; H = Ho + (dHdt + dHdto)/2.*dt; end                              % explicit update of enthalpy
H([1 end],:) = H([2 end-1],:);                                             % apply boundary conditions
H(:,[1 end]) = H(:,[2 end-1]);    
    
% update major component
advn_C = advection(rho.*m.*cm,Um,Wm,h,ADVN,'flx') ...
       + advection(rho.*x.*cx,Ux,Wx,h,ADVN,'flx');
                          
qcz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(c,h);  % major component diffusion z-flux
qcx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(c,h);  % major component diffusion x-flux
diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                      % major component diffusion
                          - ddx(qcx(2:end-1,:),h);
    
bndC = zeros(size(c));
if ~isnan(cwall); bndC = bndC + rho.*(cwall-c)./tau_a .* bndshape; end     % impose boundary layer

dCdt = - advn_C + diff_c + bndC;                                           % total rate of change
    
if step>0; C = Co + (dCdt + dCdto)/2.*dt; end                              % explicit update of major component density
C = min(rho.*cphs1-TINY,max(rho.*cphs0+TINY, C ));
C([1 end],:) = C([2 end-1],:);                                             % apply boundary conditions
C(:,[1 end]) = C(:,[2 end-1]);  
    
% update volatile component
bndV = zeros(size(v));
if any(v(:)>1e-6)
    advn_V = advection(rho.*m.*vm,Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*f.*vf,Uf,Wf,h,ADVN,'flx');
    
    qvz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(v,h);  % volatile component diffusion z-flux
    qvx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(v,h);  % volatile component diffusion x-flux
    diff_v(2:end-1,2:end-1) = - ddz(qvz(:,2:end-1),h) ...                      % volatile component diffusion
        - ddx(qvx(2:end-1,:),h);
    
    if ~isnan(vwall); bndV = bndV + rho.*(vwall-v)./tau_a .* bndshape; end     % impose boundary layer
    
    dVdt = - advn_V + diff_v + bndV;                                           % total rate of change
    
    if step>0; V = Vo + (dVdt + dVdto)/2.*dt; end                              % explicit update of volatile component density
    V = min(rho.*1-TINY,max(rho.*0+TINY, V ));
    V([1 end],:) = V([2 end-1],:);                                             % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);
end

% convert enthalpy and component densities to temperature and concentrations
if step>0
    T = H./(rhoCp + rhoDs);
    c = C./rho;
    v = V./rho;
else
    H = (rhoCp + rhoDs).*T;
    C = rho.*(m.*cm + x.*cx);
    V = rho.*(m.*vm + f.*vf);
end


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
if react
    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
end

% update crystal fraction
if diseq || ~react
    
    if react
        Gx = alpha.*Gx + (1-alpha).*((xq-x).*rho./max(4.*dt,tau_r));
        Gx = max(-0.001*rhox/dt,min(0.001*rho/dt,Gx));  % limit to 0.1 vol% phase change per time step
    end
    
    advn_x = advection(rho.*x,Ux,Wx,h,ADVN,'flx');                         % get advection term
    
    qxz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2 .* ddz(x,h);             % crystal fraction diffusion z-flux
    qxx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2 .* ddx(x,h);             % crystal fraction diffusion x-flux
    diff_x(2:end-1,2:end-1) = - ddz(qxz(:,2:end-1),h) ...                  % crystal fraction diffusion
                              - ddx(qxx(2:end-1,:),h);
    
    dxdt   = - advn_x + diff_x + Gx;                                       % total rate of change
    
    if step>0; x = (rhoo.*xo + (dxdt + dxdto)/2.*dt)./rho; end             % explicit update of crystal fraction
    x = min(1-f-TINY,max(TINY,x));                                         % enforce [0,1] limit
    x([1 end],:) = x([2 end-1],:);                                         % apply boundary conditions
    x(:,[1 end]) = x(:,[2 end-1]);
    
else
    
    x  =  alpha.*x + (1-alpha).*xq;
    if step>0; Gx = (rho.*x-rhoo.*xo)./dt + advection(rho.*x,Ux,Wx,h,ADVN,'flx'); end  % reconstruct crystallisation rate
    
end

% update bubble fraction
if (diseq && any(v(:)>1e-6)) || ~react
    
    if react
        Gf = alpha.*Gf + (1-alpha).*((fq-f).*rho./max(4.*dt,tau_r));
        Gf = max(-0.001*rhof/dt,min(0.001*rho/dt,Gf));  % limit to 0.1 vol% phase change per time step
    end
    
    advn_f = advection(rho.*f,Uf,Wf,h,ADVN,'flx');                         % get advection term
    
    qfz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2 .* ddz(f,h);             % bubble fraction diffusion z-flux
    qfx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2 .* ddx(f,h);             % bubble fraction diffusion x-flux
    diff_f(2:end-1,2:end-1) = - ddz(qfz(:,2:end-1),h) ...                  % bubble fraction diffusion
                              - ddx(qfx(2:end-1,:),h);
                          
    dfdt   = - advn_f + diff_f + Gf;                                       % total rate of change
    
    if step>0; f = (rhoo.*fo + (dfdt + dfdto)/2.*dt)./rho; end             % explicit update of bubble fraction
    f = min(1-x-TINY,max(TINY,f));                                         % enforce [0,1-x] limit
    f([1 end],:) = f([2 end-1],:);                                         % apply boundary conditions
    f(:,[1 end]) = f(:,[2 end-1]);
    
else
    
    f  =  alpha.*f + (1-alpha).*fq;
    if step>0; Gf = (rho.*f-rhoo.*fo)./dt + advection(rho.*f,Uf,Wf,h,ADVN,'flx'); end  % reconstruct exsolution rate
    
end

% update melt fraction
m = min(1-x-TINY,max(TINY,1-f-x));

% update phase compositions
if react && step>0

    % major component
    Kc = cxq./cmq;
    cm = c./(m + x.*Kc);
    cx = c./(m./Kc + x);
    
    % volatile component
    if any(v(:)>1e-6)
        Kf = vfq./vmq;
        vm = v./(m + f.*Kf);
        vf = v./(m./Kf + f);
    end

end

% get residual of thermochemical equations from iterative update
resnorm_TC = norm(T - Ti,2)./(norm(T,2)+TINY) ...
           + norm(x - xi,2)./(norm(x,2)+TINY);


%% ***** TRACE & ISOTOPE GEOCHEMISTRY  ************************************

% *****  Incompatible Trace Element  **************************************

% update incompatible trace element phase compositions
if react
    itm = it./(m + x.*KIT);
    itx = it./(m./KIT + x);
end

% update incompatible trace element composition
advn_IT = advection(rho.*m.*itm,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*itx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(it,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(it,h);
diff_it(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(itwall); assim = assim + (itwall-it).*rho./tau_a .* bndshape; end % impose wall assimilation
    
dITdt = - advn_IT + diff_it + assim;                                       % total rate of change

if step>0; IT = ITo + (dITdt + dITdto)/2.*dt; end                          % explicit update
IT = max(0+TINY, IT );
IT([1 end],:) = IT([2 end-1],:);                                           % boundary conditions
IT(:,[1 end]) = IT(:,[2 end-1]);


% *****  COMPATIBLE TRACE ELEMENT  ****************************************

% update compatible trace element phase compositions
if react
    ctm = ct./(m + x.*KCT);
    ctx = ct./(m./KCT + x);
end

% update compatible trace element composition
advn_CT = advection(rho.*m.*ctm,Um,Wm,h,ADVN,'flx') ...
        + advection(rho.*x.*ctx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(ct,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(ct,h);
diff_ct(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(ct));
if ~isnan(ctwall); assim = assim + (ctwall-ct).*rho./tau_a .* bndshape; end % impose wall assimilation

dCTdt = - advn_CT + diff_ct + assim;                                       % total rate of change

if step>0; CT = CTo + (dCTdt + dCTdto)/2.*dt; end                          % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);


% *****  STABLE ISOTOPE RATIO  ********************************************

% reactive transfer of stable isotope ratio
trns_si = Gx.*(sim.*double(Gx<0) + six.*double(Gx>=0));

% update stable isotope ratio in melt
advn_sim = advection(rho.*m.*sim,Um,Wm,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(sim,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(sim,h);
diff_sim(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-sim).*rho.*m./tau_a .* bndshape; end % impose wall assimilation

dSImdt = - advn_sim + diff_sim + assim - trns_si;                          % total rate of change

if step>0; SIm = SImo + (dSImdt + dSImdto)/2.*dt; end                      % explicit update
SIm = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho.*m,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho.*m,SIm));
SIm(m < TINY)  = 0;
SIm([1 end],:) = SIm([2 end-1],:);                                         % boundary conditions
SIm(:,[1 end]) = SIm(:,[2 end-1]);


% update stable isotope ratio in xtals
advn_six = advection(rho.*x.*six,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(six,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(six,h);
diff_six(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion when melt present
                            - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-six).*rho.*x./tau_a .* bndshape; end % impose wall assimilation

dSIxdt = - advn_six + diff_six + assim + trns_si;                          % total rate of change

if step>0; SIx = SIxo + (dSIxdt + dSIxdto)/2.*dt; end                      % explicit update
SIx = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho.*x,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho.*x,SIx));
SIx(x < TINY)  = 0;
SIx([1 end],:) = SIx([2 end-1],:);                                         % boundary conditions
SIx(:,[1 end]) = SIx(:,[2 end-1]);

SI = SIm + SIx;
SI = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho,SI));


% *****  RADIOGENIC ISOTOPES  *********************************************

% decay rate of radiogenic isotope
dcy_rip = rho.*rip./HLRIP.*log(2);
dcy_rid = rho.*rid./HLRID.*log(2);

% update radiogenic parent isotope phase compositions
if react
    ripm = rip./(m + x.*KRIP);
    ripx = rip./(m./KRIP + x);
end

% update radiogenic parent isotope composition
advn_RIP = advection(rho.*m.*ripm,Um,Wm,h,ADVN,'flx') ...
         + advection(rho.*x.*ripx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(rip,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(rip,h);
diff_rip(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(riwall); assim = assim + (riwall-rip).*rho./tau_a .* bndshape; end % impose wall assimilation
                                       % secular equilibrium!
dRIPdt = - advn_RIP + diff_rip + assim - dcy_rip + dcy_rip;                % total rate of change
                                       
if step>0; RIP = RIPo + (dRIPdt + dRIPdto)/2.*dt; end                      % explicit update
RIP = max(0+TINY, RIP );
RIP([1 end],:) = RIP([2 end-1],:);                                         % boundary conditions
RIP(:,[1 end]) = RIP(:,[2 end-1]);


% update radiogenic daughter isotope phase compositions
if react
    ridm = rid./(m + x.*KRID);
    ridx = rid./(m./KRID + x);
end

% update radiogenic daughter isotope composition
advn_RID = advection(rho.*m.*ridm,Um,Wm,h,ADVN,'flx') ...
         + advection(rho.*x.*ridx,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(rho(1:end-1,:)+rho(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(rid,h);
qx   = - kc.*(rho(:,1:end-1)+rho(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(rid,h);
diff_rid(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(riwall); assim = assim + (riwall.*HLRID./HLRIP-rid).*rho./tau_a .* bndshape; end % impose wall assimilation
    
dRIDdt = - advn_RID + diff_rid + assim - dcy_rid + dcy_rip;                % total rate of change

if step>0; RID = RIDo + (dRIDdt + dRIDdto)/2.*dt; end                      % explicit update
RID = max(0+TINY, RID );
RID([1 end],:) = RID([2 end-1],:);                                         % boundary conditions
RID(:,[1 end]) = RID(:,[2 end-1]);

if step>0
    it  = IT./rho;
    ct  = CT./rho;
    sim = SIm./rho./max(TINY,m);
    six = SIx./rho./max(TINY,x);
    si  = SI./rho./max(TINY,1-f);
    rip = RIP./rho;
    rid = RID./rho;
else
    IT  = rho.*(m.*itm + x.*itx);
    CT  = rho.*(m.*ctm + x.*ctx);
    SIm = rho.* m.*sim;
    SIx = rho.* x.*six;
    SI  = SIm + SIx;
    RIP = rho.*(m.*ripm + x.*ripx);
    RID = rho.*(m.*ridm + x.*ridx);
end