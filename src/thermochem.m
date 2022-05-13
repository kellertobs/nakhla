%% *****  THERMO-CHEMICAL EVOLUTION  **************************************

% store previous iteration
Ti = T; xi = x; fi = f;

% update temperature
advn_H = advection(rho.*m.*(Cpm.*T + 0 ),Um,Wm,h,ADVN,'flx') ...
       + advection(rho.*x.*(Cpx.*T + Lx),Ux,Wx,h,ADVN,'flx') ...
       + advection(rho.*f.*(Cpf.*T + Lf),Uf,Wf,h,ADVN,'flx');
                           
qTz    = - (kT(1:end-1,:)+kT(2:end,:))./2 .* ddz(T,h);                     % heat diffusion z-flux
qTx    = - (kT(:,1:end-1)+kT(:,2:end))./2 .* ddx(T,h);                     % heat diffusion x-flux
diff_T(2:end-1,2:end-1) = (- ddz(qTz(:,2:end-1),h) ...                     % heat diffusion
                           - ddx(qTx(2:end-1,:),h));
    
bndH = zeros(size(T));
if ~isnan(Twall); bndH = bndH + rhoCp.*(Twall-T)./tau_T .* bndshape; end   % impose boundary layer

dHdt = - advn_H + diff_T + bndH;                                           % total rate of change
    
H = Ho + (THETA.*dHdt + (1-THETA).*dHdto).*dt;            % explicit update of enthalpy
H([1 end],:) = H([2 end-1],:);                                             % apply boundary conditions
H(:,[1 end]) = H(:,[2 end-1]);    
    
% update major component
advn_C = advection(rho.*m.*cm,Um,Wm,h,ADVN,'flx') ...
       + advection(rho.*x.*cx,Ux,Wx,h,ADVN,'flx');
                          
qcz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(c,h);  % major component diffusion z-flux
qcx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(c,h);  % major component diffusion x-flux
diff_c(2:end-1,2:end-1) = - ddz(qcz(:,2:end-1),h) ...                      % major component diffusion
                          - ddx(qcx(2:end-1,:),h);
    
bndC = zeros(size(c));
if ~isnan(cwall); bndC = bndC + rho.*(cwall-c)./tau_a .* bndshape; end     % impose boundary layer

dCdt = - advn_C + diff_c + bndC;                                           % total rate of change
    
C = Co + (THETA.*dCdt + (1-THETA).*dCdto).*dt;                             % explicit update of major component density
C = min(rho.*cphs1-TINY,max(rho.*cphs0+TINY, C ));
C([1 end],:) = C([2 end-1],:);                                             % apply boundary conditions
C(:,[1 end]) = C(:,[2 end-1]);  
    
% update volatile component
bndV = zeros(size(v));
if any([v0;v1;vwall;v(:)]>1e-6)
    advn_V = advection(rho.*m.*vm,Um,Wm,h,ADVN,'flx') ...
           + advection(rho.*f.*vf,Uf,Wf,h,ADVN,'flx');
    
    qvz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(v,h);  % volatile component diffusion z-flux
    qvx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(v,h);  % volatile component diffusion x-flux
    diff_v(2:end-1,2:end-1) = - ddz(qvz(:,2:end-1),h) ...                      % volatile component diffusion
        - ddx(qvx(2:end-1,:),h);
    
    if ~isnan(vwall); bndV = bndV + rho.*(vwall-v)./tau_a .* bndshape; end     % impose boundary layer
    
    dVdt = - advn_V + diff_v + bndV;                                           % total rate of change
    
    V = Vo + (THETA.*dVdt + (1-THETA).*dVdto).*dt;                              % explicit update of volatile component density
    V = min(rho.*1-TINY,max(rho.*0+TINY, V ));
    V([1 end],:) = V([2 end-1],:);                                             % apply boundary conditions
    V(:,[1 end]) = V(:,[2 end-1]);
end

% convert enthalpy and component densities to temperature and concentrations
T = (H - rhoLh)./rhoCp;
c = C./rho;
v = V./rho;


%% *****  UPDATE PHASE PROPORTIONS  ***************************************

% update local phase equilibrium
if react && ~mod(iter-1,1)
    [xq,cxq,cmq,fq,vfq,vmq] = equilibrium(x,f,T,c,v,Pt,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg,beta);
end

% update crystal fraction
if diseq || ~react
    
    if react
        Gxi = (xq-x).*rho./max(5.*dt,tau_r);
        for i = 1:round(delta)
            Gxi(2:end-1,2:end-1) = Gxi(2:end-1,2:end-1) + diff(Gxi(:,2:end-1),2,1)./8 + diff(Gxi(2:end-1,:),2,2)./8;
            Gxi([1 end],:) = Gxi([2 end-1],:);
            Gxi(:,[1 end]) = Gxi(:,[2 end-1]);
        end
        Gx = ALPHA.*Gx + (1-ALPHA).*Gxi;
    end
    
    advn_x = advection(rho.*x,Ux,Wx,h,ADVN,'flx');                         % get advection term
    
    dxdt   = - advn_x + Gx;                                                % total rate of change
    
    x = (rhoo.*xo + (THETA.*dxdt + (1-THETA).*dxdto).*dt)./rho;            % explicit update of crystal fraction
    x = min(1-f-TINY,max(TINY,x));                                         % enforce [0,1] limit
    x([1 end],:) = x([2 end-1],:);                                         % apply boundary conditions
    x(:,[1 end]) = x(:,[2 end-1]);
    
else
    
    x  =  ALPHA.*x + (1-ALPHA).*xq;
    Gx = (rho.*x-rhoo.*xo)./dt + advection(rho.*x,Ux,Wx,h,ADVN,'flx');     % reconstruct crystallisation rate
    
end

% update bubble fraction
if (diseq && any([v0;v1;vwall;v(:)]>1e-6)) || ~react
    
    if react
        Gf = max(-0.5,min(0.5,ALPHA.*Gf + (1-ALPHA).*((fq-f).*rho./max(5.*dt,tau_r))));
    end
    
    advn_f = advection(rho.*f,Uf,Wf,h,ADVN,'flx');                         % get advection term
                          
    dfdt   = - advn_f + Gf;                                                % total rate of change
    
    f = (rhoo.*fo + (THETA.*dfdt + (1-THETA).*dfdto).*dt)./rho;            % explicit update of bubble fraction
    f = min(1-TINY,max(TINY,f));                                           % enforce [0,1-x] limit
    f([1 end],:) = f([2 end-1],:);                                         % apply boundary conditions
    f(:,[1 end]) = f(:,[2 end-1]);
    
else
    
    f  =  ALPHA.*f + (1-ALPHA).*fq;
    Gf = (rho.*f-rhoo.*fo)./dt + advection(rho.*f,Uf,Wf,h,ADVN,'flx');     % reconstruct exsolution rate
    
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
    if any([v0;v1;vwall;v(:)]>1e-6)
        Kf = vfq./vmq;
        vm = v./(m + f.*Kf);
        vf = v./(m./Kf + f);
        vf(v<1e-6) = vfq(v<1e-6);
    end

end

% get residual of thermochemical equations from iterative update
resnorm_TC = norm(T - Ti,2)./(norm(T,2)+TINY) ...
           + norm((x - xi).*(x>1e-6),2)./(norm(x,2)+TINY) ...
           + norm((f - fi).*(f>1e-6),2)./(norm(f,2)+TINY);


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

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(it,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(it,h);
diff_it(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(itwall); assim = assim + (itwall-it).*rho./tau_a .* bndshape; end % impose wall assimilation
    
dITdt = - advn_IT + diff_it + assim;                                       % total rate of change

IT = ITo + (THETA.*dITdt + (1-THETA).*dITdto).*dt;          % explicit update
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

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(ct,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(ct,h);
diff_ct(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                      % diffusion in melt
                           - ddx(qx(2:end-1,:),h);

assim = zeros(size(ct));
if ~isnan(ctwall); assim = assim + (ctwall-ct).*rho./tau_a .* bndshape; end % impose wall assimilation

dCTdt = - advn_CT + diff_ct + assim;                                       % total rate of change

CT = CTo + (THETA.*dCTdt + (1-THETA).*dCTdto).*dt;          % explicit update
CT = max(0+TINY, CT );
CT([1 end],:) = CT([2 end-1],:);                                           % boundary conditions
CT(:,[1 end]) = CT(:,[2 end-1]);


% *****  STABLE ISOTOPE RATIO  ********************************************

% reactive transfer of stable isotope ratio
trns_si = Gx.*(sim.*double(Gx<0) + six.*double(Gx>0));

% update stable isotope ratio in melt
advn_sim = advection(rho.*m.*sim,Um,Wm,h,ADVN,'flx');

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(sim,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(sim,h);
diff_sim(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-sim).*rho.*m./tau_a .* bndshape; end % impose wall assimilation

dSImdt = - advn_sim + diff_sim + assim - trns_si;                          % total rate of change

SIm = SImo + (THETA.*dSImdt + (1-THETA).*dSImdto).*dt;      % explicit update
SIm = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho.*m,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho.*m,SIm));
SIm(m <= TINY) = 0;
SIm([1 end],:) = SIm([2 end-1],:);                                         % boundary conditions
SIm(:,[1 end]) = SIm(:,[2 end-1]);


% update stable isotope ratio in xtals
advn_six = advection(rho.*x.*six,Ux,Wx,h,ADVN,'flx');

qz   = - kc.*(x(1:end-1,:)+x(2:end,:))/2.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(six,h);
qx   = - kc.*(x(:,1:end-1)+x(:,2:end))/2.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(six,h);
diff_six(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion when melt present
                            - ddx(qx(2:end-1,:),h);
                       
assim = zeros(size(si));
if ~isnan(siwall); assim = assim + (siwall-six).*rho.*x./tau_a .* bndshape; end % impose wall assimilation

dSIxdt = - advn_six + diff_six + assim + trns_si;                          % total rate of change

SIx = SIxo + (THETA.*dSIxdt + (1-THETA).*dSIxdto).*dt;      % explicit update
SIx = min(max([max(0,siwall),si0-dsi,si1+dsi]).*rho.*x,max(min([min(0,siwall),si0-dsi,si1+dsi]).*rho.*x,SIx));
SIx(x <= TINY) = 0;
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

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(rip,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(rip,h);
diff_rip(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(riwall); assim = assim + (riwall-rip).*rho./tau_a .* bndshape; end % impose wall assimilation
                                       % secular equilibrium!
dRIPdt = - advn_RIP + diff_rip + assim - dcy_rip + dcy_rip;                % total rate of change
                                       
RIP = RIPo + (THETA.*dRIPdt + (1-THETA).*dRIPdto).*dt;                     % explicit update
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

qz   = - kc.*(m(1:end-1,:)+m(2:end,:))/2 .* ddz(rid,h);
qx   = - kc.*(m(:,1:end-1)+m(:,2:end))/2 .* ddx(rid,h);
diff_rid(2:end-1,2:end-1) = - ddz(qz(:,2:end-1),h) ...                     % diffusion in melt
                            - ddx(qx(2:end-1,:),h);

assim = zeros(size(it));
if ~isnan(riwall); assim = assim + (riwall.*HLRID./HLRIP-rid).*rho./tau_a .* bndshape; end % impose wall assimilation
    
dRIDdt = - advn_RID + diff_rid + assim - dcy_rid + dcy_rip;                % total rate of change

RID = RIDo + (THETA.*dRIDdt + (1-THETA).*dRIDdto).*dt;                     % explicit update
RID = max(0+TINY, RID );
RID([1 end],:) = RID([2 end-1],:);                                         % boundary conditions
RID(:,[1 end]) = RID(:,[2 end-1]);

it  = IT./rho;
ct  = CT./rho;
sim = SIm./rho./max(TINY,m);
six = SIx./rho./max(TINY,x);
si  = SI ./rho./max(TINY,1-f);
rip = RIP./rho;
rid = RID./rho;