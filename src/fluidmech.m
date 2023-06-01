tic;

if ~bnchm && step>0

%***  update mixture mass density
drhodt  = advn_rho;

% residual of mixture mass evolution
res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);

% volume source and background velocity passed to fluid-mechanics solver
VolSrc  = Div_V - res_rho./rho/2;  % correct volume source term by scaled residual

UBG     = - 0*mean(VolSrc,'all')./2 .* (L/2-XXu);
WBG     = - 2*mean(VolSrc,'all')./2 .* (D/2-ZZw);

end

%% 0-D run does not require fluidmech solve
if Nz==1 && Nx==1  
    W  = WBG; Wm = W;  Wx = W;  Wf = W;
    U  = UBG; Um = U;  Ux = U;  Uf = U;
    P  = zeros(Nz+2,Nx+2);
    Pt = Ptop.*ones(Nz,Nx);
    resnorm_VP = 0;
else


%% assemble coefficients for matrix velocity diagonal and right-hand side

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R


% assemble coefficients of z-stress divergence
    
% left boundary
ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,end-1);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,2);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% top boundary
ii  = MapW(1,:); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii)) + WBG(1,:);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapW(end,:); jj = ii;
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa  = zeros(size(ii)) + WBG(end,:);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 =  etaco(2:end-1,1:end-1);   EtaC2 =  etaco(2:end-1,2:end);
EtaP1 =  eta  (1:end-1,:      );   EtaP2 =  eta  (2:end,:      );

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa  = - a1.*rhofz(2:end-1,:)./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = - 1/2*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 1/2*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 1/2*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = -ddz(rho,h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; (1/2*EtaC1(:)-1/2*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/2*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/2*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; (1/2*EtaC2(:)-1/2*EtaP2(:))/h^2];  % U one to the bottom and right


% z-RHS vector
rr  = - (rhofz(2:end-1,:) - mean(rhofz(2:end-1,:),2)) .* g0 - (a2.*rhoWo(2:end-1,:)+a3.*rhoWoo(2:end-1,:))/dt;
if bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii  = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+top];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii  = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+bot];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% internal points
ii    = MapU(2:end-1,:);
EtaC1 = etaco(1:end-1,:    );  EtaC2 = etaco(2:end,:    );
EtaP1 = eta  (:,[end,1:end]);  EtaP2 = eta  (:,[1:end,1]);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = MapU(2:end-1,[end-1,1:end-1]); jj2 = MapU(2:end-1,[2:end,2]); jj3 = MapU(1:end-2,:); jj4 = MapU(3:end,:);

aa  = - a1.*rhofx./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = - 1/2*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 1/2*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 1/2*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % U one below

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = -ddx(rho(:,[end,1:end,1]),h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,1:end-1); jj2 = MapW(1:end-1,2:end); jj3 = MapW(2:end,1:end-1); jj4 = MapW(2:end,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; (1/2*EtaC1(:)-1/2*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/2*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/2*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; (1/2*EtaC2(:)-1/2*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
rr  = - (a2.*rhoUo+a3.*rhoUoo)/dt;
if bnchm; rr = rr + src_U_mms(2:end-1,:); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV  = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RV  = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for gradient operator

if ~exist('GG','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    
    % coefficients for z-gradient
    ii  = MapW(2:end-1,2:end-1);
    
    %         top              ||          bottom
    jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
    
    
    % coefficients for x-gradient
    ii  = MapU(2:end-1,:);
    
    %         left             ||           right
    jj1 = MapP(2:end-1,1:end-1); jj2 = MapP(2:end-1,2:end);
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
    
    
    % assemble coefficient matrix
    GG  = sparse(IIL,JJL,AAL,NW+NU,NP);
end


%% assemble coefficients for divergence operator

if ~exist('DD','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    %internal points
    ii  = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying velocities U, W
    %          left U          ||           right U       ||           top W           ||          bottom W
    jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/h];  % W one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/h];  % W one below

    % assemble coefficient matrix
    DD  = sparse(IIL,JJL,AAL,NP,NW+NU);
end


%% assemble coefficients for matrix pressure diagonal and right-hand side

if ~exist('KP','var') || bnchm
    IIL = [];       % equation indeces into A
    JJL = [];       % variable indeces into A
    AAL = [];       % coefficients for A
    
    % boundary points
    ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
    jj1 = ii;
    jj2 = [MapP(2,:).'; MapP(end-1,:).'];
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    ii  = [MapP(:,1    ); MapP(:,end)]; % left & right
    jj1 = ii;
    if bndmode < 4
        jj2 = [MapP(:,end-1); MapP(:,2    )];
    else
        jj2 = [MapP(:,2    ); MapP(:,end-1)];
    end
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    % internal points
    ii  = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying matrix pressure P
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre
    
    % assemble coefficient matrix
    KP  = sparse(IIL,JJL,AAL,NP,NP);
end

% RHS
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

ii  = MapP(2:end-1,2:end-1);

rr  = - VolSrc;
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RP  = sparse(IIR,ones(size(IIR)),AAR,NP,1);

% set P = 0 in fixed point
nzp = round(Nz/2);
nxp = round(Nx/2);
DD(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),MapP(nzp,nxp)) = 1;
RP(MapP(nzp,nxp),:) = 0;

if bnchm; RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp); end


%% assemble and scale global coefficient matrix and right-hand side vector

LL  = [ KV  -GG  ; ...
       -DD   KP ];

RR  = [RV; RP];

SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse(1./(SCL+1)));

LL  = SCL*LL*SCL;
RR  = SCL*RR;


%% Solve linear system of equations for vx, vz, P

SOL = SCL*(LL\RR);  % update solution

% map solution vector to 2D arrays
W   = full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2)); % matrix z-velocity
U   = full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1)); %U = U-mean(U(2:end-1,:      ),'all');     % matrix x-velocity
P   = full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2)); %P = P-mean(P(2:end-1,2:end-1),'all');     % matrix dynamic pressure

% magma velocity magnitude
Vel = sqrt(((W(1:end-1,2:end-1)+W(2:end,2:end-1))/2).^2 ...
         + ((U(2:end-1,1:end-1)+U(2:end-1,2:end))/2).^2);

if ~bnchm

    % phase segregation speeds
    wm(2:end-1,2:end-1) = ((rhom(1:end-1,:)+rhom(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_m(1:end-1,:).*Ksgr_m(2:end,:)).^0.5; % melt segregation speed
    wm([1,end],:) = min(1,1-[top;bot]).*wm([2,end-1],:);
    wm(:,[1 end]) = wm(:,[end-1 2]);

    wx(2:end-1,2:end-1) = ((rhox(1:end-1,:)+rhox(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5; % solid segregation speed
    wx([1,end],:) = min(1,1-[top;bot]).*wx([2,end-1],:);
    wx(:,[1 end]) = wx(:,[end-1 2]);

    wf(2:end-1,2:end-1) = ((rhof(1:end-1,:)+rhof(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_f(1:end-1,:).*Ksgr_f(2:end,:)).^0.5; % fluid segregation speed
    wf([1,end],:) = min(1,1-[top-fout;bot-fin]).*wf([2,end-1],:);
    wf(:,[1 end]) = wf(:,[end-1 2]);

    % phase diffusion fluxes and speeds
    [~,qxz,qxx] = diffus(chi,kx,h,[1,2],BCD);
    [~,qfz,qfx] = diffus(phi,kf,h,[1,2],BCD);
    qmz  = -qxz-qfz;  qmx = -qxx-qfx;

    wqf = qfz./max(TINY^0.5,(phi([1,1:end],[end,1:end,1])+phi([1:end,end],[end,1:end,1]))./2);
    uqf = qfx./max(TINY^0.5,(phi([1,1:end,end],[end,1:end])+phi([1,1:end,end],[1:end,1]))./2);

    wqx = qxz./max(TINY^0.5,(chi([1,1:end],[end,1:end,1])+chi([1:end,end],[end,1:end,1]))./2);
    uqx = qxx./max(TINY^0.5,(chi([1,1:end,end],[end,1:end])+chi([1,1:end,end],[1:end,1]))./2);

    wqm = qmz./max(TINY^0.5,(mu ([1,1:end],[end,1:end,1])+mu ([1:end,end],[end,1:end,1]))./2);
    uqm = qmx./max(TINY^0.5,(mu ([1,1:end,end],[end,1:end])+mu ([1,1:end,end],[1:end,1]))./2);

    % update phase velocities
    Wf   = W + wf + wqf;  % mvp z-velocity
    Uf   = U + 0. + uqf;  % mvp x-velocity
    Wx   = W + wx + wqx;  % xtl z-velocity
    Ux   = U + 0. + uqx;  % xtl x-velocity
    Wm   = W + wm + wqm;  % mlt z-velocity
    Um   = U + 0. + uqm;  % mlt x-velocity

    
    %% update time step
    dtk = (h/2)^2./max(kT(:)./rho(:)./cP)/2;                                    % diffusive time step size
    dta = CFL*h/2/max(abs([Um(:).*any(m(:)>10*TINY);Wm(:).*any(m(:)>10*TINY); ... % advective time step size
                           Ux(:).*any(x(:)>10*TINY);Wx(:).*any(x(:)>10*TINY); ...
                           Uf(:).*any(f(:)>10*TINY);Wf(:).*any(f(:)>10*TINY)]+TINY));   
    dt  = min([1.05*dto,min(dtk,dta),dtmax]);                                      % time step size
end

end

FMtime = FMtime + toc;