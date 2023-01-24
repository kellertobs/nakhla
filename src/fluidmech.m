% store previous iteration
SOLi = SOL; tic;

%% 0-D run does not require fluidmech solve
if Nz==3 && Nx==3  
    W  = WBG; Wm = W;  Wx = W;  Wf = W;
    U  = UBG; Um = U;  Ux = U;  Uf = U;
    P  = zeros(Nz,Nx);
    Pt = Ptop.*ones(Nz,Nx);
    resnorm_VP = 0;
else


%% assemble coefficients for matrix velocity diagonal and right-hand side

IIL  = [];       % equation indeces into L
JJL  = [];       % variable indeces into L
AAL  = [];       % coefficients for L
IIR  = [];       % equation indeces into R
AAR  = [];       % forcing entries for R


% assemble coefficients of z-stress divergence
    
% left boundary
ii = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
ii = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% top boundary
ii = MapW(1,2:end-1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + WBG(1,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii = MapW(end,2:end-1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + WBG(end,2:end-1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
 EtaC1 =  etaco(2:end-1,1:end-1);   EtaC2 =  etaco(2:end-1,2:end  );
 EtaP1 =  eta  (2:end-2,2:end-1);   EtaP2 =  eta  (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with the drunken sailor...
aa = -ddz(rho(2:end-1,2:end-1),h).*g0.*dt/2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
    
% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right


% z-RHS vector
rr = - (rhofz(2:end-1,2:end-1) - mean(rhofz(2:end-1,2:end-1),2)) .* g0;
if bnchm; rr = rr + src_W_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+top];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% bottom boundary
ii = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+bot];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% left side boundary
ii = MapU(2:end-1,1); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + UBG(2:end-1,1);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right side boundary
ii = MapU(2:end-1,end); jj = ii;
aa = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
aa = zeros(size(ii)) + UBG(2:end-1,end);
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];


% internal points
ii    = MapU(2:end-1,2:end-1);
 EtaC1 =  etaco(1:end-1,2:end-1);   EtaC2 =  etaco(2:end  ,2:end-1);
 EtaP1 =  eta  (2:end-1,2:end-2);   EtaP2 =  eta  (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; 2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; 2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; 1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; 1/2*EtaC2(:)/h^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii));
if bnchm; rr = rr + src_U_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)];  AAR = [AAR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV = sparse(IIL,JJL,AAL,NW+NU,NW+NU);
RV = sparse(IIR,ones(size(IIR)),AAR);


%% assemble coefficients for gradient operator

if ~exist('GG','var') || bnchm
    IIL  = [];       % equation indeces into A
    JJL  = [];       % variable indeces into A
    AAL  = [];       % coefficients for A
    
    
    % coefficients for z-gradient
    ii = MapW(2:end-1,2:end-1);
    
    %         top              ||          bottom
    jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the top
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the bottom
    
    
    % coefficients for x-gradient
    ii = MapU(2:end-1,2:end-1);
    
    %         left             ||           right
    jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];     % one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];     % one to the right
    
    
    % assemble coefficient matrix
    GG = sparse(IIL,JJL,AAL,NW+NU,NP);
end


%% assemble coefficients for divergence operator

if ~exist('DD','var') || bnchm
    IIL  = [];       % equation indeces into A
    JJL  = [];       % variable indeces into A
    AAL  = [];       % coefficients for A
    
    %internal points
    ii = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying velocities U, W
    %          left U          ||           right U       ||           top W           ||          bottom W
    jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)-1/h];  % U one to the left
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+1/h];  % U one to the right
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; aa(:)-1/h];  % W one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; aa(:)+1/h];  % W one below

    % assemble coefficient matrix
    DD = sparse(IIL,JJL,AAL,NP,NW+NU);
end


%% assemble coefficients for matrix pressure diagonal and right-hand side

if ~exist('KP','var') || bnchm
    IIL  = [];       % equation indeces into A
    JJL  = [];       % variable indeces into A
    AAL  = [];       % coefficients for A
    
    % boundary points
    ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
    jj1 = ii;
    jj2 = [MapP(2,:).'; MapP(end-1,:).'];
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    ii  = [MapP(:,1); MapP(:,end  )]; % left & right
    jj1 = ii;
    jj2 = [MapP(:,2); MapP(:,end-1)];
    
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    % internal points
    ii = MapP(2:end-1,2:end-1);
    
    % coefficients multiplying matrix pressure P
    aa = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre
    
    % assemble coefficient matrix
    KP = sparse(IIL,JJL,AAL,NP,NP);
end

% RHS
IIR  = [];       % equation indeces into R
AAR  = [];       % forcing entries for R

ii = MapP(2:end-1,2:end-1);

rr = - VolSrc;
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RP = sparse(IIR,ones(size(IIR)),AAR,NP,1);

nzp = round((Nz-2)/2)+1;
nxp = round((Nx-2)/2)+1;
DD(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),MapP(nzp,nxp)) = 1;
RP(MapP(nzp,nxp),:) = 0;
if bnchm; RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp); end


%% assemble and scale global coefficient matrix and right-hand side vector

LL  = [ KV  -GG  ; ...
       -DD   KP ];

RR  = [RV; RP];

SCL = sqrt(abs(diag(LL)));
SCL = diag(sparse(1./(SCL+1)));

LL  = SCL*LL*SCL;
RR  = SCL*RR;


%% Solve linear system of equations for vx, vz, P

SOL = SCL*(LL\RR);  % update solution

% map solution vector to 2D arrays
W  = full(reshape(SOL(MapW(:))        ,(Nz-1), Nx   ));                    % matrix z-velocity
U  = full(reshape(SOL(MapU(:))        , Nz   ,(Nx-1)));                    % matrix x-velocity
P  = full(reshape(SOL(MapP(:)+(NW+NU)), Nz   , Nx   ));                    % matrix dynamic pressure

% magma velocity magnitude
Vel  = sqrt(((W([1,1:end],:)+W([1:end,end],:))/2).^2 ...
          + ((U(:,[1,1:end])+U(:,[1:end,end]))/2).^2);

if ~bnchm
    Pt( 2:end, :) = repmat(cumsum(mean(rhofz,2).*g0.*h),1,Nx);
    Pt            = Pt - Pt(2,:)/2 + Ptop;
    Pt([1 end],:) = Pt([2 end-1],:);
    Pt(:,[1 end]) = Pt(:,[2 end-1]);
    if Nz<=10; Pt = Ptop.*ones(size(Tp)); end

    % get residual of fluid mechanics equations from iterative update
    resSOL = SOL - SOLi;
    resW  = full(reshape(resSOL(MapW(:))        ,(Nz-1), Nx   ));          % matrix z-velocity
    resU  = full(reshape(resSOL(MapU(:))        , Nz   ,(Nx-1)));          % matrix x-velocity
    resP  = full(reshape(resSOL(MapP(:)+(NW+NU)), Nz   , Nx   ));          % matrix dynamic pressure
%     resnorm_VP = norm(resW(:  ,inx),'fro')./(norm(W(:  ,inx),'fro')+TINY) ...
%                + norm(resU(inz,:  ),'fro')./(norm(U(inz,:  ),'fro')+TINY) ...
%                + norm(resP(inz,inx),'fro')./(norm(P(inz,inx),'fro')+TINY);
    resnorm_VP = 0;

    qxz    = - (kx(1:end-1,:)+kx(2:end,:))./2 .* ddz(chi,h);  % z-flux
    qxx    = - (kx(:,1:end-1)+kx(:,2:end))./2 .* ddx(chi,h);  % x-flux

    qfz    = - (kf(1:end-1,:)+kf(2:end,:))./2 .* ddz(phi,h);  % z-flux
    qfx    = - (kf(:,1:end-1)+kf(:,2:end))./2 .* ddx(phi,h);  % x-flux

    qmz    = -qxz-qfz;
    qmx    = -qxx-qfx;

    % update phase velocities
    Wf   = W + wf + qfz./max(1e-6,(phi(1:end-1,:)+phi(2:end,:))./2);               % mvp z-velocity
    Uf   = U + 0. + qfx./max(1e-6,(phi(:,1:end-1)+phi(:,2:end))./2);               % mvp x-velocity
    Wx   = W + wx + qxz./max(1e-6,(chi(1:end-1,:)+chi(2:end,:))./2);               % xtl z-velocity
    Ux   = U + 0. + qxx./max(1e-6,(chi(:,1:end-1)+chi(:,2:end))./2);               % xtl x-velocity
    Wm   = W + wm + qmz./max(1e-6,(mu (1:end-1,:)+mu (2:end,:))./2);               % mlt z-velocity
    Um   = U + 0. + qmx./max(1e-6,(mu (:,1:end-1)+mu (:,2:end))./2);               % mlt x-velocity

    % update mixture volume flux
    Wbar = (mu (1:end-1,:)+mu (2:end,:))/2 .* Wm ...
         + (chi(1:end-1,:)+chi(2:end,:))/2 .* Wx ...
         + (phi(1:end-1,:)+phi(2:end,:))/2 .* Wf;
    Ubar = (mu (:,1:end-1)+mu (:,2:end))/2 .* Um ...
         + (chi(:,1:end-1)+chi(:,2:end))/2 .* Ux ...
         + (phi(:,1:end-1)+phi(:,2:end))/2 .* Uf;

    % mixture volume flux magnitude
    Vbar = sqrt(((Wbar([1,1:end],:)+Wbar([1:end,end],:))/2).^2 ...
              + ((Ubar(:,[1,1:end])+Ubar(:,[1:end,end]))/2).^2);

    
    %% update time step
    dtk = (h/2)^2./max([kT(:)./rho(:)./cP;kx(:);kf(:)])/2;                 % diffusive time step size
    dta = CFL*h/2/max(abs([Ux(:);Wx(:);Uf(:);Wf(:);Um(:);Wm(:)]+1e-16));   % advective time step size
    dt  = min([2*dto,dtmax,min(dtk,dta)]);                                 % physical time step size
end

end

FMtime = FMtime + toc;
