tic;

if ~bnchm && step>0 && ~restart

%***  update mixture mass density
drhodt  = advn_rho;% + (RHO-rho)/dt;

% residual of mixture mass evolution
res_rho = (a1*rho-a2*rhoo-a3*rhooo)/dt - (b1*drhodt + b2*drhodto + b3*drhodtoo);

% volume source and background velocity passed to fluid-mechanics solver
upd_rho = - res_rho./b1./rho; % + beta*upd_rho;
VolSrc  = VolSrc + upd_rho;  % correct volume source term by scaled residual

UBG     = - 0*mean(VolSrc,'all')./2 .* (L/2-XXu);
WBG     = - 2*mean(VolSrc,'all')./2 .* (D/2-ZZw);

dPchmbdt  = mod_wall*mean(VolSrc,'all') - mod_wall/eta_wall*Pchmb;
res_Pchmb = (a1*Pchmb-a2*Pchmbo-a3*Pchmboo)/dt - (b1*dPchmbdt + b2*dPchmbdto + b3*dPchmbdtoo);

upd_Pchmb = - alpha*res_Pchmb*dt/a1/3 + beta*upd_Pchmb;
Pchmb     = Pchmb + upd_Pchmb;

end

%% 0-D run does not require fluidmech solve
if Nz==1 && Nx==1  
    W  = WBG; Wm = W;  Wx = W;  Wf = W;
    U  = UBG; Um = U;  Ux = U;  Uf = U;
    P  = zeros(Nz+2,Nx+2);
    resnorm_VP = 0;
else

if Nx==1
    % update 1D velocity
    W(:,2) = -flipud(cumsum(flipud([VolSrc*h;-WBG(end)])));

    % update 1D pressure
    Div_tz = ddz(tzz(icz,:),h);           % z-stress divergence
    PrsSrc = - Div_tz;
    P(:,2) =  flipud(cumsum(flipud([PrsSrc*h;0])));
    P(:,2) = P(:,2) - P(Nz/2,2);

    W(:,1) = W(:,2); W(:,end) = W(:,2);
    P(:,1) = P(:,2); P(:,end) = P(:,2);
else


%% assemble coefficients for matrix velocity diagonal and right-hand side

IIL = [];       % equation indeces into L
JJL = [];       % variable indeces into L
AAL = [];       % coefficients for L
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R


% assemble coefficients of z-stress divergence
    
% left boundary
if periodic
    ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,end-1);
else
    ii  = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
end
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

% right boundary
if periodic
    ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,2);
else
    ii  = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
end
aa  = zeros(size(ii));
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)+sds];
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

aa  = a1.*rhofz(2:end-1,:)./dt;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % W on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % W one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % W one below
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % W one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % W one to the right

% what shall we do with the drunken sailor...
if ~bnchm
    aa  = ddz(rho,h).*g0.*dt;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % U one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % U one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % U one to the bottom and right


% z-RHS vector
advn_mz = advect(rhofz(2:end-1,:).*W(2:end-1,2:end-1),(U(2:end-2,:)+U(3:end-1,:))/2,(W(1:end-1,2:end-1)+W(2:end,2:end-1))/2,h,{ADVN,''},[1,2],BCA);
rr  = + (rhofz(2:end-1,:) - mean(rhofz(2:end-1,:),2)) .* g0 ...
      + (a2.*rhoWo(2:end-1,:)+a3.*rhoWoo(2:end-1,:))/dt ...
      - advn_mz;
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

if ~periodic
    % left boundary
    ii  = MapU(:,1); jj = ii;
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
    aa  = zeros(size(ii));
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];

    % right boundary
    ii  = MapU(:,end); jj = ii;
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj(:)];   AAL = [AAL; aa(:)+1];
    aa  = zeros(size(ii));
    IIR = [IIR; ii(:)]; AAR = [AAR; aa(:)];
end

% internal points
if periodic
    ii    = MapU(2:end-1,:);
    EtaC1 = etaco(1:end-1,:    );  EtaC2 = etaco(2:end,:    );
    EtaP1 = eta  (:,icx(1:end-1));  EtaP2 = eta  (:,icx(2:end));
else
    ii    = MapU(2:end-1,2:end-1);
    EtaC1 = etaco(1:end-1,2:end-1);  EtaC2 = etaco(2:end,2:end-1);
    EtaP1 = eta  (:      ,1:end-1);  EtaP2 = eta  (:      ,2:end);
end

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
if periodic
    jj1 = MapU(2:end-1,ifx(1:end-2)); jj2 = MapU(2:end-1,ifx(3:end)); jj3 = MapU(1:end-2,ifx(2:end-1)); jj4 = MapU(3:end,ifx(2:end-1));
else
    jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);
end
if periodic
    aa  = (a1+gamma).*rhofx./dt;
else
    aa  = a1.*rhofx(:,2:end-1)./dt;
end
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % inertial term

aa  = 2/3*(EtaP1+EtaP2)/h^2 + 1/2*(EtaC1+EtaC2)/h^2;
IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)           ];      % U on stencil centre
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-2/3*EtaP1(:)/h^2];      % U one to the left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;-2/3*EtaP2(:)/h^2];      % U one to the right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;-1/2*EtaC1(:)/h^2];      % U one above
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-1/2*EtaC2(:)/h^2];      % U one below

% what shall we do with the drunken sailor...
if ~bnchm
    if periodic
        aa  = ddx(rho(:,icx),h).*g0.*dt;
    else
        aa  = ddx(rho,h).*g0.*dt;
    end
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL; aa(:)];
end

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
if periodic
    jj1 = MapW(1:end-1,1:end-1); jj2 = MapW(1:end-1,2:end); jj3 = MapW(2:end,1:end-1); jj4 = MapW(2:end,2:end);
else
    jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);
end
IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL;-(1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL;+(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL;+(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL;-(1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
if periodic
    advn_mx = advect(rhofx.*U(2:end-1,:),(U(2:end-1,ifx(1:end-1))+U(2:end-1,ifx(2:end)))/2,(W(:,1:end-1)+W(:,2:end))/2,h,{ADVN,''},[1,2],BCA);
    advn_mx(:,[1 end]) = repmat((advn_mx(:,1)+advn_mx(:,end))/2,1,2);
    rr  = + (a2.*rhoUo+a3.*rhoUoo)/dt ...
          - advn_mx;
else
    advn_mx = advect(rhofx(:,2:end-1).*U(2:end-1,2:end-1),(U(2:end-1,1:end-1)+U(2:end-1,2:end))/2,(W(:,2:end-2)+W(:,3:end-1))/2,h,{ADVN,''},[1,2],BCA);
    rr  = + (a2.*rhoUo(:,2:end-1)+a3.*rhoUoo(:,2:end-1))/dt ...
          - advn_mx;
end
if bnchm
    if periodic
        rr = rr + src_U_mms(2:end-1,:); 
    else
        rr = rr + src_U_mms(2:end-1,2:end-1); 
    end
end

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
    if periodic
        ii  = MapU(2:end-1,:);
    else
        ii  = MapU(2:end-1,2:end-1);
    end
    
    %         left             ||           right
    if periodic
        jj1 = MapP(2:end-1,1:end-1); jj2 = MapP(2:end-1,2:end);
    else
        jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);
    end
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

if ~exist('KP','var') || bnchm || lambda1+lambda2>0
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
    
    ii  = [MapP(2:end-1,1    ); MapP(2:end-1,end)]; % left & right
    jj1 = ii;
    if periodic
        jj2 = [MapP(2:end-1,end-1); MapP(2:end-1,2    )];
    else
        jj2 = [MapP(2:end-1,2    ); MapP(2:end-1,end-1)];
    end
    
    aa  = zeros(size(ii));
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; aa(:)+1];
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; aa(:)-1];
    
    % internal points
    ii  = MapP(2:end-1,2:end-1);
    jj1 = MapP(1:end-2,2:end-1);
    jj2 = MapP(3:end-0,2:end-1);
    jj3 = MapP(2:end-1,1:end-2);
    jj4 = MapP(2:end-1,3:end-0);

    % coefficients multiplying matrix pressure P
    aa  = zeros(size(ii)) + lambda1*eps*h^2./eta;
    IIL = [IIL; ii(:)]; JJL = [JJL; ii(:)];    AAL = [AAL; aa(:)];  % P on stencil centre
    
    kP  = lambda2*h^2./eta;
    kP1 = (kP(icz(1:end-2),:).*kP(icz(2:end-1),:)).^0.5;   kP2 = (kP(icz(2:end-1),:).*kP(icz(3:end-0),:)).^0.5;
    kP3 = (kP(:,icx(1:end-2)).*kP(:,icx(2:end-1))).^0.5;   kP4 = (kP(:,icx(2:end-1)).*kP(:,icx(3:end-0))).^0.5;

    aa  = (kP1+kP2+kP3+kP4)/h^2;
    IIL = [IIL; ii(:)]; JJL = [JJL;  ii(:)];   AAL = [AAL;-aa(:)     ];      % P on stencil centre
    IIL = [IIL; ii(:)]; JJL = [JJL; jj1(:)];   AAL = [AAL; kP1(:)/h^2];      % P one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj2(:)];   AAL = [AAL; kP2(:)/h^2];      % P one below
    IIL = [IIL; ii(:)]; JJL = [JJL; jj3(:)];   AAL = [AAL; kP3(:)/h^2];      % P one above
    IIL = [IIL; ii(:)]; JJL = [JJL; jj4(:)];   AAL = [AAL; kP4(:)/h^2];      % P one below
    
    % assemble coefficient matrix
    KP  = sparse(IIL,JJL,AAL,NP,NP);
end

% RHS
IIR = [];       % equation indeces into R
AAR = [];       % forcing entries for R

ii  = MapP(2:end-1,2:end-1);

rr  = VolSrc;
if bnchm; rr = rr + src_P_mms(2:end-1,2:end-1); end

IIR = [IIR; ii(:)]; AAR = [AAR; rr(:)];

% assemble right-hand side vector
RP  = sparse(IIR,ones(size(IIR)),AAR,NP,1);

% set P = 0 in fixed point
nzp = round((Nz+2)/2);
nxp = round((Nx+2)/2);
np0 = MapP(nzp,nxp);
% DD(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),:) = 0;
KP(MapP(nzp,nxp),MapP(nzp,nxp)) = 1;
KP(MapP(end,nxp),MapP(nzp,nxp)) = 1;
% RP(MapP(nzp,nxp),:) = 0;

if bnchm; RP(MapP(nzp,nxp),:) = P_mms(nzp,nxp); end

if bnchm
    % set U = 0 in fixed point
    nzu = 1;
    nxu = round((Nx+2)/2);
    KV(MapU(nzu,nxu),:) = 0;
    GG(MapU(nzu,nxu),:) = 0;
    KV(MapU(nzu,nxu),MapU(nzu  ,nxu)) = 1;
    KV(MapU(nzu,nxu),MapU(nzu+1,nxu)) = 1;
    RV(MapU(nzp,nxp),:) = 0;
end

%% assemble and scale global coefficient matrix and right-hand side vector

LL  = [ KV   GG  ; ...
        DD   KP ];

RR  = [RV; RP];

SCL = (abs(diag(LL))).^0.5;
SCL = diag(sparse( 1./(SCL + sqrt(h^2./geomean(eta(:)))) ));

FF  = LL*[W(:);U(:);P(:)] - RR;

LL  = SCL*LL*SCL;
FF  = SCL*FF;
RR  = SCL*RR;


%% Solve linear system of equations for vx, vz, P

SOL = SCL*(LL\RR);  % update solution

% map solution vector to 2D arrays
W = full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2));  % matrix z-velocity
U = full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1));  % matrix x-velocity
P = full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2));  % matrix dynamic pressure

% upd_W = - alpha*full(reshape(SOL(MapW(:))        ,Nz+1,Nx+2)) + beta*upd_W;  % matrix z-velocity
% upd_U = - alpha*full(reshape(SOL(MapU(:))        ,Nz+2,Nx+1)) + beta*upd_U;  % matrix x-velocity
% upd_P = - alpha*full(reshape(SOL(MapP(:)+(NW+NU)),Nz+2,Nx+2)) + beta*upd_P;  % matrix dynamic pressure
% 
% W = W + upd_W;
% U = U + upd_U;
% P = P + upd_P;

end


if ~bnchm

    % phase segregation speeds
    wm(2:end-1,2:end-1) = ((rhom(1:end-1,:)+rhom(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_m(1:end-1,:).*Ksgr_m(2:end,:)).^0.5; % melt segregation speed
    wm([1,end],:) = min(1,1-[top;bot]).*wm([2,end-1],:);
    if periodic
        wm(:,[1 end]) = wm(:,[end-1 2]);
    else
        wm(:,[1 end]) = wm(:,[2 end-1]);
    end

    wx(2:end-1,2:end-1) = ((rhox(1:end-1,:)+rhox(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_x(1:end-1,:).*Ksgr_x(2:end,:)).^0.5; % solid segregation speed
    wx([1,end],:) = min(1,1-[top;bot]).*wx([2,end-1],:);
    if periodic
        wx(:,[1 end]) = wx(:,[end-1 2]);
    else
        wx(:,[1 end]) = wx(:,[2 end-1]);
    end

    wf(2:end-1,2:end-1) = ((rhof(1:end-1,:)+rhof(2:end,:))/2-mean(rhofz(2:end-1,:),2)).*g0.*(Ksgr_f(1:end-1,:).*Ksgr_f(2:end,:)).^0.5; % fluid segregation speed
    wf([1,end],:) = min(1,1-[top-fout;bot-fin]).*wf([2,end-1],:);
    if periodic
        wf(:,[1 end]) = wf(:,[end-1 2]);
    else
        wf(:,[1 end]) = wf(:,[2 end-1]);
    end

    % phase diffusion fluxes and speeds
    % [~,qmz,qmx] = diffus(mu ,km,h,[1,2],BCD);
    % [~,qxz,qxx] = diffus(chi,kx,h,[1,2],BCD);
    % [~,qfz,qfx] = diffus(phi,kf,h,[1,2],BCD);
    % qmz  = -qxz-qfz;  qmx = -qxx-qfx;

    dmudz = ddz( mu(icz,icx),h); dmudx = ddx( mu(icz,icx),h);
    dchdz = ddz(chi(icz,icx),h); dchdx = ddx(chi(icz,icx),h);
    dphdz = ddz(phi(icz,icx),h); dphdx = ddx(phi(icz,icx),h);

    kmz = (km(icz(1:end-1),icx)+km(icz(2:end),icx))./2;
    kmx = (km(icz,icx(1:end-1))+km(icz,icx(2:end)))./2;

    kxz = (kx(icz(1:end-1),icx)+kx(icz(2:end),icx))./2;
    kxx = (kx(icz,icx(1:end-1))+kx(icz,icx(2:end)))./2;

    kfz = (kf(icz(1:end-1),icx)+kf(icz(2:end),icx))./2;
    kfx = (kf(icz,icx(1:end-1))+kf(icz,icx(2:end)))./2;

    sumkz = kmz + kxz + kfz;
    sumkx = kmx + kxx + kfx;

    dstdz = kmz./sumkz .* dmudz + kxz./sumkz .* dchdz + kfz./sumkz .* dphdz;
    dstdx = kmx./sumkx .* dmudx + kxx./sumkx .* dchdx + kfx./sumkx .* dphdx;

    qmz = - kmz .* (dmudz-dstdz);  qmx = - kmx .* (dmudx-dstdx);
    qxz = - kxz .* (dchdz-dstdz);  qxx = - kxx .* (dchdx-dstdx);
    qfz = - kfz .* (dphdz-dstdz);  qfx = - kfx .* (dphdx-dstdx);

    phiz = (phi(icz(1:end-1),icx)+phi(icz(2:end),icx))./2;
    phix = (phi(icz,icx(1:end-1))+phi(icz,icx(2:end)))./2;

    chiz = (chi(icz(1:end-1),icx)+chi(icz(2:end),icx))./2;
    chix = (chi(icz,icx(1:end-1))+chi(icz,icx(2:end)))./2;

    muz  = (mu (icz(1:end-1),icx)+mu (icz(2:end),icx))./2;
    mux  = (mu (icz,icx(1:end-1))+mu (icz,icx(2:end)))./2;

    wqf = qfz./max(eps^0.5,phiz);  wqf([1,end],:) = min(1,1-[top-fout;bot-fin]).*wqf([2,end-1],:);
    uqf = qfx./max(eps^0.5,phix);

    wqx = qxz./max(eps^0.5,chiz);
    uqx = qxx./max(eps^0.5,chix);

    wqm = qmz./max(eps^0.5,muz);
    uqm = qmx./max(eps^0.5,mux);

    % update phase velocities
    Wf  = W + wf + wqf;  % mvp z-velocity
    Uf  = U + 0. + uqf;  % mvp x-velocity
    Wx  = W + wx + wqx;  % xtl z-velocity
    Ux  = U + 0. + uqx;  % xtl x-velocity
    Wm  = W + wm + wqm;  % mlt z-velocity
    Um  = U + 0. + uqm;  % mlt x-velocity

    
    %% update time step
    dtk = (h/2)^2/max([kc(:);kwm(:);kwx(:);kwf(:);(kT(:)+ks(:).*T(:))./rho(:)./cP(:)]); % diffusive time step size  
    dta =  h/2   /max(abs([Um(:).* mux(:);Wm(:).* muz(:); ...  % advective time step size
                           Ux(:).*chix(:);Wx(:).*chiz(:); ...
                           Uf(:).*phix(:);Wf(:).*phiz(:)]+eps));
    dtc = maxcmp./max(abs([advn_X(:)./rho(:);advn_M(:)./rho(:);advn_F(:)./rho(:)]));
    dt  = min([1.1*dto,min(CFL*[dtk,dta,dtc]),dtmax,tau_T/10]);                         % time step size
end

end

FMtime = FMtime + toc;