%% assemble coefficients for matrix velocity diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R

% set cooling boundaries to no slip, else to free slip
if coolmode==3; sds = +1;      % no slip
else;           sds = -1; end  % free slip
if coolmode >0; top = +1;      % no slip
else;           top = -1; end  % free slip
if coolmode >1; bot = +1;      % no slip
else;           bot = -1; end  % free slip

% assemble coefficients of z-stress divergence
    
% left boundary
ii = MapW(:,1); jj1 = ii; jj2 = MapW(:,2);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+sds];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right boundary
ii = MapW(:,end); jj1 = ii; jj2 = MapW(:,end-1);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+sds];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% top boundary
ii = MapW(1,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - WBG(1,2:end-1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapW(end,2:end-1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - WBG(end,2:end-1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = MapW(2:end-1,2:end-1);
EtaC1 = etac(2:end-1,1:end-1); EtaC2 = etac(2:end-1,2:end  );
EtaP1 = eta (2:end-2,2:end-1); EtaP2 = eta (3:end-1,2:end-1);

% coefficients multiplying z-velocities W
%             top          ||         bottom          ||           left            ||          right
jj1 = MapW(1:end-2,2:end-1); jj2 = MapW(3:end,2:end-1); jj3 = MapW(2:end-1,1:end-2); jj4 = MapW(2:end-1,3:end);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)           ];      % W on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/h^2];      % W one above
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/h^2];      % W one below
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/h^2];      % W one to the left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/h^2];      % W one to the right

% coefficients multiplying x-velocities U
%         top left         ||        bottom left          ||       top right       ||       bottom right
jj1 = MapU(2:end-2,1:end-1); jj2 = MapU(3:end-1,1:end-1); jj3 = MapU(2:end-2,2:end); jj4 = MapU(3:end-1,2:end);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% z-RHS vector

rr = rhoBF .* g0;

IR = [IR; ii(:)];  RR = [RR; rr(:)];


%  assemble coefficients of x-stress divergence

% top boundary
ii = MapU(1,:); jj1 = ii; jj2 = MapU(2,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+top];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% bottom boundary
ii = MapU(end,:); jj1 = ii; jj2 = MapU(end-1,:);
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+bot];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% left side boundary
ii = MapU(2:end-1,1); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - UBG(2:end-1,1);
IR = [IR; ii(:)]; RR = [RR; aa(:)];

% right side boundary
ii = MapU(2:end-1,end); jj = ii;
aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj(:)];   AA = [AA; aa(:)+1];
aa = zeros(size(ii)) - UBG(2:end-1,end);
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii    = MapU(2:end-1,2:end-1);
EtaC1 = etac(1:end-1,2:end-1); EtaC2 = etac(2:end  ,2:end-1);
EtaP1 = eta (2:end-1,2:end-2); EtaP2 = eta (2:end-1,3:end-1);

% coefficients multiplying x-velocities U
%            left          ||          right          ||           top             ||          bottom
jj1 = MapU(2:end-1,1:end-2); jj2 = MapU(2:end-1,3:end); jj3 = MapU(1:end-2,2:end-1); jj4 = MapU(3:end,2:end-1);

aa = - 2/3*(EtaP1+EtaP2)/h^2 - 1/2*(EtaC1+EtaC2)/h^2;
II = [II; ii(:)]; JJ = [JJ;  ii(:)];   AA = [AA; aa(:)           ];      % U on stencil centre
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; 2/3*EtaP1(:)/h^2];      % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; 2/3*EtaP2(:)/h^2];      % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; 1/2*EtaC1(:)/h^2];      % U one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; 1/2*EtaC2(:)/h^2];      % U one below

% coefficients multiplying z-velocities W
%         top left         ||        top right          ||       bottom left       ||       bottom right
jj1 = MapW(1:end-1,2:end-2); jj2 = MapW(1:end-1,3:end-1); jj3 = MapW(2:end,2:end-2); jj4 = MapW(2:end,3:end-1);

II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; (1/2*EtaC1(:)-1/3*EtaP1(:))/h^2];  % W one to the top and left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA;-(1/2*EtaC1(:)-1/3*EtaP2(:))/h^2];  % W one to the top and right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA;-(1/2*EtaC2(:)-1/3*EtaP1(:))/h^2];  % W one to the bottom and left
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; (1/2*EtaC2(:)-1/3*EtaP2(:))/h^2];  % W one to the bottom and right


% x-RHS vector
rr = zeros(size(ii));
IR = [IR; ii(:)];  RR = [RR; rr(:)];


% assemble coefficient matrix & right-hand side vector
KV = sparse(II,JJ,AA,NW+NU,NW+NU);
RV = sparse(IR,ones(size(IR)),RR);


%% assemble coefficients for gradient operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A


% coefficients for z-gradient
ii = MapW(2:end-1,2:end-1);

%         top              ||          bottom
jj1 = MapP(2:end-2,2:end-1); jj2 = MapP(3:end-1,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];     % one to the top
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];     % one to the bottom


% coefficients for x-gradient
ii = MapU(2:end-1,2:end-1);

%         left             ||           right
jj1 = MapP(2:end-1,2:end-2); jj2 = MapP(2:end-1,3:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];     % one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];     % one to the right


% assemble coefficient matrix
GG = sparse(II,JJ,AA,NW+NU,NP);


%% assemble coefficients for divergence operator
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A

%internal points
ii = MapP(2:end-1,2:end-1);

% coefficients multiplying velocities U, W
%          left U          ||           right U       ||           top W           ||          bottom W
jj1 = MapU(2:end-1,1:end-1); jj2 = MapU(2:end-1,2:end); jj3 = MapW(1:end-1,2:end-1); jj4 = MapW(2:end,2:end-1);

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)-1/h];  % U one to the left
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)+1/h];  % U one to the right
II = [II; ii(:)]; JJ = [JJ; jj3(:)];   AA = [AA; aa(:)-1/h];  % W one above
II = [II; ii(:)]; JJ = [JJ; jj4(:)];   AA = [AA; aa(:)+1/h];  % W one below

% assemble coefficient matrix
DD = sparse(II,JJ,AA,NP,NW+NU);


%% assemble coefficients for matrix pressure diagonal and right-hand side
II  = [];       % equation indeces into A
JJ  = [];       % variable indeces into A
AA  = [];       % coefficients for A
IR  = [];       % equation indeces into R
RR  = [];       % forcing entries for R


% boundary points
ii  = [MapP(1,:).'; MapP(end  ,:).']; % top & bottom
jj1 = ii;
jj2 = [MapP(2,:).'; MapP(end-1,:).'];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];

ii  = [MapP(:,1); MapP(:,end  )]; % left & right
jj1 = ii;
jj2 = [MapP(:,2); MapP(:,end-1)];

aa = zeros(size(ii));
II = [II; ii(:)]; JJ = [JJ; jj1(:)];   AA = [AA; aa(:)+1];
II = [II; ii(:)]; JJ = [JJ; jj2(:)];   AA = [AA; aa(:)-1];
IR = [IR; ii(:)]; RR = [RR; aa(:)];


% internal points
ii = MapP(2:end-1,2:end-1);

% coefficients multiplying matrix pressure P
aa = gamma.*h^2./eta(2:end-1,2:end-1);
II = [II; ii(:)]; JJ = [JJ; ii(:)];    AA = [AA; aa(:)];  % P on stencil centre


% RHS
% rr = (theta.*VolSrc(2:end-1,2:end-1)+(1-theta).*VolSrco(2:end-1,2:end-1));
rr = VolSrc(2:end-1,2:end-1);
IR = [IR; ii(:)];
RR = [RR; rr(:)];


% assemble coefficient matrix and right-hand side vector
KP = sparse(II,JJ,AA,NP,NP);
RP = sparse(IR,ones(size(IR)),RR,NP,1);

np = round(N-2)/2+1;
KP(MapP(np,np),:) = 0;
KP(MapP(np,np),MapP(np,np)) = 1;
RP(MapP(np,np),:) = 0;


%% assemble global coefficient matrix and right-hand side vector
Pscale = 2*geomean(eta(:))/h;
LL = [-KV          Pscale.*GG  ; ...
       Pscale.*DD  Pscale.*KP  ];

RR = [RV; RP.*Pscale];


%% Scale system of equations (diagonal preconditioning)
% CC  =  sqrt(abs(diag(LL)));
% CC  =  diag(sparse(1./CC));
% 
% LL  =  CC*LL*CC;
% RR  =  CC*RR;


%% get residual
% get non-linear residual
% FF = LL*(CC\S) - RR;
FF      = LL*S - RR;
resnorm = norm(FF(:),2)./norm(RR(:),2);

% map residual vector to 2D arrays
res_W  = full(reshape(FF(MapW(:))        ,(Nz-1), Nx   ));                 % z-velocity
res_U  = full(reshape(FF(MapU(:))        , Nz   ,(Nx-1)));                 % x-velocity
res_P  = full(reshape(FF(MapP(:)+(NW+NU)), Nz   , Nx   ));                 % dynamic pressure
        

%% Solve linear system of equations for vx, vz, P
% S = CC*(LL\RR);  % update solution
S = LL\RR;  % update solution

% Read out solution
% map solution vector to 2D arrays
W  = full(reshape(S(MapW(:))        ,(Nz-1), Nx   ));                      % matrix z-velocity
U  = full(reshape(S(MapU(:))        , Nz   ,(Nx-1)));                      % matrix x-velocity
P  = full(reshape(S(MapP(:)+(NW+NU)), Nz   , Nx   )).*Pscale;              % matrix dynamic pressure

