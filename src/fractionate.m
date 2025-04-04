if ~exist('CML','var'); CML.M = 0; end

% fractional crystallisation

if fractxtl && Nx==1 && Nz==1

    % update P,T-path
    Ptop = Ptop + (T-To).*dPdT;

    % record cumulate mass increments
    CML.r    (step)    = x-min(x,fractres);
    CML.M    (step)    = CML.r(step) .* (hist.sumB(1) - sum(CML.M(:)));
    CML.f    (step)    = CML.M(step)./hist.sumB(1);

    % record cumulate compositions
    CML.c    (step,:)  = cx;
    CML.c_oxd(step,:)  = cx_oxd;
    CML.c_mem(step,:)  = cx_mem;
    CML.c_msy(step,:)  = cx_msy;

    % record cumulate properties
    CML.rho  (step)    = rhox;
    CML.d     (step)   = CML.M(step)./CML.rho(step);  % relative cumulate depth
    CML.eta(step)      = etax;

    % remove cumulate fraction from system
    x = min(x,fractres); f = min(f,0.01); SUM = x+m+f;
    x = x./SUM;  m = m./SUM;  f = f./SUM;
    c = x.*cx + m.*cm + f.*cf;
    s = x.*sx + m.*sm + f.*sf;

    % update coefficients after cumulate removal
    update;
    C = c.*rho;
    X = x.*rho;  M = m.*rho;  F = f.*rho;
    S = s.*rho;

    
% fractional melting  

elseif fractmlt && Nx==1 && Nz==1

    % update P,T-path
    Ptop = Ptop + (T-To).*dPdT;

    % record cumulate mass increments
    CML.r    (step)    = m-min(m,fractres);
    CML.M    (step)    = CML.r(step) .* (hist.sumB(1) - sum(CML.M(:)));
    CML.f    (step)    = CML.M(step)./hist.sumB(1);

    % record cumulate compositions
    CML.c    (step,:)  = cm;
    CML.c_oxd(step,:)  = cm_oxd;
    CML.c_mem(step,:)  = cm_mem;
    CML.c_msy(step,:)  = cm_msy;

    % record cumulate properties
    CML.rho  (step)    = rhom;
    CML.d     (step)   = CML.M(step)./CML.rho(step);  % relative cumulate depth
    CML.eta(step)      = etam;

    % remove cumulate fraction from system
    m = min(m,fractres); SUM = x+m+f;
    x = x./SUM;  m = m./SUM;  f = f./SUM;
    c = x.*cx + m.*cm + f.*cf;
    s = x.*sx + m.*sm + f.*sf;
    
    % update coefficients after cumulate removal
    update;
    C = c.*rho;
    X = x.*rho;  M = m.*rho;  F = f.*rho;
    S = s.*rho;

end