% fractional crystallisation

if fractxtl && Nx==1 && Nz==1

    x = min(x,fractres); f = min(f,0.01); SUM = x+m+f;
    x = x./SUM;  m = m./SUM;  f = f./SUM;
    c = x.*cx + m.*cm + f.*cf;
    s = x.*sx + m.*sm + f.*sf;

    update;
    C = c.*rho;
    X = x.*rho;  M = m.*rho;  F = f.*rho;
    S = s.*rho;

    
% fractional melting  

elseif fractmlt && Nx==1 && Nz==1

    m = min(m,fractres); SUM = x+m+f;
    x = x./SUM;  m = m./SUM;  f = f./SUM;
    c = x.*cx + m.*cm + f.*cf;
    s = x.*sx + m.*sm + f.*sf;
    
    update;
    C = c.*rho;
    X = x.*rho;  M = m.*rho;  F = f.*rho;
    S = s.*rho;

end