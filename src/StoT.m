function [T,si] = StoT(T, s, P, f, cP, aT, bP, rho0, s0, T0, P0)

% Input dimensions
[nz, nx, ~] = size(f);
s0 = permute(repmat(s0,1,nz,nx),[2,3,1]);
aT = permute(repmat(aT,1,nz,nx),[2,3,1]);
bP = permute(repmat(bP,1,nz,nx),[2,3,1]);
cP = permute(repmat(cP,1,nz,nx),[2,3,1]);

% Convergence parameters
resnorm  = 1;
tol      = 1e-9;
max_iter = 100;
iter     = 0;

while resnorm >= tol && iter < max_iter
    iter = iter + 1;

    % Intermediate terms with broadcasting
    a = aT .* (T-T0);            % 3D: (nz, nx, nphs)
    b = bP .* (P-P0);            % 3D: (nz, nx, nphs)

    % Phase-specific entropy s^i
    si = s0 + cP.*log(T/T0) - (aT./(bP.*rho0)) .* log((1-a+b)./(1-a));

    % Residual f(T)
    r = sum(f .* si, 3) - s;  % 2D: (nz, nx)

    % Derivative ds^i/dT        
    dsi_dT = cP./T - (aT.^2.*(P-P0))./(rho0.*(1-a).*(1-a+b));

    % Derivative f'(T)
    dr_dT = sum(f .* dsi_dT, 3);    % 2D: (nz, nx)

    % Newton update
    T = T - r ./ dr_dT;

    % residual norm
    resnorm = norm(r./dr_dT)./norm(T);

end

Tp = T;

while resnorm >= tol && iter < max_iter
    iter = iter + 1;

    % Intermediate terms with broadcasting
    a = aT .* (Tp-T0);           % 3D: (nz, nx, nphs)

    % Phase-specific entropy s^i
    spi = s0 + cP.*log(Tp/T0);

    % Residual f(T)
    r = sum(f .* si, 3) - sum(f .* spi, 3);  % 2D: (nz, nx)

    % Derivative ds^i/dT        
    dspi_dT = cP./Tp;

    % Derivative f'(T)
    dr_dT = sum(f .* dspi_dT, 3);    % 2D: (nz, nx)

    % Newton update
    T = T - r ./ dr_dT;

    % residual norm
    resnorm = norm(r./dr_dT)./norm(T);

end

% Convergence feedback
if iter >= max_iter
    warning('Maximum iterations reached. Solution may not have converged.');
end
end