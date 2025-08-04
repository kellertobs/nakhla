% Identify convective layer boundaries
laybnd = zeros(Nz+1,Nx);
laybnd(2:end-1,:) = sqrt(eta(1:end-1,:).*eta(2:end,:))>etamax/10;
laybnd([1 end],:) = 1;
laybnd(2:end-1,:) = min(1,laybnd(2:end-1,:) + islocalmax(diff(rhop,1,1),1,'MinProminence',5,'MinSeparation',10));

% Create a matrix of row indices matching A.
rowind = repmat((1:(Nz+1))', [1, Nx]);

% Forward fill: for each column, record the last row index where A is nonzero.
lower = nan(size(laybnd));
lower(laybnd~=0) = rowind(laybnd~=0);
lower = fillmissing(lower, 'previous');

% Backward fill: for each column, record the next row index where laybnd is nonzero.
upper = nan(size(laybnd));
upper(laybnd~=0) = rowind(laybnd~=0);
upper = fillmissing(upper, 'next');

% B cell centers are at positions x = i + 0.5, for i = 1:nz.
ccpos = ((1:Nz)' + 0.5);

% For each column j, the lower candidate distance is:
%   d_lower = x - L(1:nz, j)
% and the upper candidate is:
%   d_upper = R(2:(Nz+1), j) - x.
d_lower = ccpos - lower(1:Nz, :);
d_upper = upper(2:(Nz+1), :) - ccpos;

% The nearest distance (in index units) is the minimum of the two.
dist = min(d_lower, d_upper);

% Convert to physical distance and cap at Delta_cnv0.
Delta_cnv = min(h * dist, Delta_cnv0);

% smooth correlation length to avoid sharp contrasts in regularisation
for i=1:10
    Delta_cnv = Delta_cnv + diffus(Delta_cnv,1/8*ones(size(Delta_cnv)),1,[1,2],BCD);
    Delta_cnv([1 end],:) = h/2;
end