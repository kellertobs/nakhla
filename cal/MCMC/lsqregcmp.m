% lsqregcmp Linear least squares with regularisation and nonnegativity and
%   unit sum constraints
%   x = lsqregcmp(A,b,lambda) returns the vector x that minimizes 
%   the regularised problem NORM(A*x-b) + lambda(1)*norm(x) + lambda(2)*D*x
%   subject to constraints X >= 0, and sum(X) = 1. 
% 
%   If b is a matrix, the solution above is returned line by line to produce 
%   a solution matrix x. To introduce regularisation of the solution along 
%   the second dimension of x, a Laplacian smoothing step is applied as
%   x = x + lambda(3)/4*diff(x,2)

function [x] = lsqregcmp(A,b,lambda)

% fill in unspecified regularisation parameters
if     length(lambda)<2
    lambda(2:4) = 0;
elseif length(lambda)<3
    lambda(3:4) = 0;
elseif length(lambda)<4
    lambda(  4) = 0;
end

% get size of arrays
[n,m] = size(A);
[p,q] = size(b);

% perform least squares fit
x = zeros(p,m);
for ip = 1:p
    bp      = b(ip,:).';
    D       = diff(eye(m), 1, 1);
    A_aug   = [A; ones(1,m); sqrt(lambda(1)) * eye(m); sqrt(lambda(2)) * D];
    b_aug   = [bp; 1; zeros(m, 1); zeros(m-1, 1)];
    x(ip,:) = lsqnonneg(A_aug,b_aug);
end

% additional regularisation and final normalisation
x = max(lambda(4), x + lambda(3)/4*diff(x([1 1:end end],:),2,1));
x = x./(sum(x,2)+eps);

end