function [dff,qz,qx] = diffus(f, k, h, dim, BC)
%
% [dff,qz,qx] = advect(f, k, h, dim, BC)
%
% example:
% dff = diffus(f, h, [1,2], {'closed','periodic'});
%
% diffuses a quantity f with a diffusion parameters k, for 2D problems
%
%
% the stencil covers 2 cell in each direction
%%
%                  i-1/2     i+1/2
%                             |         |
%         |   i-1   |    i    |   i+1   |
%       ..|----o----|----o----|----o----|...
%         |    fm   |   fcc   |    fp   |
%                             |         |
%
%
%
%
% INPUTS
% f         quantity you want to diffuse [Nz x Nx]
% k         diffusion parameter, collocated with f [Nz x Nx]
% h         grid spacing
% dim       which dimensions correspond to z,x [2-element array]
% BC        boundary conditions on z,x [2-element cell]
%
% OUTPUT
% dff       - div(- k grad(f))
% qz        - k df/dz
% qx        - k df/dx
%
% TK, 10 Feb 2023
% -------------------------------------------------------------------------

% collect information on the dimensions and BCs corresponding to (z,x)
zdim = dim(1);  zBC = BC{1};
xdim = dim(2);  xBC = BC{2};

% cell-centered values
fcc = f;
kcc = k;

% centered differences, prone to num dispersion
[fxm, fxp] = makestencil(f, xdim, xBC);
[fzm, fzp] = makestencil(f, zdim, zBC);

[kxm, kxp] = makestencil(k, xdim, xBC);
[kzm, kzp] = makestencil(k, zdim, zBC);

% get diffusive fluxes  q = - k grad(f)
qxp = - (kxp + kcc)/2 .* (fxp - fcc)/h;
qxm = - (kcc + kxm)/2 .* (fcc - fxm)/h;

qzp = - (kzp + kcc)/2 .* (fzp - fcc)/h;
qzm = - (kcc + kzm)/2 .* (fcc - fzm)/h;

% now calculate div(q)
dff = - (qxp - qxm)./h  ...
      - (qzp - qzm)./h;

if nargout>1
    % return diffusive fluxes on grid faces
    qz = zeros(size(f,dim(1))+1,size(f,dim(2))+2);
    qx = zeros(size(f,dim(1))+2,size(f,dim(2))+1);
    qz(1:end-1,2:end-1) = qzm;
    qz(end    ,2:end-1) = qzp(end,:);
    if strcmp(xBC,'periodic') && size(f,xdim)>1
        qz(:      ,[1,end]) = qz(:,[end-1,2]);
    else
        qz(:      ,[1,end]) = qz(:,[2,end-1]);
    end
    qx(2:end-1,1:end-1) = qxm;
    qx(2:end-1,end    ) = qxp(:,end);
    if strcmp(zBC,'periodic') && size(f,zdim)>1
        qx([1,end],:      ) = qx([end-1,2],:);
    else
        qx([1,end],:      ) = qx([2,end-1],:);
    end
end

end


%%  utility function

function [fm, fp] = makestencil (f, dim, BC)
%
% makes stencil for calculating differences
% use circshift which is faster than slicing
% default assumes periodic BCs, anything else just repeats boundary nodes
%
% |   i-1   |    i    |   i+1   |
% |----o----|----o----|----o----|
% |    fm   |   fcc   |    fp   |

% use these variables to handle arbitrary dim
shift = circshift([1, 0, 0], [0, dim-1]);

% 3 point stencil (used by most schemes)
fm  = circshift(f,    shift);
fp  = circshift(f, -  shift);


% if the boundary is closed or constant value
if ~strcmp(BC,'periodic') && size(f,dim)>1
    if ischar(BC)
        if dim==1
            fm (1  ,:,:) =  f( 1 ,:,:);
            fp (end,:,:) =  f(end,:,:);

        elseif dim==2
            fm (:,1  ,:) =  f(:, 1 ,:);
            fp (:,end,:) =  f(:,end,:);

        elseif dim==3
            fm (:,:,1  ) =  f(:,:, 1 );
            fp (:,:,end) =  f(:,:,end);

        end
    else
        if dim==1
            fm (1  ,:,:) = BC(1);
            fp (end,:,:) = BC(2);

        elseif dim==2
            fm (:,1  ,:) = BC(1);
            fp (:,end,:) = BC(2);

        elseif dim==3
            fm (:,:,1  ) = BC(1);
            fp (:,:,end) = BC(2);

        end
    end
end

end


