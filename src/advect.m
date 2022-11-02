function [adv, advscl] = advect (f, u, w, h, scheme, dim, BC)
%
% [adv, advscl] = advect(f, u, w, h, scheme, dim, BC)
%
% example:
% adv = advect(f, u, w, h, {'quick','vdf'}, [2,3], {'periodic','periodic'});
%
% possible advection schemes :
%       'centr', 'upwd1', 'quick', 'fromm', 'weno3', 'weno5', 'tvdim'
%
% advects a quantity f on a velocity field (u,w), for 2D problems
%
%
% Depending on the scheme, the stencil may cover 2--5 cells, and we split
% the flux into positive and negative components to ensure upwind bias, so
% that the equation of div(v f) on the i-th node is
%
% div(v f) = 1/h x { [uppos fppos + upneg fpneg] -  [umpos fmpos + umneg fmneg] }
%
%                            i-1/2     i+1/2
%                             |         |
%         |   i-2   |   i-1   |    i    |   i+1   |   i+2   |
%       ..|----o----|----o----|----o----|----o----|----o----|...
%         |   fmm   |    fm   |   fcc   |    fp   |   fpp   |
%                             |         |
% +flux:           ---------umpos-----umpos-------->
%                           fmpos     fppos
%                             |         |
% -flux:           <--------umneg-----upneg---------
%                           fmneg     fpneg
%
%
%
%
% INPUTS
% f         quantity you want to advect [Nz x Nx, or Nphs x Nz x Nx]
% u         STAGGERED horiz velocity field [Nz x Nx+1, or Nphs x Nz x Nx+1]
% w         STAGGERED vert  velocity field [Nz+1 x Nx, or Nphs x Nz+1 x Nx]
% h         grid spacing
% scheme    info about advection scheme [1 or 2-element cell]
%               1st element is the advection scheme name
%               if 2nd element='vdf', returns v grad(f), otherwise returns div(f v)
% dim       which dimensions correspond to z,x [2-element array]
% BC        boundary conditions on z,x [2-element cell]
%
% OUTPUT
% adv       div(v x f) or v grad(f) depending on scheme{2}, same size as f
%
% YQW, 11 Oct 2022
% -------------------------------------------------------------------------

% collect information on the dimensions and BCs corresponding to (z,x)
zdim = dim(1);  zBC = BC{1};
xdim = dim(2);  xBC = BC{2};

% collect velocities on - (m) and + (p) faces
[umpos, umneg, uppos, upneg] = facevels(u, xdim);
[wmpos, wmneg, wppos, wpneg] = facevels(w, zdim);

% cell-centered value
fcc = f;

% in this switch, define phase fractions at the cell faces, split into
% positive and negative flux components
switch scheme{1}
    case 'centr'
        % centered differences, prone to num dispersion
        [fxm, fxp] = makestencil(f, xdim, xBC);
        [fzm, fzp] = makestencil(f, zdim, zBC);
        
        fxppos = (fcc+fxp)./2;      fxpneg = (fcc+fxp)./2;
        fxmpos = (fcc+fxm)./2;      fxmneg = (fcc+fxm)./2;
        fzppos = (fcc+fzp)./2;      fzpneg = (fcc+fzp)./2;
        fzmpos = (fcc+fzm)./2;      fzmneg = (fcc+fzm)./2;
        
    case 'upwd1'
        % upwind differences, prone to num diffusion
        % flux conservative approach, split velocities into + and -
        [fxm, fxp] = makestencil(f, xdim, xBC);
        [fzm, fzp] = makestencil(f, zdim, zBC);
        
        fxppos = fcc;      fxpneg = fxp;
        fxmpos = fxm;      fxmneg = fcc;
        fzppos = fcc;      fzpneg = fzp;
        fzmpos = fzm;      fzmneg = fcc;
        
    case 'quick'
        % quick scheme == 3rd order upwind
        % flux conservative approach, split velocities into + and -
        [fxm, fxp, fmmx, fppx] = makestencil(f, xdim, xBC);
        [fzm, fzp, fmmz, fppz] = makestencil(f, zdim, zBC);
        
        fxppos = (2*fxp + 5*fcc - fxm )./6;      fxpneg = (2*fcc + 5*fxp - fppx)./6;
        fxmpos = (2*fcc + 5*fxm - fmmx)./6;      fxmneg = (2*fxm + 5*fcc - fxp )./6;
        fzppos = (2*fzp + 5*fcc - fzm )./6;      fzpneg = (2*fcc + 5*fzp - fppz)./6;
        fzmpos = (2*fcc + 5*fzm - fmmz)./6;      fzmneg = (2*fzm + 5*fcc - fzp )./6;
        
    case 'fromm'
        % Fromm scheme
        [fxm, fxp, fmmx, fppx] = makestencil(f, xdim, xBC);
        [fzm, fzp, fmmz, fppz] = makestencil(f, zdim, zBC);
        
        fxppos = fcc + (fxp-fxm )./4;      fxpneg = fxp + (fcc-fppx)./4;
        fxmpos = fxm + (fcc-fmmx)./4;      fxmneg = fcc + (fxm-fxp )./4;
        fzppos = fcc + (fzp-fzm )./4 ;     fzpneg = fzp + (fcc-fppz)./4;
        fzmpos = fzm + (fcc-fmmz)./4;      fzmneg = fcc + (fzm-fzp )./4;
        
    case 'weno3'
        % 3rd order WENO from Jiang & Shu 1996, J Comp Physics
        [fxppos, fxpneg, fxmpos, fxmneg] = weno3(fcc, xdim, xBC);
        [fzppos, fzpneg, fzmpos, fzmneg] = weno3(fcc, zdim, zBC);
        
    case 'weno5'
        % 5th order WENO from Jiang & Shu 1996, J Comp Physics
        [fxppos, fxpneg, fxmpos, fxmneg] = weno5(fcc, xdim, xBC);
        [fzppos, fzpneg, fzmpos, fzmneg] = weno5(fcc, zdim, zBC);
        
    case 'tvdim'
        % total variation diminishing approach from Sramek et al. 2010, GJI
        [fxppos, fxpneg, fxmpos, fxmneg] = tvd(fcc, umpos, umneg, uppos, upneg, xdim, xBC);
        [fzppos, fzpneg, fzmpos, fzmneg] = tvd(fcc, wmpos, wmneg, wppos, wpneg, zdim, zBC);
end


% now calculate div(f v)
adv = (uppos.*fxppos + upneg.*fxpneg - umpos.*fxmpos - umneg.*fxmneg)./h + ...
      (wppos.*fzppos + wpneg.*fzpneg - wmpos.*fzmpos - wmneg.*fzmneg)./h;

% if you only want the advection term, remove f x div(v)
if strcmp(scheme{2}, 'vdf')
    % calculate f x div(v)
    fdv = f.*(diff(u,1,xdim)./h + diff(w,1,zdim)./h);
    
    % v x grad(f) = div (f v) - f x div(v)
    adv = adv - fdv;
end

if nargout>1
    % return advection scale for calculating time step
    advscl = max(abs([f.*(uppos+upneg); f.*(wppos+wpneg)]), [], 'all');
end
end



%%  utility functions used for many schemes

function [vmpos, vmneg, vppos, vpneg] = facevels (v, dim)
% return - and + components of cell face velocities separated into vm, vp

if     dim==1, vm = v(1:end-1,:,:);    vp = v(2:end,:,:);
elseif dim==2, vm = v(:,1:end-1,:);    vp = v(:,2:end,:);
elseif dim==3, vm = v(:,:,1:end-1);    vp = v(:,:,2:end);
end

% now split the velocities into positive and negative
vmpos = 0.5*(vm + abs(vm));    % positive velocity
vmneg = 0.5*(vm - abs(vm));    % negative velocity

vppos = 0.5*(vp + abs(vp));    % positive velocity
vpneg = 0.5*(vp - abs(vp));    % negative velocity
end



function [fm, fp, fmm, fpp, fmmm, fppp] = makestencil (f, dim, BC)
%
% makes stencil for calculating differences
% use circshift which is faster than slicing
% default assumes periodic BCs, anything else just repeats boundary nodes
%
% |   i-3   |   i-2   |   i-1   |    i    |   i+1   |   i+2   |   i+3   |
% |----o----|----o----|----o----|----o----|----o----|----o----|----o----|
% |   fmmm  |   fmm   |    fm   |   fcc   |    fp   |   fpp   |   fppp  |

sten7 = (nargout>4) ;

% use these variables to handle arbitrary dim
shift = circshift([1, 0, 0], [0, dim-1]);

% 5 point stencil (used by most schemes)
fmm = circshift(f,  2*shift);
fm  = circshift(f,    shift);
fp  = circshift(f, -  shift);
fpp = circshift(f, -2*shift);

if (sten7)
    % extras for 7 point stencil (only used by WENO5)
    fmmm = circshift(f,  3*shift);
    fppp = circshift(f, -3*shift);
end


% if the boundary is closed or open (could use more elegance)
% but i wanted to be able to handle any specified dimension
if ~strcmp(BC,'periodic') && size(f,dim)>1
    if ischar(BC)
        if dim==1
            fmm(    1:2  ,:,:) = repmat(f( 1 ,:,:),2,1,1);
            fm (    1    ,:,:) =        f( 1 ,:,:);
            fp (      end,:,:) =        f(end,:,:);
            fpp(end-1:end,:,:) = repmat(f(end,:,:),2,1,1);

            if (sten7)
                fmmm(    1:3  ,:,:) = repmat(f( 1 ,:,:),3,1,1);
                fppp(end-2:end,:,:) = repmat(f(end,:,:),3,1,1);
            end
        elseif dim==2
            fmm(:,    1:2  ,:) = repmat(f(:, 1 ,:),1,2,1);
            fm (:,    1    ,:) =        f(:, 1 ,:);
            fp (:,      end,:) =        f(:,end,:);
            fpp(:,end-1:end,:) = repmat(f(:,end,:),1,2,1);

            if (sten7)
                fmmm(:,    1:3  ,:) = repmat(f(:, 1 ,:),1,3,1);
                fppp(:,end-2:end,:) = repmat(f(:,end,:),1,3,1);
            end
        elseif dim==3
            fmm(:,:,    1:2  ) = repmat(f(:,:, 1 ),1,1,2);
            fm (:,:,    1    ) =        f(:,:, 1 );
            fp (:,:,      end) =        f(:,:,end);
            fpp(:,:,end-1:end) = repmat(f(:,:,end),1,1,2);

            if (sten7)
                fmmm(:,:,    1:3  ) = repmat(f(:,:, 1 ),1,1,3);
                fppp(:,:,end-2:end) = repmat(f(:,:,end),1,1,3);
            end
        end
    else
        if dim==1
            fmm(    1:2  ,:,:) = BC(1);
            fm (    1    ,:,:) = BC(1);
            fp (      end,:,:) = BC(2);
            fpp(end-1:end,:,:) = BC(2);

            if (sten7)
                fmmm(    1:3  ,:,:) = BC(1);
                fppp(end-2:end,:,:) = BC(2);
            end
        elseif dim==2
            fmm(:,    1:2  ,:) = BC(1);
            fm (:,    1    ,:) = BC(1);
            fp (:,      end,:) = BC(2);
            fpp(:,end-1:end,:) = BC(2);

            if (sten7)
                fmmm(:,    1:3  ,:) = BC(1);
                fppp(:,end-2:end,:) = BC(2);
            end
        elseif dim==3
            fmm(:,:,    1:2  ) = BC(1);
            fm (:,:,    1    ) = BC(1);
            fp (:,:,      end) = BC(2);
            fpp(:,:,end-1:end) = BC(2);

            if (sten7)
                fmmm(:,:,    1:3  ) = BC(1);
                fppp(:,:,end-2:end) = BC(2);
            end
        end
    end
end
end




%% weno3 functions

function [fppos, fpneg, fmpos, fmneg] = weno3 (f, dim, BC)

[fm, fp, fmm, fpp] = makestencil(f, dim, BC);

% for the p (i+1/2) face
fppos = makeweno3poly( fm, f, fp);   % +flux upwind polynomials
fpneg = makeweno3poly(fpp, fp, f);   % -flux upwind polynomials

% for the m (i-1/2) face
fmpos = makeweno3poly(fmm, fm, f);   % +flux upwind polynomials
fmneg = makeweno3poly( fp, f, fm);   % -flux upwind polynomials

end

function [fface] = makeweno3poly (fm, fc, fp)
% 3rd order polynomials, also defined in Jiang & Shu, 1996, J Comp Physics

% linear polynomials
p1 = 0.5*(  fc + fp);
p2 = 0.5*(3*fc - fm);

% smoothness measure
b1 = (p1 - fc).^2;
b2 = (p2 - fc).^2;

% weights
g   = [1/3, 2/3];
eps = 1e-16;    %stabiliser
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);

fface = ( wp1.*p1 + wp2.*p2 ) ./ (wp1 + wp2);

end


%% weno5 functions

function [fppos, fpneg, fmpos, fmneg] = weno5 (f, dim, BC)

[fm, fp, fmm, fpp, fmmm, fppp] = makestencil(f, dim, BC);

% for the p (i+1/2) face
fppos = makeweno5poly(fmm ,  fm, f, fp, fpp);    % +flux upwind polynomials
fpneg = makeweno5poly(fppp, fpp, fp, f, fm );    % -flux upwind polynomials

% for the m (i-1/2) face
fmpos = makeweno5poly(fmmm, fmm, fm, f , fp );   % +flux upwind polynomials
fmneg = makeweno5poly( fpp,  fp, f , fm, fmm);   % -flux upwind polynomials
end

function [fface] = makeweno5poly (fmm, fm, fc, fp, fpp)
% 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics

% 5th order polynomials
p1 = (2*fmm - 7*fm + 11*fc )/6;
p2 = ( -fm  + 5*fc +  2*fp )/6;
p3 = (2*fc  + 5*fp -    fpp)/6;

% smoothness measure (special from this paper)
b1 = 13/12*(fmm - 2*fm + fc ).^2 + 1/4*(  fmm - 4*fm + 3*fc ).^2;
b2 = 13/12*(fm  - 2*fc + fp ).^2 + 1/4*(  fm  -          fp ).^2;
b3 = 13/12*(fc  - 2*fp + fpp).^2 + 1/4*(3*fc  - 4*fp +   fpp).^2;

% weights
g   = [1/10, 6/10, 3/10];
eps = 1e-16;    %stabiliser
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);
wp3 = g(3)./(b3.^2 + eps);

% get flux (normalize weights at the same time)
fface = (wp1.*p1 + wp2.*p2 + wp3.*p3) ./ (wp1 + wp2 + wp3) ;

end



%% tvd function
% total variation diminishing approach as described in Sramek 2010

function [fppos, fpneg, fmpos, fmneg] = tvd (f, vmpos, vmneg, vppos, vpneg, dim, BC)

[fm, fp, fmm, fpp] = makestencil(f, dim, BC);

% we need extra velocity terms from one cell upwind
[vmmpos   ] = makestencil(vmpos, dim, BC);
[~, vppneg] = makestencil(vpneg, dim, BC);

% for the p (i+1/2) face
fppos = gettvdflux( fm, f , fp,  vmpos, vppos);    % +flux
fpneg = gettvdflux(fpp, fp, f , vppneg, vpneg);    % -flux

% for the m (i-1/2) face
fmpos = gettvdflux(fmm, fm, f , vmmpos, vmpos);    % +flux
fmneg = gettvdflux( fp, f , fm, vpneg , vmneg);    % -flux

end

function [fface] = gettvdflux (fm, fc, fp, vm, vp)

% ratio of flux in one cell upwind / current cell
R = ( vm.*(fc-fm))./( vp.*(fp-fc) );

% choose lambda (flux limiter)
% l = max(0, min(1,Rpos));  % minmod approach
l = max( zeros(size(R)), max(min(1,2*R), min(2,R)) );  % superbee scheme

% get flux on cell face
fface = fc  + 0.5*l.*(fp - fc);

end

