function [adv, advscl] = advect (f, u, w, h, scheme, dim, BC)
%
% [adv] = advect (f, u, w, h, scheme, dim, BC)
%
% example:
% adv = advect(f, u, w, h, {'quick','vdf'}, [2,3], {'periodic','periodic'});
% 
% advects a quantity f on a velocity field represented by (u,w),
% for 2D problems
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
% INPUTS
% f         quantity you want to advect [Nz x Nx, or Nphs x Nz x Nx]
% u         STAGGERED horiz velocity field [Nz x Nx+1, or Nphs x Nz x Nx+1]
% w         STAGGERED vert  velocity field [Nz+1 x Nx, or Nphs x Nz+1 x Nx]
% h         grid spacing
% scheme    info about advection scheme [2-element cell]
%               1st element is the advection scheme name
%               2nd element is '' (returns div(v x f)) or 'vdf' (returns v grad(f))
% dim       which dimensions correspond to z, x [2-element array]
% BC        boundary conditions on z,x [2-element cell]
% 
% OUTPUT
% adv       div(v x f) or v grad(f) depending on scheme{2}, same size as f
% 
% YQW, 11 Oct 2022
% 



% collect information on the dimensions and BCs corresponding to (z,x)
zdim = dim(1);  zBC = BC{1};
xdim = dim(2);  xBC = BC{2};

% calculate f x div(v)
dudx = 0; dwdz = 0;
if size(f,xdim)>1, dudx = diff(u,1,xdim)./h; end
if size(f,zdim)>1, dwdz = diff(w,1,zdim)./h; end
fdv = f.*(dudx + dwdz);   

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

          
    % the schemes below here definitely work for periodic BCs, 
    % still not sure about closed BCs
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
        [fxppos, fxpneg, fxmpos, fxmneg] = tvd(fcc, uppos, upneg, xdim, xBC);
        [fzppos, fzpneg, fzmpos, fzmneg] = tvd(fcc, wppos, wpneg, zdim, zBC);
end


% now calculate div(v x f)
duf = 0; dwf = 0;
if size(f,xdim)>1, duf = (uppos.*fxppos + upneg.*fxpneg - umpos.*fxmpos - umneg.*fxmneg)./h; end
if size(f,zdim)>1, dwf = (wppos.*fzppos + wpneg.*fzpneg - wmpos.*fzmpos - wmneg.*fzmneg)./h; end
adv = duf + dwf;

if strcmp(scheme{2}, 'vdf')
    % v x grad(f) = div (v x f) - f x div(v)
    adv = adv - fdv;
end

if nargout>1
    % return advection scale for calculating time step
    advscl = max(abs([f.*(uppos+upneg); f.*(wppos+wpneg)]), [], 'all');
end
end



%%  utility functions used for many schemes

function [vmpos, vmneg, vppos, vpneg] = facevels (v, dim)
% return cell face velocities - (vm) and + (vp) of cell center

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

function [flux] = shiftflux (flux, shift, dim, BC)
% function needed to deal with circshift when applied to closed BCs
% fluxes = 0 at closed boundaries 

flux = circshift(flux, shift);

if strcmp(BC, 'closed')
    % which direction is the shift? + is right, - is left
    dir = sign(sum(shift)); 
    if dir>0
        if     dim==1, flux(1,:,:) = 0;
        elseif dim==2, flux(:,1,:) = 0;
        elseif dim==3, flux(:,:,1) = 0;
        end
    else
        if     dim==1, flux(end,:,:) = 0;
        elseif dim==2, flux(:,end,:) = 0;
        elseif dim==3, flux(:,:,end) = 0;
        end
    end
end

end

function [fm, fp, fmm, fpp, fppp] = makestencil (f, dim, BC)
% 
% makes stencil for calculating differences
% use circshift which is faster than slicing
% 
%         |   i-2   |   i-1   |    i    |   i+1   |   i+2   |
%       ..|----o----|----o----|----o----|----o----|----o----|...
%         |   fmm   |    fm   |   fcc   |    fp   |   fpp   |)

shift = circshift([1, 0, 0], [0, dim-1]);
sten5 = (nargout>3) ;

% 3 point stencil
fm = circshift(f,  shift);
fp = circshift(f, -shift);

if (sten5) 
    % extras for 5 point stencil
    fmm  = circshift(f,  2*shift);
    fpp  = circshift(f, -2*shift);
    fppp = circshift(f, -3*shift);
end

% if the boundary is closed (could use more elegance) but i wanted to be
% able to handle any specified dimension
if strcmp(BC,'closed') && size(f,dim)>1
    if dim==1 
        fm( 1 ,:,:) = f( 1 ,:,:);
        fp(end,:,:) = f(end,:,:);
        
        if (sten5)
            fmm (    1:2  ,:,:) = repmat(f( 1 ,:,:),2,1,1);
            fpp (end-1:end,:,:) = repmat(f(end,:,:),2,1,1);
            fppp(end-2:end,:,:) = repmat(f(end,:,:),3,1,1);
        end
    elseif dim==2
        fm(:, 1 ,:) = f(:, 1 ,:);
        fp(:,end,:) = f(:,end,:);
        
        if (sten5)
            fmm (:,    1:2  ,:) = repmat(f(:, 1 ,:),1,2,1);
            fpp (:,end-1:end,:) = repmat(f(:,end,:),1,2,1);
            fppp(:,end-2:end,:) = repmat(f(:,end,:),1,3,1);
        end
    elseif dim==3
        fm(:,:, 1 ) = f(:,:, 1 );
        fp(:,:,end) = f(:,:,end);
        
        if (sten5)
            fmm (:,:,    1:2  ) = repmat(f(:,:, 1 ),1,1,2);
            fpp (:,:,end-1:end) = repmat(f(:,:,end),1,1,2);
            fppp(:,:,end-2:end) = repmat(f(:,:,end),1,1,3);
        end
    end
    
end
end




%% weno3 functions

function [fppos, fpneg, fmpos, fmneg] = weno3 (f, dim, BC)

% define shifting dimension 
shift = circshift([1, 0, 0], [0, dim-1]);

% smaller stencil than weno 5
[fm, fp, ~, fpp, ~] = makestencil(f, dim, BC);
fppos = makeweno3poly( fm, f, fp);   %  left upwind polynomials
fpneg = makeweno3poly(fpp, fp, f);   % right upwind polynomials

% get the fractions for left cell face 
% (cheat, just set boundaries to 0 for closed BCs)
fmpos = shiftflux(fppos, shift, dim, BC);
fmneg = shiftflux(fpneg, shift, dim, BC);

end

function [fhalf] = makeweno3poly (fw, fc, fe)
% 3rd order polynomials, also defined in Jiang & Shu, 1996, J Comp Physics

% linear polynomials
p1 = 0.5*(  fc + fe);
p2 = 0.5*(3*fc - fw);

% smoothness measure
b1 = (p1 - fc).^2;
b2 = (p2 - fc).^2;

% weights
g   = [1/3, 2/3]; 
eps = 1e-16;    %stabiliser
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);

fhalf = ( wp1.*p1 + wp2.*p2 ) ./ (wp1 + wp2);

end


%% weno5 functions
function [fppos, fpneg, fmpos, fmneg] = weno5 (f, dim, BC)

% define shifting dimension 
shift = circshift([1, 0, 0], [0, dim-1]);

[fm, fp, fmm, fpp, fppp] = makestencil(f, dim, BC);
fppos = makeweno5poly(fmm ,  fm, f, fp, fpp);   %  left upwind polynomials
fpneg = makeweno5poly(fppp, fpp, fp, f, fm );   % right upwind polynomials

% get the fractions for left cell face 
% (cheat, just set boundaries to 0 for closed BCs)
fmpos = shiftflux(fppos, shift, dim, BC);
fmneg = shiftflux(fpneg, shift, dim, BC);
end

function [fhalf] = makeweno5poly (fww, fw, fc, fe, fee)
% 5th order WENO polynomials from Jiang & Shu, 1996, J Comp Physics

% 5th order polynomials
p1 = (2*fww - 7*fw + 11*fc )/6;
p2 = ( -fw  + 5*fc +  2*fe )/6;
p3 = (2*fc  + 5*fe -    fee)/6;

% smoothness measure (special from this paper)
b1 = 13/12*(fww - 2*fw + fc ).^2 + 1/4*(  fww - 4*fw + 3*fc ).^2; 
b2 = 13/12*(fw  - 2*fc + fe ).^2 + 1/4*(  fw  -          fe ).^2;
b3 = 13/12*(fc  - 2*fe + fee).^2 + 1/4*(3*fc  - 4*fe +   fee).^2;

% weights
g   = [1/10, 6/10, 3/10]; 
eps = 1e-16;    %stabiliser
wp1 = g(1)./(b1.^2 + eps);
wp2 = g(2)./(b2.^2 + eps);
wp3 = g(3)./(b3.^2 + eps);

% get flux (normalize weights at the same time)
fhalf = (wp1.*p1 + wp2.*p2 + wp3.*p3) ./ (wp1 + wp2 + wp3) ;

end



%% tvd function

function [fppos, fpneg, fmpos, fmneg] = tvd (f, vpos, vneg, dim, BC)
% total variation diminishing approach as described in Sramek 2010
% flux is split to ensure upwind-bias.
% vpos for positive fluxes; vneg for negative fluxes

% define shifting dimension 
shift = circshift([1, 0, 0], [0, dim-1]);

[fm, fp, ~, fpp] = makestencil(f, dim, BC);

% out/influx ratios on the cells that share the east face
Rpos = ( shiftflux(vpos, shift,dim,BC).*(f  -fm))./( vpos.*(fp-f) );
Rneg = ( shiftflux(vneg,-shift,dim,BC).*(fpp-fp))./( vneg.*(fp-f) );

% use the flux ratios to retrieve weights
% minmod approach
% lpos = max(0, min(1,Rpos));
% lneg = max(0, min(1,Rneg));

% superbee scheme
lpos = max( zeros(size(Rpos)), max(min(1,2*Rpos), min(2, Rpos)) );
lneg = max( zeros(size(Rneg)), max(min(1,2*Rneg), min(2, Rneg)) );

% calculate + and - phase fractions on the east face
fppos = f  + 0.5*lpos.*(fp - f );
fpneg = fp + 0.5*lneg.*(f  - fp);

% get the fractions for left cell face 
% (cheat, just set boundaries to 0 for closed BCs)
fmpos = shiftflux(fppos, shift, dim, BC);
fmneg = shiftflux(fpneg, shift, dim, BC);

end

