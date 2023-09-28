function [Tsol, Tliq] = solidusliquidus(model, Psol, Pliq)

if nargin<3, Pliq = Psol; end

switch model

    case 'johnson2021'
        % Johnson, T. E., Morrissey, L. J., Nemchin, A. A., Gardiner, N. J.,
        % & Snape, J. F. (2021). 
        % The phases of the Moon: Modelling crystallisation of the lunar 
        % magma ocean through equilibrium thermodynamics. 
        % Earth and Planetary Science Letters, 556, 116721. 
        % https://doi.org/10.1016/j.epsl.2020.116721

        load("fitsolliq/solliq_johnson2021.mat", 'cfit_liq', 'cfit_sol', 'simonlaw');
        Tsol = simonlaw(cfit_sol, Psol);
        Tliq = simonlaw(cfit_liq, Pliq);

    %case other paper

end

end
