% script to fit a polynomial function to the solidus and liquidus data from 
% Johnson et al. 2021:
% Johnson, T. E., Morrissey, L. J., Nemchin, A. A., Gardiner, N. J., & Snape, J. F. (2021). 
% The phases of the Moon: Modelling crystallisation of the lunar magma ocean through equilibrium thermodynamics. 
% Earth and Planetary Science Letters, 556, 116721. https://doi.org/10.1016/j.epsl.2020.116721

clear all; 

%% load data

sol = table2array( readtable("sol.xlsx",ReadVariableNames=false) ); % P and T for solidus curve
liq = table2array( readtable("liq.xlsx",ReadVariableNames=false) ); % P and T for liquidus curve

Tsol = sol(:,1);    Psol = sol(:,2);
Tliq = liq(:,1);    Pliq = liq(:,2);

%% do it using Simon's Law 
% preferred method for consistency with meltmodel function

simonlaw = @(x, P) x(1).*(1 + P./x(2)).^(1./x(3));
cfit_liq = lsqcurvefit(simonlaw, [Tliq(end),1,1], Pliq, Tliq);
cfit_sol = lsqcurvefit(simonlaw, [Tsol(end),1,1], Psol, Tsol);

figure;
set(gcf,'defaultaxescolororder', repelem(lines(2), 2, 1));
plot(Tsol, Psol, '+', 'markersize', 10); hold on; axis manual
plot(simonlaw(cfit_sol, Psol), Psol); 
plot(Tliq, Pliq, '+', 'markersize', 10); 
plot(simonlaw(cfit_liq, Pliq), Pliq);  hold off;
set(gca,'ydir', 'reverse'); xlim([1000,1900])
xlabel('Temperature [degC]');
ylabel('Pressure [GPa]');

save('solliq_johnson2021.mat', 'cfit_liq', 'cfit_sol', 'simonlaw');

%% other options I tested

%% use polynomial fitting
% but we need to split into deep and shallow parts because there are kinks
% in the curves

solcurve_deep = polyfit(Psol(Psol>1), Tsol(Psol>1), 2);
solcurve_shal = polyfit(Psol(Psol<1), Tsol(Psol<1), 2);
liqcurve_deep = polyfit(Pliq(Pliq>4), Tliq(Pliq>4), 2);
liqcurve_shal = polyfit(Pliq(Pliq<4), Tliq(Pliq<4), 2);

figure;
set(gcf,'defaultaxescolororder', repelem(lines(2), 3, 1));
plot(Tsol, Psol, '+', 'markersize', 10); hold on; axis manual
plot(polyval(solcurve_deep, Psol), Psol); 
plot(polyval(solcurve_shal, Psol), Psol); 
plot(Tliq, Pliq, '+', 'markersize', 10); 
plot(polyval(liqcurve_deep, Pliq), Pliq); 
plot(polyval(liqcurve_shal, Pliq), Pliq); hold off;
set(gca,'ydir', 'reverse'); xlim([1000,1900])
xlabel('Temperature [degC]');
ylabel('Pressure [GPa]');

%% use pchip to approximate piecewise polynomials - this is better

ppsol = pchip(Psol, Tsol);
ppliq = pchip(Pliq, Tliq);

figure;
set(gcf,'defaultaxescolororder', repelem(lines(2), 2, 1));
plot(Tsol, Psol, '+', 'markersize', 10); hold on; axis manual
plot(ppval(ppsol, Psol), Psol); 
plot(Tliq, Pliq, '+', 'markersize', 10); 
plot(ppval(ppliq, Pliq), Pliq);  hold off;
set(gca,'ydir', 'reverse'); xlim([1000,1900])
xlabel('Temperature [degC]');
ylabel('Pressure [GPa]');

% save('solliq_johnson2021.mat', 'ppsol', 'ppliq')


