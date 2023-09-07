function [datafit,var] = OxdFromMeltTemp(model,T,P,H2Osat,c,cal)

cal.T0 = model(1:cal.ncmp-1).';
cal.r  = model(cal.ncmp:end).';
cal.A  = (cal.T0+273.15)./300;

var.m = 1.0; var.x = 0; var.f = 0;

% update local phase equilibrium
var.c      = c;             % component fractions [wt]
var.T      = T;             % temperature [C]
var.P      = P/1e9;         % pressure [GPa]
var.H2O    = c(:,end);      % water concentration [wt]
cal.H2Osat = H2Osat/100;
[var,cal]  = meltmodel(var,cal,'E');

% get fitted phase oxide compositions
MLTfit = var.cm*cal.cmp_oxd;
SOLfit = var.cx*cal.cmp_oxd;

datafit = [var.m*100];%[MLTfit(:);SOLfit(:)];

end