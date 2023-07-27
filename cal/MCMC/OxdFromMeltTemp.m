function [datafit] = OxdFromMeltTemp(model,T,P,c,cal)

np = length(T);

cal.T0 = model(1:cal.ncmp-1).';
cal.r  = model(cal.ncmp:end).';
cal.A  = (cal.T0+273.15)./300;

var.m = 1.0; var.x = 0; var.f = 0;

% update local phase equilibrium
var.c      = c;        % component fractions [wt]
var.T      = T;          % temperature [C]
var.P      = P/1e9;      % pressure [GPa]
var.H2O    = c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P*1e9,0,cal);
[var,cal]  = meltmodel(var,cal,'E');

% get fitted phase oxide compositions
oxdLIQfit = var.cm*cal.cmp_oxd;
oxdSOLfit = var.cx*cal.cmp_oxd;

% get fitted phase fractions
phsfit = zeros(np,cal.nmsy+1);
phsfit(:,1) = var.m*100;

datafit = [oxdLIQfit(:);oxdSOLfit(:)];

end