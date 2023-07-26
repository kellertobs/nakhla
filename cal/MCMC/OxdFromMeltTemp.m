function [oxdfit] = OxdFromMeltTemp(mlttmp,T,P,c,cal)

Nz = length(T); Nx = 1;

cal.T0 = mlttmp(1:cal.ncmp-1).';
cal.r  = mlttmp(cal.ncmp:end).';
cal.A  = (cal.T0+273.15)./300;

var.m = 1.0; var.x = 0; var.f = 0;

% update local phase equilibrium
var.c      = c;        % component fractions [wt]
var.T      = T;          % temperature [C]
var.P      = P/1e9;      % pressure [GPa]
var.H2O    = c(:,end);      % water concentration [wt]
cal.H2Osat = fluidsat(var.T,var.P*1e9,0,cal);
[var,cal]  = meltmodel(var,cal,'E');

cm_oxd = var.cm*cal.cmp_oxd;
cx_oxd = var.cx*cal.cmp_oxd;

oxdfit = [cm_oxd(:);cx_oxd(:)];

end