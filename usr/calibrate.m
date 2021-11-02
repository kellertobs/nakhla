% calibrate phase diagram
addpath('../src');


% set parameters
cphs0    =  0.35;                % phase diagram lower bound composition [wt SiO2]
cphs1    =  0.70;                % phase diagram upper bound composition [wt SiO2]
Tphs0    =  750;                 % phase diagram lower bound temperature [degC]
Tphs1    =  1750;                % phase diagram upper bound temperature [degC]
PhDg     =  4.0;                 % Phase diagram curvature factor (> 1)
perCm    =  0.57;                % peritectic liquidus composition [wt SiO2]
perCx    =  0.52;                % peritectic solidus  composition [wt SiO2]
perT     =  1050;                % peritectic temperature [degC]
clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
dTH2O    =  1300;                % solidus shift from water content [degC/wt^0.75]

% plot T-dependence of melting
T = linspace(800,1500,1e4);  % temperature [degC]
P = 1e9.*ones(size(T));         % pressure [Pa]
c = 0.42.*ones(size(T));        % major component [wt SiO2]
v = 0.001.*ones(size(T));        % volatile component [wt H2O]

xq = ones(size(T)).*0.5;
fq = ones(size(T)).*0.0;

[xq,cxq,cmq,fq,vfq,vmq]  =  equilibrium(xq,fq,T,c,v,P,Tphs0,Tphs1,cphs0,cphs1,perT,perCx,perCm,clap,dTH2O,PhDg);