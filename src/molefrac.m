function [nox,mf] = molefrac(wt)

% [nox,xmf] = MOLEFRAC(wt) Calculates mole fraction oxides from wt % 
% input : weight fraction oxides in order as SiO2 2 TiO2  3 Al2O3  4 FeO  5 MnO  6 MgO  7CaO 8 Na2O  9 K2O  10 P2O5 11 H2O 12 F2O-1  
% output: mole fractions of equivalent

% modified and optimised by Tobias Keller, 10. June, 2022

% molar weights
mw  = [ 60.0843, 79.8658, 101.961276, 71.8444, 70.937449,40.3044,56.0774, 61.97894, 94.1960, 141.9446,18.01528, 18.9984];

[~,nox] = size(wt);

% note all data are normalized on anhydrous components first
wtn = [wt(:,1:10).*(100-wt(:,11))./(sum(wt(:,1:10),2)+wt(:,12)) wt(:,11) 0.5.*wt(:,12).*(100-wt(:,11))./(sum(wt(:,1:10),2)+wt(:,12))];
mp  = wtn./mw;
mf  = 100.*(mp./sum(mp,2));

