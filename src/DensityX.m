% DensityX: a programme for calculating silicate melt density after
%           Iacovino & Till (2019)
%
% Usage:     [rho] = DensityX(oxd_wt,T,P)
%
% Input: 
%
%     oxd_wt Melt composition as 9-element vector containing concentrations 
%            in [wt%] of the following oxides ordered in the exact sequence 
%            [SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O]
%             1    2    3     4   5   6   7    8   9  
%     T      Melt emperature in [degC]
%     P      Pressure in [kbar]
%
% Output:
%
%     rho    Melt density in [kg/m3]
%
% Reference:
% Iacovino K & Till C (2019) 
% DensityX: A program for calculating the densities of magmatic liquids up 
% to 1,627 °C and 30 kbar. Volcanica 2(1), p 1-10.
% doi: 10.30909/vol.02.01.0110
%
% Translated to Matlab by Tobias Keller (github.com/kellertobs), 03/2023
%
% Original Python code by K Iacovino: github.com/kaylai/DensityX


function [rho] = DensityX(oxd_wt,T,P)


% Set parameter values

% Molecular Weights [g/mol]
MW  = [60.0855, 79.88, 101.96, 71.85, 40.3, 56.08, 61.98, 94.2, 18.02];

% Partial Molar Volumes
% Volumes for SiO2, Al2O3, MgO, CaO, Na2O, K2O at Tref = 1773 K (Lange, 1997; CMP)
% Volume for H2O at Tref = 1273 K (Ochs and Lange, 1999)
% Volume for FeO at Tref = 1723 K (Guo et al., 2014)
% Volume for TiO2 at Tref = 1773 K (Lange and Carmichael, 1987)
MV  = [26.86, 28.32, 37.42, 12.68, 12.02, 16.90, 29.65, 47.28, 22.9];

% dV/dT values
% MgO, CaO, Na2O, K2O Table 4 (Lange, 1997)
% SiO2, TiO2, Al2O3 Table 9 (Lange and Carmichael, 1987)
% H2O from Ochs & Lange (1999)
% FeO from Guo et al (2014)
dVdT  = [0.0, 0.00724, 0.00262, 0.00369, 0.00327, 0.00374, 0.00768, 0.01208, 0.0095];

% dV/dP values
% Anhydrous component data from Kess and Carmichael (1991)
% H2O data from Ochs & Lange (1999)
dVdP  = [-0.000189, -0.000231, -0.000226, -0.000045, 0.000027, 0.000034, -0.00024, -0.000675, -0.00032];

% Tref values
Tref  = [1773.15, 1773.15, 1773.15, 1723.15, 1773.15, 1773.15, 1773.15, 1773.15, 1273.15];


% Prepare composition, temperature, and pressure variables

% Normalize original wt% values to 100% sum
norm_WP = oxd_wt./sum(oxd_wt,2).*100;

% Convert temperature to [K]
T_K     = T	+ 273.15;

% Convert pressure to [bar]
P_bar   = P*1e3;


% Calculate melt density

% Divide normalized wt% values by molecular weights
part_MP   = norm_WP ./ MW;
sum_MP    = sum(part_MP,2);

% convert to mol fraction
norm_MP   = part_MP ./ sum_MP;

% Calculate the partial Vliq for each oxide
part_Vliq = (MV  + (dVdT .* (T_K - Tref)) + (dVdP .* (P_bar-1))) .* norm_MP;
sum_Vliq  = sum(part_Vliq,2);

% Calculate partial X*MW
part_XMW  = norm_MP .* MW;
sum_XMW	  = sum(part_XMW,2);

% Calculate melt density in [kg/m3]
rho 	  = sum_XMW ./ sum_Vliq * 1000;
