% Giordano08: a programme for calculating silicate melt viscosity after
%             Iacovino & Till (2019)
%
% Usage:     [rho] = Giordano08(oxd_wt,T,P)
%
% Input: 
%
%     oxd_wt Melt composition as 9-element vector containing concentrations 
%            in [wt%] of the following oxides ordered in the exact sequence 
%            [SiO2 TiO2 Al2O3 FeO MgO CaO Na2O K2O H2O]
%             1    2    3     4   5   6   7    8   9  
%     T      Melt emperature in [degC]
%
% Output:
%
%     eta    Melt viscosity in [Pas]
%
% Reference:
% Giordano D, Russell JK, & Dingwell DB (2008)
% Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
% Letters, 271, 123-134.
% doi: 10.1016/j.epsl.2008.03.038
%
% Modified by Tobias Keller (github.com/kellertobs), 03/2023
%
% Original Matlab code by JK Russel: www.eoas.ubc.ca/~krussell/ftp/jkrftp.html


function   [eta] = Giordano08(oxd_wt,T)

% Fitting parameter values
AT  =  -4.55;
bb  =  [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  =  [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

% Molar weights
MW  = [60.0855, 79.88, 101.96, 71.85, 40.3, 56.08, 61.98, 94.2, 18.02];

% Convert from wt % to mol%
% Note all data are normalized on anhydrous components first
norm_oxd_wt = [oxd_wt(:,1:8).*(100-oxd_wt(:,9))./(sum(oxd_wt(:,1:8),2)) oxd_wt(:,9)];
mp          = norm_oxd_wt./MW;
oxd_mol     = mp./sum(mp,2)*100;

% Load composition-basis matrix for multiplication against model-coefficients
siti  =  oxd_mol(:,1) + oxd_mol(:,2);
tial  =  oxd_mol(:,2) + oxd_mol(:,3);
fmg   =  oxd_mol(:,4) + oxd_mol(:,5);
nak   =  oxd_mol(:,7) + oxd_mol(:,8);
b1    =  siti;
b2    =  oxd_mol(:,3);
b3    =  oxd_mol(:,4);
b4    =  oxd_mol(:,5);
b5    =  oxd_mol(:,6);
b6    =  oxd_mol(:,7) + oxd_mol(:,9);
b7    =  oxd_mol(:,9) + log(1+oxd_mol(:,9));
b12   =  siti.*fmg;
b13   =  (siti + oxd_mol(:,3)).*(nak + oxd_mol(:,9));
b14   =  oxd_mol(:,3).*nak;

c1    =  oxd_mol(:,1);
c2    =  tial;
c3    =  fmg;
c4    =  oxd_mol(:,6);
c5    =  nak;
c6    =  log(1+oxd_mol(:,9));
c11   =  oxd_mol(:,3) + fmg + oxd_mol(:,6);
c11   =  c11.*(nak + oxd_mol(:,9));
bcf   =  [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf   =  [c1 c2 c3 c4 c5 c6 c11];

BT    =  sum(bb.*bcf,2);
CT    =  sum(cc.*ccf,2);

% Convert temperature to [K]
TK    =  max(600,T) + 273.15;

% Calculate melt viscosity in [Pas]
eta   =  10.^(min(12,max(-6,AT + BT./(TK(:)-CT))));

