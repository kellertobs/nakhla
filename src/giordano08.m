% MATLAB script to compute silicate melt viscosity.
%
%  Citation: Giordano D, Russell JK, & Dingwell DB (2008)
%  Viscosity of Magmatic Liquids: A Model. Earth & Planetary Science
%  Letters, v. 271, 123-134.
%
% ________________________________________________________________
% INPUT: Chemical compositions of silicate melts (wt% oxides) as:
% [SiO2 TiO2 Al2O3 FeO(T) MgO CaO Na2O K2O H2O]
%  1    2    3     4      5   6   7    8   9   
% _________________________________________________________________

% ________________________________________________________________
% OUTPUT: eta values at T(C) values (1 line per melt)
% ________________________________________________________________

% VFT Multicomponent-Model Coedfficients
% -----------------------------------------------------------------

% Modified and optimised by Tobias Keller (github.com/kellertobs), 06/2022

function   [eta] = giordano08(wt,TC)

AT  =  -4.55;
bb  =  [159.56  -173.34 72.13 75.69 -38.98 -84.08 141.54 -2.43 -0.91 17.62];
cc  =  [2.75 15.72 8.32 10.2 -12.29 -99.54 0.3 ];

% molar weights
mw  = [60.0855, 79.88, 101.96, 71.85, 40.3, 56.08, 61.98, 94.2, 18.02];

% convert from wt % to mol%
% note all data are normalized on anhydrous components first
% wtn = [wt(:,1:10).*(100-wt(:,11))./(sum(wt(:,1:10),2)+wt(:,12)) wt(:,11) 0.5.*wt(:,12).*(100-wt(:,11))./(sum(wt(:,1:10),2)+wt(:,12))];
wtn = [wt(:,1:8).*(100-wt(:,9))./(sum(wt(:,1:8),2)) wt(:,9)];
mp  = wtn./mw;
mol = mp./sum(mp,2)*100;

% Load composition-basis matrix for multiplication against model-coefficients
siti  =  mol(:,1) + mol(:,2);
tial  =  mol(:,2) + mol(:,3);
fmg   =  mol(:,4) + mol(:,5);
nak   =  mol(:,7) + mol(:,8);
b1    =  siti;
b2    =  mol(:,3);
b3    =  mol(:,4);
b4    =  mol(:,5);
b5    =  mol(:,6);
b6    =  mol(:,7) + mol(:,9);
b7    =  mol(:,9) + log(1+mol(:,9));
b12   =  siti.*fmg;
b13   =  (siti + mol(:,3)).*(nak + mol(:,9));
b14   =  mol(:,3).*nak;

c1    =  mol(:,1);
c2    =  tial;
c3    =  fmg;
c4    =  mol(:,6);
c5    =  nak;
c6    =  log(1+mol(:,9));
c11   =  mol(:,3) + fmg + mol(:,6);
c11   =  c11.*(nak + mol(:,9));
bcf   =  [b1 b2 b3 b4 b5 b6 b7 b12 b13 b14];
ccf   =  [c1 c2 c3 c4 c5 c6 c11];

BT    =  sum(bb.*bcf,2);
CT    =  sum(cc.*ccf,2);

% Convert temperature to [K]
TK    =  TC + 273.15;

% Calculate melt viscosity in [Pas]
eta   =  10.^(min(12,max(-6,AT + BT./(TK(:)-CT))));

