% calibrate phase diagram
clear all; close all;

addpath(genpath('../'));
addpath('../../cal')
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
TINY = 1e-16;
FS = {'FontSize',15};
TX = {'Interpreter','latex'};

rng(15);    % for reproducibility

% set phase diagram parameters
cal_MORB;  % load melt model calibration

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

%       Si Ti Al Cr Fe Mg Ca Na K O H
ioxd = [1   2  3     5  6  7  8     11 ]; % oxide indices from MAGEMin to standard
Si = 1; Ti = 2; Al = 3; FeO = 4; Mg = 5; Ca = 6; Na = 7; H = 8;

%% load MAGEMin results

filename = 'MORB_fract10_out.mat';
load(filename);

% lump in free O to FeO, Cr2O3 to Al2O3, normalise to anhydrous unit sum
phs = fieldnames(OUT.PhaseProps);
phs = [phs(:)',{'SYS'},{'sol'},{'cum'}];
for iph = 1:length(phs)
    OUT.OxideFract.(phs{iph}) = zeros(size(OUT.OxideFractions.(phs{iph})));
    OUT.OxideFractions.(phs{iph})(:,4 ) = 0; % no Cr
    OUT.OxideFractions.(phs{iph})(:,9 ) = 0; % no K
    OUT.OxideFractions.(phs{iph})(:,10) = 0; % no O
    OUT.OxideFract.(phs{iph}) = OUT.OxideFractions.(phs{iph})(:,ioxd);
    OUT.OxideFract.(phs{iph}) = OUT.OxideFract.(phs{iph})./sum(OUT.OxideFract.(phs{iph})+1e-16,2).*100;
end

% combine all feldspar instances
if isfield(OUT.PhaseProps,'pl4T2')
    OUT.OxideFract.pl4T = (OUT.OxideFract.pl4T.*OUT.PhaseFractions.pl4T + OUT.OxideFract.pl4T2.*OUT.PhaseFractions.pl4T2) ./ (OUT.PhaseFractions.pl4T+OUT.PhaseFractions.pl4T2+1e-16);
    OUT.EMFractions.pl4T = (OUT.EMFractions.pl4T.*OUT.PhaseFractions.pl4T + OUT.EMFractions.pl4T2.*OUT.PhaseFractions.pl4T2) ./ (OUT.PhaseFractions.pl4T+OUT.PhaseFractions.pl4T2+1e-16);
    OUT.PhaseFractions.pl4T = OUT.PhaseFractions.pl4T+OUT.PhaseFractions.pl4T2;
    OUT.OxideFract = rmfield(OUT.OxideFract,'pl4T2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'pl4T2');
end

% combine all orthopyroxene instances
if isfield(OUT.PhaseFractions,'opx2')
    OUT.OxideFract.olv = (OUT.OxideFract.olv.*OUT.PhaseFractions.olv + OUT.OxideFract.opx2.*OUT.PhaseFractions.opx2) ./ (OUT.PhaseFractions.olv+OUT.PhaseFractions.opx2+1e-16);
    OUT.EMFractions.olv = (OUT.EMFractions.olv.*OUT.PhaseFractions.olv + OUT.EMFractions.opx2.*OUT.PhaseFractions.opx2) ./ (OUT.PhaseFractions.olv+OUT.PhaseFractions.opx2+1e-16);
    OUT.PhaseFractions.olv = OUT.PhaseFractions.olv+OUT.PhaseFractions.opx2;
    OUT.OxideFract = rmfield(OUT.OxideFract,'opx2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'opx2');
end

% combine all clinopyroxene instances
if isfield(OUT.PhaseFractions,'cpx2')
    OUT.OxideFract.cpx = (OUT.OxideFract.cpx.*OUT.PhaseFractions.cpx + OUT.OxideFract.cpx2.*OUT.PhaseFractions.cpx2) ./ (OUT.PhaseFractions.cpx+OUT.PhaseFractions.cpx2+1e-16);
    OUT.EMFractions.cpx = (OUT.EMFractions.cpx.*OUT.PhaseFractions.cpx + OUT.EMFractions.cpx2.*OUT.PhaseFractions.cpx2) ./ (OUT.PhaseFractions.cpx+OUT.PhaseFractions.cpx2+1e-16);
    OUT.PhaseFractions.cpx = OUT.PhaseFractions.cpx+OUT.PhaseFractions.cpx2;
    OUT.OxideFract = rmfield(OUT.OxideFract,'cpx2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'cpx2');
end

% % combine all pyroxene instances
% if isfield(OUT.PhaseFractions,'opx')
%     OUT.OxideFract.pxn = (OUT.OxideFract.cpx.*OUT.PhaseFractions.cpx + OUT.OxideFract.opx.*OUT.PhaseFractions.opx) ./ (OUT.PhaseFractions.cpx+OUT.PhaseFractions.opx+1e-16);
%     OUT.PhaseFractions.pxn = OUT.PhaseFractions.cpx+OUT.PhaseFractions.opx;
%     OUT.OxideFract = rmfield(OUT.OxideFract,'opx');
%     OUT.PhaseFractions = rmfield(OUT.PhaseFractions,'opx');
% end

% combine all spinel instances
if isfield(OUT.PhaseFractions,'spn2')
    OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseFractions.spn + OUT.OxideFract.spn2.*OUT.PhaseFractions.spn2) ./ (OUT.PhaseFractions.spn+OUT.PhaseFractions.spn2+1e-16);
    OUT.EMFractions.spn = (OUT.EMFractions.spn.*OUT.PhaseFractions.spn + OUT.EMFractions.spn2.*OUT.PhaseFractions.spn2) ./ (OUT.PhaseFractions.spn+OUT.PhaseFractions.spn2+1e-16);
    OUT.PhaseFractions.spn = OUT.PhaseFractions.spn+OUT.PhaseFractions.spn2;
    OUT.OxideFract = rmfield(OUT.OxideFract,'spn2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'spn2');
end

% % lump in ilmenite with spinel
% if isfield(OUT.PhaseFractions,'ilm')
%     OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseFractions.spn + OUT.OxideFract.ilm.*OUT.PhaseFractions.ilm) ./ (OUT.PhaseFractions.spn+OUT.PhaseFractions.ilm+1e-16);
%     OUT.PhaseFractions.spn = OUT.PhaseFractions.spn+OUT.PhaseFractions.ilm;
%     % OUT.OxideFract = rmfield(OUT.OxideFract,'ilm');
%     % OUT.PhaseFractions = rmfield(OUT.PhaseFractions,'ilm');
% end


%% collate and reduce mineral and melt oxide compositions

% detect where phases are stable
hasMLT = OUT.PhaseFractions.liq >=1e-4 & OUT.PhaseFractions.sol_wt>=1e-4 & OUT.PhaseFractions.q<=1e-4;
hasCUM = OUT.PhaseFractions.cum >=1e-4 & hasMLT;
hasOLV = OUT.PhaseFractions.ol  >=1e-4 & hasMLT;
hasFSP = OUT.PhaseFractions.pl4T>=1e-4 & hasMLT;
hasCPX = OUT.PhaseFractions.cpx >=1e-4 & hasMLT;
hasOPX = OUT.PhaseFractions.opx >=1e-4 & hasMLT;
hasSPN = OUT.PhaseFractions.spn >=1e-4 & hasMLT;
hasILM = OUT.PhaseFractions.ilm >=1e-4 & hasMLT;
hasQTZ = OUT.PhaseFractions.q   >=1e-4 & hasMLT;

% set oxides present in phases
oxdSYS = [Si,Ti,Al,FeO,Mg,Ca,Na,H]; noxd = length(oxdSYS);
oxdMLT = [Si,Ti,Al,FeO,Mg,Ca,Na,H];
oxdCUM = [Si,Ti,Al,FeO,Mg,Ca,Na,H]; 
oxdOLV = [Si,      FeO,Mg        ];
oxdFSP = [Si,   Al,       Ca,Na  ];
oxdCPX = [Si,   Al,FeO,Mg,Ca,Na  ];
oxdOPX = [Si,   Al,FeO,Mg,Ca     ];
oxdSPN = [   Ti,Al,FeO,Mg        ];
oxdILM = [   Ti,   FeO           ];
oxdQTZ = [Si                     ];

% extract oxide composition of phases
SOL = OUT.OxideFract.sol (hasMLT,oxdMLT).*100; SOL = SOL./sum(SOL,2)*100;
MLT = OUT.OxideFract.liq (hasMLT,oxdMLT).*100; MLT = MLT./sum(MLT,2)*100; nMLT = size(MLT,1);
SYS = (SOL.*OUT.PhaseFractions.sol_wt(hasMLT) + MLT.*OUT.PhaseFractions.liq_wt(hasMLT)) ...
    ./(     OUT.PhaseFractions.sol_wt(hasMLT) +      OUT.PhaseFractions.liq_wt(hasMLT));
CUM = zeros(length(hasMLT(hasCUM)),length(oxdMLT)); CUM(:,oxdCUM) = OUT.OxideFract.cum (hasCUM,oxdCUM); CUM = CUM./sum(CUM,2)*100; nCUM = size(CUM,1);
OLV = zeros(length(hasMLT(hasOLV)),length(oxdMLT)); OLV(:,oxdOLV) = OUT.OxideFract.ol  (hasOLV,oxdOLV); OLV = OLV./sum(OLV,2)*100; nOLV = size(OLV,1);
FSP = zeros(length(hasFSP(hasFSP)),length(oxdMLT)); FSP(:,oxdFSP) = OUT.OxideFract.pl4T(hasFSP,oxdFSP); FSP = FSP./sum(FSP,2)*100; nFSP = size(FSP,1);
CPX = zeros(length(hasMLT(hasCPX)),length(oxdMLT)); CPX(:,oxdCPX) = OUT.OxideFract.cpx (hasCPX,oxdCPX); CPX = CPX./sum(CPX,2)*100; nCPX = size(CPX,1);
OPX = zeros(length(hasMLT(hasOPX)),length(oxdMLT)); OPX(:,oxdOPX) = OUT.OxideFract.opx (hasOPX,oxdOPX); OPX = OPX./sum(OPX,2)*100; nOPX = size(OPX,1);
SPN = zeros(length(hasMLT(hasSPN)),length(oxdMLT)); SPN(:,oxdSPN) = OUT.OxideFract.spn (hasSPN,oxdSPN); SPN = SPN./sum(SPN,2)*100; nSPN = size(SPN,1);
ILM = zeros(length(hasMLT(hasILM)),length(oxdMLT)); ILM(:,oxdILM) = OUT.OxideFract.ilm (hasILM,oxdILM); ILM = ILM./sum(ILM,2)*100; nILM = size(ILM,1);
QTZ = zeros(length(hasQTZ(hasQTZ)),length(oxdMLT)); QTZ(:,oxdQTZ) = OUT.OxideFract.q   (hasQTZ,oxdQTZ); QTZ = QTZ./sum(QTZ,2)*100; nQTZ = size(QTZ,1);
T   = OUT.T(hasMLT);
P   = OUT.P(hasMLT)*1e8;
H2O = OUT.OxideFract.liq(hasMLT,H)*100;

PHS    = [OUT.PhaseFractions.liq( hasMLT  )./(OUT.PhaseFractions.sol_wt(hasMLT)+1*OUT.PhaseFractions.liq_wt(hasMLT)), ...
          OUT.PhaseFractions.ol(  hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT)), ...
          OUT.PhaseFractions.pl4T(hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT)), ...
          OUT.PhaseFractions.cpx( hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT)), ...
          OUT.PhaseFractions.spn( hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT))+ ...
          OUT.PhaseFractions.ilm( hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT)), ...
          OUT.PhaseFractions.opx( hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT)), ...
          OUT.PhaseFractions.q(   hasMLT,1)./(OUT.PhaseFractions.sol_wt(hasMLT)+0*OUT.PhaseFractions.liq_wt(hasMLT))];

Tsol = 925.*ones(size(T)); % solidus estimate from MAGEMin
Tliq = 1335.*ones(size(T)); % liquidus estimate from MAGEMin
Psl  = 0.2.*ones(size(T)); % P [GPa]

%% olivine system
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdOLV);
DATA.SNAMES = {};
DATA.X = OLV(:,oxdOLV);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
OLV_PCA = DGN;

OLVp           = zeros(size(OLV));
OLVp(:,oxdOLV) = max(0,Xp)./sum(max(0,Xp),2)*100;
OLV_PCA.OLVp   = OLVp;

OLV_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
OLV_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%% plot olivine system
cal_MORB; % load melt model calibration

figure(102); clf;
subplot(1,2,1);
scatter(OLV (:,cal.Si),OLV (:,cal.Fe),25,T(hasOLV(hasMLT))); colormap('copper'); hold on
scatter(OLVp(:,cal.Si),OLVp(:,cal.Fe),25,T(hasOLV(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(1,2,2);
scatter(OLV (:,cal.Si),OLV (:,cal.Mg),25,T(hasOLV(hasMLT))); colormap('copper'); hold on
scatter(OLVp(:,cal.Si),OLVp(:,cal.Mg),25,T(hasOLV(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})

sgtitle('OLV PCA',FS{:},TX{:})
drawnow


%% feldspar system
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdFSP);
DATA.SNAMES = {};
DATA.X = FSP(:,oxdFSP);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
FSP_PCA = DGN;

FSPp           = zeros(size(FSP));
FSPp(:,oxdFSP) = max(0,Xp)./sum(max(0,Xp),2)*100;
FSP_PCA.FSPp   = FSPp;

FSP_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
FSP_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%% plot feldspar system
cal_MORB; % load melt model calibration

figure(105); clf;
subplot(1,3,1);
scatter(FSP (:,cal.Si),FSP (:,cal.Al),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.Al),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(1,3,2);
scatter(FSP (:,cal.Si),FSP (:,cal.Ca),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.Ca),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(1,3,3);
scatter(FSP (:,cal.Si),FSP (:,cal.Na),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.Na),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})

sgtitle('FSP PCA',FS{:},TX{:})
drawnow


%% clinopyroxene system
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdCPX);
DATA.SNAMES = {};
DATA.X = CPX(:,oxdCPX);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
CPX_PCA = DGN;

CPXp           = zeros(size(CPX));
CPXp(:,oxdCPX) = max(0,Xp)./sum(max(0,Xp),2)*100;
CPX_PCA.CPXp   = CPXp;

CPX_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
CPX_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%%
cal_MORB; % load melt model calibration

figure(103); clf;
subplot(2,3,1);
scatter(CPX (:,cal.Si),CPX (:,cal.Al),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Al),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,2);
scatter(CPX (:,cal.Si),CPX (:,cal.Fe),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Fe),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,3);
scatter(CPX (:,cal.Si),CPX (:,cal.Mg),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Mg),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,4);
scatter(CPX (:,cal.Si),CPX (:,cal.Ca),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Ca),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,5);
scatter(CPX (:,cal.Si),CPX (:,cal.Na),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Na),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
sgtitle('CPX PCA',FS{:},TX{:})
drawnow


%% oxides system
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdSPN);
DATA.SNAMES = {};
DATA.X = SPN(:,oxdSPN);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
SPN_PCA = DGN;

SPNp           = zeros(size(SPN));
SPNp(:,oxdSPN) = max(0,Xp)./sum(max(0,Xp),2)*100;
SPN_PCA.SPNp   = SPNp;

SPN_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
SPN_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);

DATA.X = ILM(:,oxdILM);
DATA.X = DATA.X./sum(DATA.X,2);

Fe             = mean(DATA.X);
Fi             = Fe;
Xp             = repmat(Fe,nILM,1);

ILMp           = zeros(size(ILM));
ILMp(:,oxdILM) = max(0,Xp)./sum(max(0,Xp),2)*100;
ILM_PCA.ILMp   = ILMp;

ILM_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
ILM_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);

%%
cal_MORB; % load melt model calibration

figure(104); clf;
subplot(1,2,1);
scatter(SPN (:,cal.Fe),SPN (:,cal.Ti),25,T(hasSPN(hasMLT))); colormap('copper'); hold on
scatter(SPNp(:,cal.Fe),SPNp(:,cal.Ti),25,T(hasSPN(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tms,cal.Fe),cal.mem_oxd(cal.tms,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.ilm,cal.Fe),cal.mem_oxd(cal.ilm,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(1,2,2);
scatter(SPN (:,cal.Fe),SPN (:,cal.Mg),25,T(hasSPN(hasMLT))); colormap('copper'); hold on
scatter(SPNp(:,cal.Fe),SPNp(:,cal.Mg),25,T(hasSPN(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tms,cal.Fe),cal.mem_oxd(cal.tms,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.ilm,cal.Fe),cal.mem_oxd(cal.ilm,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
% colorbar;
sgtitle('SPN PCA',FS{:},TX{:})
drawnow

%% orthopyroxene system
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdOPX);
DATA.SNAMES = {};
DATA.X = OPX(:,oxdOPX);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
OPX_PCA = DGN;

OPXp           = zeros(size(OPX));
OPXp(:,oxdOPX) = max(0,Xp)./sum(max(0,Xp),2)*100;
OPX_PCA.OPXp   = OPXp;

OPX_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
OPX_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%%
cal_MORB; % load melt model calibration

figure(103); clf;
subplot(2,2,1);
scatter(OPX (:,cal.Si),OPX (:,cal.Al),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Al),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(OPX (:,cal.Si),OPX (:,cal.Fe),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Fe),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,2,3);
scatter(OPX (:,cal.Si),OPX (:,cal.Mg),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Mg),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,2,4);
scatter(OPX (:,cal.Si),OPX (:,cal.Ca),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Ca),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
sgtitle('OPX PCA',FS{:},TX{:})
drawnow



%% add up projected mineral compositions to solid composition
wt   = zeros(size(T)) + 1e-16;
SOLp = zeros(size(SOL));

SOLp(hasOLV(hasMLT),:) = SOLp(hasOLV(hasMLT),:) + OLVp.*OUT.PhaseFractions.ol(hasOLV);
wt(hasOLV(hasMLT))     = wt(hasOLV(hasMLT)) + OUT.PhaseFractions.ol(hasOLV);

SOLp(hasFSP(hasMLT),:) = SOLp(hasFSP(hasMLT),:) + FSPp.*OUT.PhaseFractions.pl4T(hasFSP);
wt(hasFSP(hasMLT))     = wt(hasFSP(hasMLT)) + OUT.PhaseFractions.pl4T(hasFSP);

SOLp(hasCPX(hasMLT),:) = SOLp(hasCPX(hasMLT),:) + CPXp.*OUT.PhaseFractions.cpx(hasCPX);
wt(hasCPX(hasMLT))     = wt(hasCPX(hasMLT)) + OUT.PhaseFractions.cpx(hasCPX);

SOLp(hasSPN(hasMLT),:) = SOLp(hasSPN(hasMLT),:) + SPNp.*OUT.PhaseFractions.spn(hasSPN);
wt(hasSPN(hasMLT))     = wt(hasSPN(hasMLT)) + OUT.PhaseFractions.spn(hasSPN);

SOLp(hasILM(hasMLT),:) = SOLp(hasILM(hasMLT),:) + ILMp.*OUT.PhaseFractions.ilm(hasILM);
wt(hasILM(hasMLT))     = wt(hasILM(hasMLT)) + OUT.PhaseFractions.ilm(hasILM);

SOLp(hasOPX(hasMLT),:) = SOLp(hasOPX(hasMLT),:) + OPXp.*OUT.PhaseFractions.opx(hasOPX);
wt(hasOPX(hasMLT))     = wt(hasOPX(hasMLT)) + OUT.PhaseFractions.opx(hasOPX);

SOLp(hasQTZ(hasMLT),:) = SOLp(hasQTZ(hasMLT),:) + QTZ.*OUT.PhaseFractions.q(hasQTZ);
wt(hasQTZ(hasMLT))     = wt(hasQTZ(hasMLT)) + OUT.PhaseFractions.q(hasQTZ);

SOLp = SOLp./wt;

% collate oxide compositions into global data array
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdMLT(1:H));
DATA.SNAMES = {};
DATA.X      = [MLT(:,1:end-1);SOLp(:,1:end-1)];%;CUM(:,1:end-1)];
DATA.X      = DATA.X./sum(DATA.X,2);

unmix;
MLT_PCA = DGN;

MLTp = MLT;  CUMp = CUM;
MLTp(:,1:end-1) = max(0,Xp(0   +(1:nMLT),:))./sum(max(0,Xp(0   +(1:nMLT),:)),2).*(100-MLT(:,end));

EMExt = max(0,Fe)./sum(max(0,Fe),2)*100;
EMInt = max(0,Fi)./sum(max(0,Fi),2)*100;

[~,is] = sort(EMExt(:,Si)-EMExt(:,Mg),'ascend');
EMExt = EMExt(is,:);

[~,is] = sort(EMInt(:,Si)-EMInt(:,Mg),'ascend');
EMInt = EMInt(is,:);

MLT_PCA.EMInt = EMInt;
MLT_PCA.EMExt = EMExt;


%% reconstruct bulk composition based on projected solid and liquid compositions
SYSp = (SOLp.*OUT.PhaseFractions.sol_wt(hasMLT) + MLTp.*OUT.PhaseFractions.liq_wt(hasMLT)) ...
     ./(      OUT.PhaseFractions.sol_wt(hasMLT) +       OUT.PhaseFractions.liq_wt(hasMLT));

close all;
save('MAGEMin_processed');


%% liquid, solid, mixture compositions
cal_MORB;  % load melt model calibration
figure(107); clf;
subplot(2,3,1);
scatter(MLT(:,Si),MLT(:,Al),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Al),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Al),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Al),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Al),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Al),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Al),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Al),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,2);
scatter(MLT(:,Si),MLT(:,FeO),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,FeO),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,FeO),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,FeO),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,FeO),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,FeO),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,FeO),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,FeO),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,3);
scatter(MLT(:,Si),MLT(:,Mg),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Mg),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Mg),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Mg),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Mg),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Mg),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Mg),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Mg),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,4);
scatter(MLT(:,Si),MLT(:,Ca),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Ca),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Ca),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Ca),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Ca),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Ca),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Ca),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Ca),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,5);
scatter(MLT(:,Si),MLT(:,Na),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Na),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Na),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Na),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Na),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Na),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Na),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Na),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
subplot(2,3,6);
scatter(MLT(:,Si),MLT(:,H),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,H),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,H),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,H),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,H),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,H),25,OUT.T(hasMLT),'d','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.H),FS{:},TX{:})
% colorbar;

sgtitle('MLT \& SOL PCA',FS{:},TX{:})
drawnow


%% load projected data and prepare for fitting routines
load('MAGEMin_processed');
cal_MORB;  % load melt model calibration
                % for fay ant alb dps aug tms mgt ilm hyp fsl qtz wat
indmem  = logical([1   1   0   0   0   0   0   0   0   0   0   0   0
                   1   1   1   0   0   0   0   0   0   0   0   0   0
                   0   1   1   1   1   0   1   0   0   0   0   0   0
                   0   0   1   1   1   1   1   1   0   1   0   0   0
                   0   0   1   1   0   1   0   1   1   1   1   0   0
                   0   0   0   1   0   0   0   0   1   0   1   1   0
                   0   0   0   0   0   0   0   0   0   0   0   0   1]);


% convert factor analysis end-member to mineral end-member proportions

cmp_oxd = MLT_PCA.EMInt;%(MLT_PCA.EMInt+MLT_PCA.EMExt)/2;
cmp_oxd = cmp_oxd./sum(cmp_oxd,2)*100;

Xp = zeros(size(cmp_oxd,1),cal.nmem);
for ip = 1:size(cmp_oxd,1)
    Xp(ip,indmem(ip,:)) = lsqnonneg(cal.mem_oxd(indmem(ip,:),oxdSYS(1:end-1)).',cmp_oxd(ip,:).');
end
cmp_mem = Xp./sum(Xp,2)*100;

cmp_mem_FINT = zeros(cal.ncmp,cal.nmem);
cmp_mem_FINT(1:end-1,:) = cmp_mem;
cmp_mem_FINT(end,cal.wat) = 100;
cmp_oxd_FINT = cmp_mem_FINT*cal.mem_oxd/100;

cmp_mem_MAP = cmp_mem_FINT;

T0_MAP = [1880  1210  1165  1090  970  820];
A_MAP  = (cal.T0+273.15)./350;
B_MAP  = [9  5  4.3  3.5  3  2.5];
r_MAP  = [24.5  2.5  3.5  10.0  14.0  10.0];


%%
cal_MORB;  % load melt model calibration

data   = [MLTp(:);SOLp(:);PHS(:);Tsol(:);Tliq(:)];

m0     = [T0_MAP.';A_MAP.';B_MAP.';r_MAP.';cmp_mem_MAP(:).*indmem(:);];
m0_lw  = m0 - [max(10,0.2*T0_MAP.');0*max(0.1,0.01*A_MAP.');0*max(0.5,0.2*B_MAP.');max(1,0.5*r_MAP.');max(100,1*cmp_mem_MAP(:)).*indmem(:)];
m0_up  = m0 + [max(10,0.2*T0_MAP.');0*max(0.1,0.01*A_MAP.');0*max(0.5,0.2*B_MAP.');max(1,0.5*r_MAP.');max(100,1*cmp_mem_MAP(:)).*indmem(:)];
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(m0==100 ,:) = 100;
mbnds(m0==T0_MAP(1),:) = T0_MAP(1);
mbnds(   cal.ncmp-1 +(1:cal.ncmp-1)       ,:) = max(2,        mbnds(   cal.ncmp-1 +(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(0,min(100,mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));

mNames = cell(cal.ncmp*cal.nmem+4*(cal.ncmp-1),1);
k = 1;
for j=1:4
    for i=1:cal.ncmp-1
        if j==1
            mNames{k} = ['T0:',cal.cmpStr{i}];
        elseif j==2
            mNames{k} = ['A:',cal.cmpStr{i}];
        elseif j==3
            mNames{k} = ['B:',cal.cmpStr{i}];
        elseif j==4
            mNames{k} = ['r:',cal.cmpStr{i}];
        end
        k = k+1;
    end
end
for j=1:cal.nmem
    for i=1:cal.ncmp
        mNames{k} = [cal.memStr{j},':',cal.cmpStr{i}];
        k = k+1;
    end
end

% set data uncertainties
sigma_MLTSOL    =  0.1 * ones(size([MLTp(:);SOLp(:)]));
sigma_PHS       =  1.0 * ones(size([PHS(:)]));
sigma_TsolTliq  =  1.0 * ones(size([Tsol(:);Tliq(:)]));
sigma = [sigma_MLTSOL;sigma_PHS;sigma_TsolTliq];

% function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,T,P,MLTp,SOLp,SYSp,PHS(:,1),Psl,cal);

% function to apply further constraints to a proposed set of model param values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% function to calculate likelihood of dhat
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.1,10,model,cal);

% run MCMC algorithm
Niter = 2e5;

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.0015 * diff(mbnds,1,2);
anneal.levels   = 3;
anneal.burnin   = max(1,Niter/5);
anneal.refine   = max(1,Niter/10);
bestfit         = m0;

tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),m0,mbnds,anneal,Niter);
RunTime(1) = toc;

plotmcmc(models, prob, [], mbnds, anneal, mNames); 

T0_MAP       = bestfit(               (1:cal.ncmp-1)).';
A_MAP        = bestfit(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
B_MAP        = bestfit(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
r_MAP        = bestfit(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem_MAP  = reshape(bestfit(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd_MAP  = cmp_mem_MAP*cal.mem_oxd/100;

[dhat,MLTfit,SOLfit,SYSfit,PHSfit,cmpSYS,Tsolfit,Tliqfit,Tm] = dhatFunc([T0_MAP.';A_MAP.';B_MAP.';r_MAP.';cmp_mem_MAP(:)]);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

%
if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(1e5,3e9,50).';
var.m      = ones(size(PP)); var.x = 0*var.m; var.f = 0*var.m;
cal.T0     = T0_MAP;
cal.A      = (cal.T0+273.15)./350;
cal.r      = r_MAP;
var.c      = repmat(mean(cmpSYS(1:5,:).*[1,1,1,1,1,1,0]./sum(cmpSYS(1:5,1:end-1),2),1),length(PP),1);   % component fractions [wt]
var.P      = PP/1e9;         % pressure [GPa]
var.T      = 1000+PP*5e-8;             % temperature [C]
var.H2O    = zeros(size(PP)); % water concentration [wt]
cal.H2Osat = var.H2O+0.001;
[~,cal]    = meltmodel(var,cal,'T');

Tm         = cal.Tm;

% retrieve distributions
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% liquid, solid, mixture compositions
% cal_MORB; % load melt model calibration
 
figure(108); clf;

subplot(2,4,1);
scatter(MLTp(:,Si),MLTp(:,Ti),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Ti),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Ti),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ti),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ti),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ti),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})

subplot(2,4,2);
scatter(MLTp(:,Si),MLTp(:,Al),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Al),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Al),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Al),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Al),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Al),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})

subplot(2,4,3);
scatter(MLTp(:,Si),MLTp(:,FeO),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,FeO),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,FeO),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,FeO),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,FeO),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,FeO),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})

subplot(2,4,4);
scatter(MLTp(:,Si),MLTp(:,Mg),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Mg),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Mg),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Mg),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Mg),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Mg),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})

subplot(2,4,5);
scatter(MLTp(:,Si),MLTp(:,Ca),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Ca),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Ca),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ca),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ca),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ca),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

subplot(2,4,6);
scatter(MLTp(:,Si),MLTp(:,Na),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Na),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Na),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Na),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Na),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Na),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})

subplot(2,4,7);
scatter(MLTp(:,Si),MLTp(:,H ),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,H ),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,H ),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,H ),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,H ),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,H ),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.H ),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.H ),FS{:},TX{:})

sgtitle('MCMC component fit',FS{:},TX{:})
drawnow

figure(109); clf; cmap = colororder;
plot(T,PHS(:,1),'-','Color',cmap(2,:),'LineWidth',1.5); axis tight; hold on % melt
plot(T,PHS(:,2),'-','Color',cmap(5,:),'LineWidth',1.5); % olv
plot(T,PHS(:,3),'-','Color',cmap(6,:),'LineWidth',1.5); % opx
plot(T,PHS(:,4),'-','Color',cmap(3,:),'LineWidth',1.5); % cpx
plot(T,PHS(:,5),'-','Color',cmap(4,:),'LineWidth',1.5); % fsp
plot(T,PHS(:,6),'-','Color',cmap(1,:),'LineWidth',1.5); % spn
plot(T,PHS(:,7),'-','Color',cmap(7,:),'LineWidth',1.5); % qtz

plot(T,PHSfit(:,1),'--','Color',cmap(2,:),'LineWidth',1.5); % mlt
plot(T,PHSfit(:,2),'--','Color',cmap(5,:),'LineWidth',1.5); % olv
plot(T,PHSfit(:,3),'--','Color',cmap(6,:),'LineWidth',1.5); % opx
plot(T,PHSfit(:,4),'--','Color',cmap(3,:),'LineWidth',1.5); % cpx
plot(T,PHSfit(:,5),'--','Color',cmap(4,:),'LineWidth',1.5); % fsp
plot(T,PHSfit(:,6),'--','Color',cmap(1,:),'LineWidth',1.5); % spn
plot(T,PHSfit(:,7),'--','Color',cmap(7,:),'LineWidth',1.5); % qtz
legend(['mlt',cal.msyStr],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]',FS{:},TX{:})
ylabel('Phase proportions [wt\%]',FS{:},TX{:})
sgtitle('Phase stability',FS{:},TX{:})
drawnow

% plot phase diagram
figure(110); clf;
subplot(2,4,1)
plot(MLTp(:,Si)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Si)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Si)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Si)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Si)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,Si)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{Si},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,2)
plot(MLTp(:,Ti)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Ti)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Ti)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Ti)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Ti)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,Ti)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{Ti},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,3)
plot(MLTp(:,Al)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Al)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Al)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Al)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Al)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,Al)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{Al},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,4)
plot(MLTp(:,FeO)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,FeO)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,FeO)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,FeO)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,FeO)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,FeO)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{FeO},' [wt]'],'Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

subplot(2,4,5)
plot(MLTp(:,Mg)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Mg)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Mg)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Mg)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Mg)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,Mg)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{Mg},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,6)
plot(MLTp(:,Ca)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Ca)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Ca)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Ca)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Ca)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,Ca)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{Ca},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,7)
plot(MLTp(:,Na)./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Na)./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Na)./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Na)./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Na)./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,Na)./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{Na},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,8)
plot(MLTp(:,H )./sum(MLTp(:,1:end-1),2)*100,T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,H )./sum(SOLp(:,1:end-1),2)*100,T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,H )./sum(SYSp(:,1:end-1),2)*100,T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,H )./sum(MLTfit(:,1:end-1),2)*100,T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,H )./sum(SOLfit(:,1:end-1),2)*100,T,'bs');
plot(SYSfit(:,H )./sum(SYSfit(:,1:end-1),2)*100,T,'kd');

xlabel([cal.oxdStr{H },' [wt]'],'Interpreter','latex','FontSize',15)
drawnow

figure(111); clf
plot(Tm,PP/1e9,'LineWidth',1); axis ij tight; hold on
plot(Tsolfit,Psl,'bo','LineWidth',2);
plot(Tliqfit,Psl,'ro','LineWidth',2);
plot(Tsol,Psl,'ko','LineWidth',1.5);
plot(Tliq,Psl,'ko','LineWidth',1.5);

legend([cal.cmpStr(1:end-1),'Tsol','Tliq'],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Pressure [GPa]','Interpreter','latex','FontSize',15)
drawnow

%%
save('MORB_calibration');


%% update material closures
Nz = length(T); Nx = 1; Ptop = min(P); Pt = P; calibrt = 1;
var.m = ones(size(T)); var.x = 0*var.m; var.f = 0*var.m;

var.c      = c;             % component fractions [wt]
var.T      = T;             % temperature [C]
var.P      = P/1e9;         % pressure [GPa]
var.H2O    = c(:,end);      % water concentration [wt]
cal.H2Osat = MLTfit(:,H)/100;
[var,cal]  = meltmodel(var,cal,'E');

c  = zeros(Nz,Nx,cal.ncmp); c (:,1,:) = var.c;
cm = zeros(Nz,Nx,cal.ncmp); cm(:,1,:) = var.cm;
cx = zeros(Nz,Nx,cal.ncmp); cx(:,1,:) = var.cx;

m = var.m; 
x = var.x; 
f = var.f;

T = T+273.15;

update;

T = T-273.15;

wm =  mu.*Ksgr_m .* (rhom-rho)*g0; % melt segregation speed
wx = chi.*Ksgr_x .* (rhox-rho)*g0; % crystal segregation speed
wf = phi.*Ksgr_f .* (rhof-rho)*g0; % fluid segregation speed

%% plot phase fractions
figure(8); clf;
plot(T,x.*100,'k',T,m.*100,'r',T,f.*1000,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',15,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Melting model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Phase fractions [wt\%]','Interpreter','latex','FontSize',15)

% plot major phase compositions
figure(9); clf;
plot(T,cx.*100,'b',T,cm.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)

% plot phase densities
figure(11); clf;
plot(T,rhox,'k',T,rhom,'r',T,rhof,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
plot(T,rho ,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','mixture','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Density model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Density [kg/m$^3$]','Interpreter','latex','FontSize',15)

% plot mixture rheology
figure(12); clf;
semilogy(T,eta,'k',T,etam,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('mixture','melt','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Viscosity model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Viscosity [log$_{10}$ Pas]','Interpreter','latex','FontSize',15)

% plot phase segregation speeds
figure(13); clf;
semilogy(T,max(1e-18,abs(chi.*wx)).*3600,'k',T,max(1e-18,abs(mu.*wm)).*3600,'r',T,max(1e-18,abs(phi.*wf)).*3600,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase segregation model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Segregation flux [m/hr]','Interpreter','latex','FontSize',15)

% plot oxide compositions
figure(14); clf;
subplot(2,1,1)
sgtitle('Phase Oxide Fractions','Interpreter','latex','FontSize',18)
for i=1:cal.noxd
    plot(T,cm_oxd(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.noxd)+1,:)); hold on; box on; axis tight;
end
legend(cal.oxdStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
for i=1:cal.noxd
    plot(T,cx_oxd(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.noxd)+1,:)); hold on; box on; axis tight;
end
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

% plot end-member component compositions
figure(15); clf;
subplot(2,1,1)
sgtitle('Phase Component Fractions','Interpreter','latex','FontSize',18)
for i=1:cal.ncmp
    plot(T,cm(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); hold on; box on; axis tight;
end
legend(cal.cmpStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
for i=1:cal.ncmp
    plot(T,cx(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); hold on; box on; axis tight;
end
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

% plot simplified mineral assemblage
figure(16); clf;
patch([T;flipud(T)],[zeros(size(T));flipud(cx_msy(:,1))],[0.6,0.8,0.5],'LineWidth',2); hold on; box on; axis tight;
patch([T;flipud(T)],[sum(cx_msy(:,1  ),2);flipud(sum(cx_msy(:,1:2),2))],0.7.*[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:2),2);flipud(sum(cx_msy(:,1:3),2))],0.9.*[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:3),2);flipud(sum(cx_msy(:,1:4),2))],1.1.*[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:4),2);flipud(sum(cx_msy(:,1:5),2))],[0.9,0.9,0.9],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_msy(:,1:5),2);flipud(sum(cx_msy(:,1:6),2))],[0.9,0.7,0.9],'LineWidth',2);
legend('~olivine','~spinel','~oreutpyroxene','~clinopyroxene','~feldspar','~quartz','Interpreter','latex','FontSize',15,'box','off','location','southeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Mineral assemblage','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Mineral fraction [wt]','Interpreter','latex','FontSize',15)


% create output directory
if ~isfolder(['../out/',runID])
    mkdir(['../out/',runID]);
end

% save output to file
if save_plot
    name = ['../out/',runID,'/',runID,'_phase_dgrm'];
    print(figure(1),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_melt_model'];
    print(figure(2),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_maj_compnt'];
    print(figure(3),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_vol_compnt'];
    print(figure(4),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_density'];
    print(figure(5),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_viscosity'];
    print(figure(6),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_segr_speed'];
    print(figure(7),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_oxides'];
    print(figure(8),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_components'];
    print(figure(9),name,'-dpng','-r300','-opengl');
    name = ['../out/',runID,'/',runID,'_modal'];
    print(figure(10),name,'-dpng','-r300','-opengl');
end
