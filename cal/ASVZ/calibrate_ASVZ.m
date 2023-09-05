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
cal_ASVZ;  % load melt model calibration

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

%       Si Ti Al Fe Mg Ca Na  K
ioxd = [ 1  8  2  5  4  3  7  6]; % oxide indices from MAGEMin to standard
Si = 1; Ti = 2; Al = 3; FeO = 4; Mg = 5; Ca = 6; Na = 7; K = 8; H = 9;

%% load MAGEMin results

filename = 'ASVZ_fract5_out.mat';
load(filename);

% lump in free O to FeO, Cr2O3 to Al2O3, normalise to anhydrous unit sum
PHS = fieldnames(OUT.PhaseProps);
PHS = [PHS(:)',{'SYS'},{'sol'}];
for iph = 1:length(PHS)
    OUT.OxideFract.(PHS{iph}) = zeros(size(OUT.OxideFractions.(PHS{iph})));
    OUT.OxideFractions.(PHS{iph})(:,10) = 0; % no Cr %OUT.OxideFractions.(phs{iph})(:,2) + OUT.OxideFractions.(phs{iph})(:,10); OUT.OxideFractions.(phs{iph})(:,10) = 0;
    OUT.OxideFractions.(PHS{iph})(:, 5) = OUT.OxideFractions.(PHS{iph})(:,5) + OUT.OxideFractions.(PHS{iph})(:,9);  OUT.OxideFractions.(PHS{iph})(:,9)  = 0; % combine FeO + O
    OUT.OxideFract.(PHS{iph}) = OUT.OxideFractions.(PHS{iph})(:,[ioxd 11]);
    OUT.OxideFract.(PHS{iph}) = OUT.OxideFract.(PHS{iph})./sum(OUT.OxideFract.(PHS{iph})+1e-16,2);
end

% combine all feldspar instances
if isfield(OUT.PhaseProps,'pl4T2')
    OUT.OxideFract.pl4T = (OUT.OxideFract.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.OxideFract.pl4T2.*OUT.PhaseProps.pl4T2(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1)+1e-16);
    OUT.EMFractions.pl4T = (OUT.EMFractions.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.EMFractions.pl4T2.*OUT.PhaseProps.pl4T2(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1)+1e-16);
    OUT.PhaseProps.pl4T(:,6) = (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1))./(OUT.PhaseProps.pl4T(:,1)./OUT.PhaseProps.pl4T(:,6)+OUT.PhaseProps.pl4T2(:,1)./OUT.PhaseProps.pl4T2(:,6));
    OUT.PhaseProps.pl4T(:,1) = OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1);
    OUT.OxideFract = rmfield(OUT.OxideFract,'pl4T2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'pl4T2');
end

% combine all orthopyroxene instances
if isfield(OUT.PhaseProps,'opx2')
    OUT.OxideFract.opx = (OUT.OxideFract.opx.*OUT.PhaseProps.opx(:,1) + OUT.OxideFract.opx2.*OUT.PhaseProps.opx2(:,1)) ./ (OUT.PhaseProps.opx(:,1)+OUT.PhaseProps.opx2(:,1)+1e-16);
    OUT.EMFractions.opx = (OUT.EMFractions.opx.*OUT.PhaseProps.opx(:,1) + OUT.EMFractions.opx2.*OUT.PhaseProps.opx2(:,1)) ./ (OUT.PhaseProps.opx(:,1)+OUT.PhaseProps.opx2(:,1)+1e-16);
    OUT.PhaseProps.opx(:,1) = OUT.PhaseProps.opx(:,1)+OUT.PhaseProps.opx2(:,1);
    OUT.OxideFract = rmfield(OUT.OxideFract,'opx2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'opx2');
end

% combine all clinopyroxene instances
if isfield(OUT.PhaseProps,'cpx2')
    OUT.OxideFract.cpx = (OUT.OxideFract.cpx.*OUT.PhaseProps.cpx(:,1) + OUT.OxideFract.cpx2.*OUT.PhaseProps.cpx2(:,1)) ./ (OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx2(:,1)+1e-16);
    OUT.EMFractions.cpx = (OUT.EMFractions.cpx.*OUT.PhaseProps.cpx(:,1) + OUT.EMFractions.cpx2.*OUT.PhaseProps.cpx2(:,1)) ./ (OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx2(:,1)+1e-16);
    OUT.PhaseProps.cpx(:,1) = OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx2(:,1);
    OUT.OxideFract = rmfield(OUT.OxideFract,'cpx2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'cpx2');
end

% combine all spinel instances
if isfield(OUT.PhaseProps,'spn2')
    OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseProps.spn(:,1) + OUT.OxideFract.spn2.*OUT.PhaseProps.spn2(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.spn2(:,1)+1e-16);
    OUT.EMFractions.spn = (OUT.EMFractions.spn.*OUT.PhaseProps.spn(:,1) + OUT.EMFractions.spn2.*OUT.PhaseProps.spn2(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.spn2(:,1)+1e-16);
    OUT.PhaseProps.spn(:,1) = OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.spn2(:,1);
    OUT.OxideFract = rmfield(OUT.OxideFract,'spn2');
    OUT.EMFractions = rmfield(OUT.EMFractions,'spn2');
end

% lump in ilmenite with spinel
if isfield(OUT.PhaseProps,'ilm')
    OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseProps.spn(:,1) + OUT.OxideFract.ilm.*OUT.PhaseProps.ilm(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.ilm(:,1)+1e-16);
    OUT.PhaseProps.spn(:,1) = OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.ilm(:,1);
    OUT.OxideFract = rmfield(OUT.OxideFract,'ilm');
    OUT.PhaseProps = rmfield(OUT.PhaseProps,'ilm');
end


%% collate and reduce mineral and melt oxide compositions

% detect where phases are stable
hasMLT = OUT.PhaseProps.liq (:,1)>=1e-4 & OUT.PhaseFractions.sol_wt>=1e-4;% & OUT.PhaseProps.hb(:,1)==0;
hasOPX = OUT.PhaseProps.opx (:,1)>=1e-4 & hasMLT;
hasCPX = OUT.PhaseProps.cpx (:,1)>=1e-4 & hasMLT;
hasOXS = OUT.PhaseProps.spn (:,1)>=1e-4 & hasMLT;
hasFSP = OUT.PhaseProps.pl4T(:,1)>=1e-4 & hasMLT;

% set oxides present in phases
oxdSYS = [Si,Ti,Al,FeO,Mg,Ca,Na,K,H]; noxd = length(oxdSYS);
oxdMLT = [Si,Ti,Al,FeO,Mg,Ca,Na,K,H];
oxdOPX = [Si   ,Al,FeO,Mg,Ca       ];
oxdCPX = [Si,Ti,Al,FeO,Mg,Ca,Na,K  ];
oxdOXS = [   Ti,Al,FeO,Mg          ];
oxdFSP = [Si,   Al,       Ca,Na,K  ];

% extract oxide composition of phases
SOL = OUT.OxideFract.sol (hasMLT,oxdMLT).*100; SOL = SOL./sum(SOL,2)*100;
MLT = OUT.OxideFract.liq (hasMLT,oxdMLT).*100; MLT = MLT./sum(MLT,2)*100; nMLT = size(MLT,1);
SYS = (SOL.*OUT.PhaseFractions.sol_wt(hasMLT) + MLT.*OUT.PhaseFractions.liq_wt(hasMLT)) ...
    ./(     OUT.PhaseFractions.sol_wt(hasMLT) +      OUT.PhaseFractions.liq_wt(hasMLT));
OPX = zeros(length(hasMLT(hasOPX)),length(oxdMLT)); OPX(:,oxdOPX) = OUT.OxideFract.opx (hasOPX,oxdOPX).*100; OPX = OPX./sum(OPX,2)*100; nOPX = size(OPX,1);
CPX = zeros(length(hasMLT(hasCPX)),length(oxdMLT)); CPX(:,oxdCPX) = OUT.OxideFract.cpx (hasCPX,oxdCPX).*100; CPX = CPX./sum(CPX,2)*100; nCPX = size(CPX,1);
OXS = zeros(length(hasMLT(hasOXS)),length(oxdMLT)); OXS(:,oxdOXS) = OUT.OxideFract.spn (hasOXS,oxdOXS).*100; OXS = OXS./sum(OXS,2)*100; nOXS = size(OXS,1);
FSP = zeros(length(hasFSP(hasFSP)),length(oxdMLT)); FSP(:,oxdFSP) = OUT.OxideFract.pl4T(hasFSP,oxdFSP).*100; FSP = FSP./sum(FSP,2)*100; nFSP = size(FSP,1);
T   = OUT.T(hasMLT);
P   = OUT.P(hasMLT)*1e8;
H2O = OUT.OxideFract.liq(hasMLT,H)*100;

PHS    = [OUT.PhaseFractions.liq_wt(hasMLT)*100, ...
                OUT.PhaseProps.opx( hasMLT,1)*100, ...
                OUT.PhaseProps.cpx( hasMLT,1)*100, ...
                OUT.PhaseProps.spn( hasMLT,1)*100, ...
                OUT.PhaseProps.pl4T(hasMLT,1)*100, ...
                zeros(size(OUT.PhaseProps.pl4T(hasMLT,1)))] ...
                ./(OUT.PhaseFractions.sol_wt(hasMLT)+OUT.PhaseFractions.liq_wt(hasMLT));


%% orthopyroxene system
DATA.PRJCT  = 'ASVZ';
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


%% plot orthopyroxene system
cal_ASVZ; % load melt model calibration

figure(102); clf;
subplot(2,3,1);
scatter(OPX (:,cal.Si),OPX (:,cal.Ti),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Ti),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(2,3,2);
scatter(OPX (:,cal.Si),OPX (:,cal.Al),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Al),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,3);
scatter(OPX (:,cal.Si),OPX (:,cal.Fe),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Fe),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,4);
scatter(OPX (:,cal.Si),OPX (:,cal.Mg),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Mg),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,5);
scatter(OPX (:,cal.Si),OPX (:,cal.Ca),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Ca),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

sgtitle('OPX PCA',FS{:},TX{:})
drawnow


%% clinopyroxene system
DATA.PRJCT  = 'ASVZ';
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
cal_ASVZ; % load melt model calibration

figure(103); clf;
subplot(2,4,1);
scatter(CPX (:,cal.Si),CPX (:,cal.Ti),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Ti),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(2,4,2);
scatter(CPX (:,cal.Si),CPX (:,cal.Al),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Al),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,4,3);
scatter(CPX (:,cal.Si),CPX (:,cal.Fe),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Fe),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,4,4);
scatter(CPX (:,cal.Si),CPX (:,cal.Mg),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Mg),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,4,5);
scatter(CPX (:,cal.Si),CPX (:,cal.Ca),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Ca),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,4,6);
scatter(CPX (:,cal.Si),CPX (:,cal.Na),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Na),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
subplot(2,4,7);
scatter(CPX (:,cal.Si),CPX (:,cal.K ),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.K ),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.K ),200,'kh','filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.K ),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.K ),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.K ),FS{:},TX{:})
% colorbar;
sgtitle('CPX PCA',FS{:},TX{:})
drawnow


%% oxides system
DATA.PRJCT  = 'ASVZ';
DATA.VNAMES = cal.oxdStr(oxdOXS);
DATA.SNAMES = {};
DATA.X = OXS(:,oxdOXS);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
OXS_PCA = DGN;

OXSp           = zeros(size(OXS));
OXSp(:,oxdOXS) = max(0,Xp)./sum(max(0,Xp),2)*100;
OXS_PCA.OXSp   = OXSp;

OXS_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
OXS_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%%
cal_ASVZ; % load melt model calibration

figure(104); clf;
subplot(1,3,1);
scatter(OXS (:,cal.Fe),OXS (:,cal.Ti),25,T(hasOXS(hasMLT))); colormap('copper'); hold on
scatter(OXSp(:,cal.Fe),OXSp(:,cal.Ti),25,T(hasOXS(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.msp,cal.Fe),cal.mem_oxd(cal.msp,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.tsp,cal.Fe),cal.mem_oxd(cal.tsp,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(1,3,2);
scatter(OXS (:,cal.Fe),OXS (:,cal.Al),25,T(hasOXS(hasMLT))); colormap('copper'); hold on
scatter(OXSp(:,cal.Fe),OXSp(:,cal.Al),25,T(hasOXS(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.msp,cal.Fe),cal.mem_oxd(cal.msp,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.tsp,cal.Fe),cal.mem_oxd(cal.tsp,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(1,3,3);
scatter(OXS (:,cal.Fe),OXS (:,cal.Mg),25,T(hasOXS(hasMLT))); colormap('copper'); hold on
scatter(OXSp(:,cal.Fe),OXSp(:,cal.Mg),25,T(hasOXS(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.msp,cal.Fe),cal.mem_oxd(cal.msp,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.tsp,cal.Fe),cal.mem_oxd(cal.tsp,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
% colorbar;
sgtitle('OXS PCA',FS{:},TX{:})
drawnow


%% feldspar system
DATA.PRJCT  = 'ASVZ';
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
cal_ASVZ; % load melt model calibration

figure(105); clf;
subplot(2,2,1);
scatter(FSP (:,cal.Si),FSP (:,cal.Al),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.Al),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(FSP (:,cal.Si),FSP (:,cal.Ca),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.Ca),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,2,3);
scatter(FSP (:,cal.Si),FSP (:,cal.Na),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.Na),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
subplot(2,2,4);
scatter(FSP (:,cal.Si),FSP (:,cal.K ),25,T(hasFSP(hasMLT))); colormap('copper'); hold on
scatter(FSPp(:,cal.Si),FSPp(:,cal.K ),25,T(hasFSP(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.K ),200,'kh','filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.K ),200,'kh','filled');
scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.K ),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.K ),FS{:},TX{:})
% colorbar;
sgtitle('FSP PCA',FS{:},TX{:})
drawnow


%% add up projected mineral compositions to solid composition
wt   = zeros(size(T)) + 1e-16;
SOLp = zeros(size(SOL));

SOLp(hasOPX(hasMLT),:) = SOLp(hasOPX(hasMLT),:) + OPXp.*OUT.PhaseProps.opx(hasOPX,1);
wt(hasOPX(hasMLT))     = wt(hasOPX(hasMLT)) + OUT.PhaseProps.opx(hasOPX,1);

SOLp(hasCPX(hasMLT),:) = SOLp(hasCPX(hasMLT),:) + CPXp.*OUT.PhaseProps.cpx(hasCPX,1);
wt(hasCPX(hasMLT))     = wt(hasCPX(hasMLT)) + OUT.PhaseProps.cpx(hasCPX,1);

SOLp(hasOXS(hasMLT),:) = SOLp(hasOXS(hasMLT),:) + OXSp.*OUT.PhaseProps.spn(hasOXS,1);
wt(hasOXS(hasMLT))     = wt(hasOXS(hasMLT)) + OUT.PhaseProps.spn(hasOXS,1);

SOLp(hasFSP(hasMLT),:) = SOLp(hasFSP(hasMLT),:) + FSPp.*OUT.PhaseProps.pl4T(hasFSP,1);
wt(hasFSP(hasMLT))     = wt(hasFSP(hasMLT)) + OUT.PhaseProps.pl4T(hasFSP,1);

SOLp = SOLp./wt;


%% collate oxide compositions into global data array
DATA.PRJCT  = 'ASVZ';
DATA.VNAMES = cal.oxdStr(oxdMLT(1:K));
DATA.SNAMES = {};
DATA.X      = [SOLp(:,1:K);MLT(:,1:K)];
DATA.X      = DATA.X./sum(DATA.X,2);

unmix;
MLT_PCA = DGN;

MLTp = MLT;
MLTp(:,1:cal.K) = max(0,Xp(1+nMLT:end,:))./sum(max(0,Xp(1+nMLT:end,:)),2).*(100-MLT(:,H));
MLT_PCA.MLTp    = MLTp;

EMExt = max(0,Fe)./sum(max(0,Fe),2)*100;
EMInt = max(0,Fi)./sum(max(0,Fi),2)*100;

EMExt(EMExt(:,Al)==max(EMExt(:,Al)),:) = [];
[~,iSi] = sort(EMExt(:,Si),'ascend');
EMExt = EMExt(iSi,:);
EMExt = [cal.mem_oxd(cal.ant,1:K);EMExt];

EMInt(EMInt(:,cal.Al)==max(EMInt(:,cal.Al)),:) = [];
[~,iSi] = sort(EMInt(:,Si),'ascend');
EMInt = EMInt(iSi,:);
EMInt = [cal.mem_oxd(cal.ant,1:K);EMInt];

MLT_PCA.EMInt = EMInt;
MLT_PCA.EMExt = EMExt;


%% reconstruct bulk composition based on projected solid and liquid compositions
SYSp = (SOLp.*OUT.PhaseFractions.sol_wt(hasMLT) + MLTp.*OUT.PhaseFractions.liq_wt(hasMLT)) ...
     ./(      OUT.PhaseFractions.sol_wt(hasMLT) +       OUT.PhaseFractions.liq_wt(hasMLT));

close all;
save('MAGEMin_processed');


%% liquid, solid, mixture compositions
cal_ASVZ;  % load melt model calibration
figure(106); clf;
subplot(2,4,1);
scatter(MLT(:,Si),MLT(:,Ti),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Ti),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Ti),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Ti),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Ti),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Ti),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(2,4,2);
scatter(MLT(:,Si),MLT(:,Al),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Al),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Al),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Al),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Al),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Al),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,4,3);
scatter(MLT(:,Si),MLT(:,FeO),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,FeO),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,FeO),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,FeO),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,FeO),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,FeO),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,FeO),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,4,4);
scatter(MLT(:,Si),MLT(:,Mg),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Mg),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Mg),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Mg),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Mg),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Mg),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,4,5);
scatter(MLT(:,Si),MLT(:,Ca),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Ca),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Ca),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Ca),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Ca),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Ca),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,4,6);
scatter(MLT(:,Si),MLT(:,Na),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Na),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,Na),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,Na),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,Na),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,Na),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
subplot(2,4,7);
scatter(MLT(:,Si),MLT(:,K),25,OUT.T(hasMLT),'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,K),25,OUT.T(hasMLT),'o','filled');
scatter(SOL(:,Si),SOL(:,K),25,OUT.T(hasMLT),'s');
scatter(SOLp(:,Si),SOLp(:,K),25,OUT.T(hasMLT),'s','filled');
scatter(SYS(:,Si),SYS(:,K),25,OUT.T(hasMLT),'d');
scatter(SYSp(:,Si),SYSp(:,K),25,OUT.T(hasMLT),'d','filled');
scatter(EMInt(:,Si),EMInt(:,K),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.K),FS{:},TX{:})
subplot(2,4,8);
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
cal_ASVZ;  % load melt model calibration

indmem  = logical([0   0   0   0   0   0   0   1   0   0   0   0
                   0   0   0   0   0   1   1   1   1   0   0   0
                   1   0   1   0   0   1   1   1   1   0   0   0
                   1   1   0   1   0   1   1   1   1   1   0   0
                   0   1   0   0   1   0   1   0   1   1   1   0
                   0   0   0   0   0   0   0   0   0   0   0   1]);


%% convert factor analysis end-member to mineral end-member proportions

cmp_oxd = MLT_PCA.EMInt;
cmp_oxd = cmp_oxd./sum(cmp_oxd,2)*100;
cmp_oxd(cmp_oxd(:,cal.Al)==max(cmp_oxd(:,cal.Al)),:) = [];
[~,iSi] = sort(cmp_oxd(:,Si),'ascend');
cmp_oxd_FINT = cmp_oxd(iSi,:);

Xp = zeros(size(cmp_oxd,1),cal.nmem);
for ip = 1:size(cmp_oxd,1)
    Xp(ip,indmem(ip+1,:)) = lsqnonneg(cal.mem_oxd(indmem(ip+1,:),oxdSYS(1:K)).',cmp_oxd_FINT(ip,:).');
end
cmp_mem = Xp./sum(Xp,2)*100;

cmp_mem_FINT = zeros(cal.ncmp,cal.nmem);
cmp_mem_FINT(2:end-1,:) = cmp_mem;
cmp_mem_FINT(1,cal.ant) = 100;
cmp_mem_FINT(end,cal.wat) = 100;
cmp_oxd_FINT = cmp_mem_FINT*cal.mem_oxd/100;

cmp_mem_MAP = cmp_mem_FINT;


%%
cal_ASVZ;  % load melt model calibration

data  = [MLTp(:);SOLp(:);PHS(:)];

m0     = cmp_mem_MAP.*indmem;
m0_lw  = max(0,floor(m0 - max(1,0.1*m0)).*indmem);
m0_up  = max(0, ceil(m0 + max(1,0.1*m0)).*indmem);
m0     = m0(:);

mNames = cell(cal.ncmp*cal.nmem,1);
k = 1;
for j=1:cal.nmem
    for i=1:cal.ncmp
        mNames{k} = [cal.cmpStr{i},':',cal.memStr{j}];
        k = k+1;
    end
end

mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(m0==100,:) = 100;

sigma  = max(0.02,0.02.*data);

% function to calculate forward model
% m --> dhat
dhatFunc  = @(model) OxdFromCmpMem(model,data,SYSp,cal);

% function to apply further constraints to a proposed set of model param values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% function to sample from the prior probability distributionn
% [] --> set of model parameter values [Niter x Nm]
PrSmpFunc = @(Niter) ProbFuncs('PriorSampFunc', Niter, mbnds, 'uniform');

% function to calculate likelihood of dhat
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,5,model,cal);

% function to calculate likelihood from model parameter values
% model --> dhat --> likelihood
LkMdFunc  = @(model) ProbFuncs('LikeFuncModel', dhatFunc, model, data, sigma);

% run MCMC algorithm
Niter = 1e5;

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.05 * diff(mbnds,1,2);
anneal.levels   = 3;
anneal.burnin   = Niter/20;
anneal.refine   = Niter/10;

tic;
[models_cmp,prob_cmp,accept_cmp,bestfit_cmp] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,m0,mbnds,anneal,Niter);
RunTime(1) = toc;

% plot mcmc outputs
plotmcmc(models_cmp, prob_cmp, [], mbnds, anneal, mNames);

cmp_mem_MAP = reshape(bestfit_cmp,cal.ncmp,cal.nmem);
cmp_oxd_MAP = cmp_mem_MAP*cal.mem_oxd/100;
dhat        = dhatFunc(bestfit_cmp);
oxdfit      = dhat(1:2*length(T)*cal.noxd);
phsfit      = reshape(dhat(2*length(T)*cal.noxd+1:end),[],cal.nmsy+1);
MLTfit      = reshape(oxdfit(1:length(T)*cal.noxd,:),[],cal.noxd);
SOLfit      = reshape(oxdfit(length(T)*cal.noxd+1:end,:),[],cal.noxd);
SYSfit      = MLTfit.*phsfit(:,1)/100 + SOLfit.*(1-phsfit(:,1)/100);

OPTIONS.TolX = 1e-11;
Xp = zeros(length(SYSfit),cal.ncmp);
for ip = 1:length(SYSfit)
    Xp(ip,:) = lsqnonneg(cmp_oxd_MAP.',SYSfit(ip,:).',OPTIONS);
end
c = Xp./sum(Xp,2);

% retrieve distributions
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% liquid, solid, mixture compositions
% cal_ASVZ; % load melt model calibration
 
figure(107); clf;
subplot(2,4,1);
scatter(MLT(:,Si),MLT(:,Ti),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,Ti),25,T,'s');
scatter(SYS(:,Si),SYS(:,Ti),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ti),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ti),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ti),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})

subplot(2,4,2);
scatter(MLT(:,Si),MLT(:,Al),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,Al),25,T,'s');
scatter(SYS(:,Si),SYS(:,Al),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Al),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Al),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Al),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})

subplot(2,4,3);
scatter(MLT(:,Si),MLT(:,FeO),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,FeO),25,T,'s');
scatter(SYS(:,Si),SYS(:,FeO),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,FeO),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,FeO),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,FeO),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})

subplot(2,4,4);
scatter(MLT(:,Si),MLT(:,Mg),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,Mg),25,T,'s');
scatter(SYS(:,Si),SYS(:,Mg),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Mg),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Mg),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Mg),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})

subplot(2,4,5);
scatter(MLT(:,Si),MLT(:,Ca),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,Ca),25,T,'s');
scatter(SYS(:,Si),SYS(:,Ca),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ca),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ca),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ca),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

subplot(2,4,6);
scatter(MLT(:,Si),MLT(:,Na),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,Na),25,T,'s');
scatter(SYS(:,Si),SYS(:,Na),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Na),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Na),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Na),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})

subplot(2,4,7);
scatter(MLT(:,Si),MLT(:,K),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,K),25,T,'s');
scatter(SYS(:,Si),SYS(:,K),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,K),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,K),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,K),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.K),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.K),FS{:},TX{:})

subplot(2,4,8);
scatter(MLT(:,Si),MLT(:,H),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOL(:,Si),SOL(:,H),25,T,'s');
scatter(SYS(:,Si),SYS(:,H),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,H),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,H),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,H),25,T,'d','filled');
% scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.H),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.H),FS{:},TX{:})

sgtitle('MCMC component fit',FS{:},TX{:})

figure(108); clf; cmap = colororder;
plot(T,PHS(:,1),'-',T,phsfit(:,1),'--','Color',cmap(2,:),'LineWidth',1.5); axis tight; hold on % melt
plot(T,PHS(:,2),'-',T,phsfit(:,2),'--','Color',cmap(6,:),'LineWidth',1.5); % opx
plot(T,PHS(:,3),'-',T,phsfit(:,3),'--','Color',cmap(3,:),'LineWidth',1.5); % cpx
plot(T,PHS(:,4),'-',T,phsfit(:,4),'--','Color',cmap(4,:),'LineWidth',1.5); % spn
plot(T,PHS(:,5),'-',T,phsfit(:,5),'--','Color',cmap(1,:),'LineWidth',1.5); % fsp
plot(T,PHS(:,6),'-',T,phsfit(:,6),'--','Color',cmap(7,:),'LineWidth',1.5); % qtz
xlabel('Temperature [$^\circ$C]',FS{:},TX{:})
ylabel('Phase proportions [wt\%]',FS{:},TX{:})
sgtitle('Phase stability',FS{:},TX{:})
drawnow

%% best fit melting temperatures
clear cal var x m f cx cm
cal_ASVZ;  % load melt model calibration

OPTIONS.TolX = 1e-16;

% MLTfitH = MLTfit;
% SYSfitH = SYSfit;
% MLTfitH(:,H) = H2O                 ; MLTfitH(:,1:end-1) = MLTfitH(:,1:end-1)./sum(MLTfitH(:,1:end-1),2).*(100-MLTfitH(:,end));
% SYSfitH(:,H) = H2O.*phsfit(:,1)/100; SYSfitH(:,1:end-1) = SYSfitH(:,1:end-1)./sum(SYSfitH(:,1:end-1),2).*(100-SYSfitH(:,end));

Xp = zeros(length(SYSfit),cal.ncmp);
for ip = 1:length(SYSfit)
    Xp(ip,:) = lsqnonneg(cal.cmp_oxd.',SYSfit(ip,:).',OPTIONS);
end
c = Xp./sum(Xp,2);

% equilibrium phase fractions and compositions

data   = [MLTfit(:);SOLfit(:)];

m0    = [1553.0  1200.0  1100.0  1000.0  700.0, 32.0    4.0    8.0   16.0   22.0].';
m0_lw  = max(4,floor(m0 - max(1,0.05*m0)));
m0_up  = max(4, ceil(m0 + max(1,0.05*m0)));
m0_lw(1) = m0(1);
m0_up(1) = m0(1);

mNames = cell((cal.ncmp-1)*2,1);
k = 1;
for j=1:2
    for i=1:cal.ncmp-1
        if j==1
            mNames{k} = [cal.cmpStr{i},':T0'];
        else
            mNames{k} = [cal.cmpStr{i},':r'];
        end
        k = k+1;
    end
end

mbnds  = [m0_lw,m0_up]; % model parameter bounds

sigma  = max(0.05,0.01.*data);
% sigma(length(data)/2+1:end) = sigma(length(data)/2+1:end)*2;

% function to calculate forward model
% m --> dhat
dhatFunc  = @(model) OxdFromMeltTemp(model,T,P,MLTfit(:,H),c,cal);

% function to apply further constraints to a proposed set of model param values
% m --> prior prob
ConstrFunc = @(model) ConstrFuncs('NoConstr', model);

% function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% function to sample from the prior probability distributionn
% [] --> set of model parameter values [Niter x Nm]
PrSmpFunc = @(Niter) ProbFuncs('PriorSampFunc', Niter, mbnds, 'uniform');

% function to calculate likelihood of dhat
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFunc',dhat,data,sigma,model);

% function to calculate likelihood from model parameter values
% model --> dhat --> likelihood
LkMdFunc  = @(model) ProbFuncs('LikeFuncModel', dhatFunc, model, data, sigma);

% run MCMC algorithm
Niter = 1e4;
Nbins = min(500,Niter/20);

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.025 * diff(mbnds,1,2);
anneal.levels   = 3;
anneal.burnin   = Niter/20;
anneal.refine   = Niter/10;

tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,m0,mbnds,anneal,Niter);
RunTime(1) = toc;

% plot mcmc outputs
xMAP = plotmcmc(models, prob, [], mbnds, anneal, mNames);

% xMAP = [cal.T0,cal.r];
T0_MAP = xMAP(1:cal.ncmp-1);
r_MAP  = xMAP(cal.ncmp:end);
dhat   = dhatFunc(xMAP.');
cm_oxd_MAP = reshape(dhat(1:length(dhat)/2),[],cal.noxd);%*cmp_oxd_MAP;
cx_oxd_MAP = reshape(dhat(length(dhat)/2+1:end),[],cal.noxd);%*cmp_oxd_MAP;
c_oxd_MAP  = c*cmp_oxd_MAP;

% retrieve distributions
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% plot phase diagram
figure(6); clf;
subplot(3,3,1)
plot(MLTp(:,cal.Si)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Si)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Si)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Si)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Si)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Si)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Si)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Si)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Si)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Si},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,2)
plot(MLTp(:,cal.Ti)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Ti)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Ti)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Ti)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Ti)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Ti)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Ti)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Ti)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Ti)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Ti},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,3)
plot(MLTp(:,cal.Al)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Al)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Al)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Al)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Al)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Al)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Al)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Al)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Al)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Al},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,4)
plot(MLTp(:,cal.Fe)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Fe)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Fe)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Fe)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Fe)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Fe)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Fe)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Fe)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Fe)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Fe},' [wt]'],'Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

subplot(3,3,5)
plot(MLTp(:,cal.Mg)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Mg)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Mg)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Mg)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Mg)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Mg)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Mg)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Mg)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Mg)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Mg},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,6)
plot(MLTp(:,cal.Ca)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Ca)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Ca)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Ca)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Ca)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Ca)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Ca)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Ca)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Ca)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Ca},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,7)
plot(MLTp(:,cal.Na)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.Na)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.Na)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.Na)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.Na)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.Na)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Na)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.Na)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.Na)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.Na},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,8)
plot(MLTp(:,cal.K)./sum(MLTp(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.K)./sum(SOLp(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.K)./sum(SYSp(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.K)./sum(MLTfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.K)./sum(SOLfit(:,1:end-1),2),T,'bs');
plot(SYSfit(:,cal.K)./sum(SYSfit(:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.K)./sum(c_oxd_MAP (:,1:end-1),2),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.K)./sum(cx_oxd_MAP(:,1:end-1),2),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.K)./sum(cm_oxd_MAP(:,1:end-1),2),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.K},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,9)
plot(MLTp(:,cal.H),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,cal.H),T,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,cal.H),T,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,cal.H),T,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,cal.H),T,'bs');
plot(SYSfit(:,cal.H),T,'kd');

plot(c_oxd_MAP (:,cal.H),T,'kd','LineWidth',2);
plot(cx_oxd_MAP(:,cal.H),T,'bs','LineWidth',2);
plot(cm_oxd_MAP(:,cal.H),T,'ro','LineWidth',2);
xlabel([cal.oxdStr{cal.H},' [wt]'],'Interpreter','latex','FontSize',15)

PlotPhaseDiagrams;


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
