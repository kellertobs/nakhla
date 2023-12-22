% calibrate phase diagram
clear all; close all;

addpath(genpath('./'));
addpath(genpath('../'));
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
TINY = 1e-16;
FS = {'FontSize',15};
TX = {'Interpreter','latex'};

rng(15);    % for reproducibility

% set phase diagram parameters
cal_LUNA;  % load melt model calibration


%% set up data for calibration

% load experimental data from Schmidt & Kraettli (2020), Table 3
Load_SKTable3;

nphs = 7; olv=1; opx=2; cpx=3; fsp=4; spn=5; qtz=6; mlt=7; blk=8;
ncmp = 11; Si=1; Ti=2; Al=3; Cr=4; FeO=5; Mn=6; Mg=7; Ca=8; Na=9; K=10; P=11;

oxd(:,:,P ) = [];  % ignore P2O5
oxd(:,:,K ) = [];  % ignore K2O
oxd(:,:,Mn) = [];  % ignore MnO
oxd(:,:,Cr) = [];  % ignore Cr2O3

phs = phs./(sum(phs,2)+TINY);      % normalise remaining phases to 1
oxd = oxd./(sum(oxd,3)+TINY)*100;  % normalise remaining oxides to 100%

nphs = 7; olv=1; opx=2; cpx=3; fsp=4; spn=5; qtz=6; mlt=7; blk=8; xtl=9;
ncmp = 7; Si=1; Ti=2; Al=3; FeO=4; Mg=5; Ca=6; Na=7; H=8;

oxd(:,:,H) = 0;


stages = 1:size(oxd,1);

for stg = stages

    % get oxide composition of xtal assemblage
    oxd(stg,xtl,:) = squeeze(phs(stg,1:nphs-1))*squeeze(oxd(stg,1:nphs-1,:))./sum(squeeze(phs(stg,1:nphs-1)));
    oxd(stg,blk,:) = squeeze(phs(stg,1:nphs-0))*squeeze(oxd(stg,1:nphs-0,:))./sum(squeeze(phs(stg,1:nphs-0)));

end

% interpolate experimental results to finer grid to expand fitting constraints
Tmp = interp1((1:10).',Tmp,(1:0.5:10).','linear','extrap');
Prs = interp1((1:10).',Prs,(1:0.5:10).','linear','extrap');
phs = interp1((1:10).',phs,(1:0.5:10).','linear','extrap');
oxd = interp1((1:10).',oxd,(1:0.5:10).','linear','extrap');
oxd = oxd./max(1e-16,sum(oxd,3))*100;

% load data for solidus and liquidus
Psl = linspace(0,4,length(Prs))'; % all across pressure space
[Tsol, Tliq] = solidusliquidus('johnson2021', Psl);


%% collate and reduce mineral and melt oxide compositions

% detect where phases are stable
hasMLT = phs(:,mlt)>=1e-4;
hasOLV = phs(:,olv)>=1e-4 & hasMLT;
hasOPX = phs(:,opx)>=1e-4 & hasMLT;
hasCPX = phs(:,cpx)>=1e-4 & hasMLT;
hasFSP = phs(:,fsp)>=1e-4 & hasMLT;
hasSPN = phs(:,spn)>=1e-4 & hasMLT;
hasQTZ = phs(:,qtz)>=1e-4 & hasMLT;

% set oxides present in phases
oxdSYS = [Si,Ti,Al,FeO,Mg,Ca,Na,H]; noxd = length(oxdSYS);
oxdMLT = [Si,Ti,Al,FeO,Mg,Ca,Na,H];
oxdOLV = [Si,      FeO,Mg,Ca     ];
oxdOPX = [Si   ,Al,FeO,Mg,Ca     ];
oxdCPX = [Si,Ti,Al,FeO,Mg,Ca,Na  ];
oxdFSP = [Si,   Al,       Ca,Na  ];
oxdSPN = [   Ti,Al,FeO,Mg        ];
oxdQTZ = [Si                     ];

% extract oxide composition of phases
SOL = squeeze(oxd(hasMLT,xtl,oxdMLT)); SOL = SOL./sum(SOL,2)*100;
MLT = squeeze(oxd(hasMLT,mlt,oxdMLT)); MLT = MLT./sum(MLT,2)*100; nMLT = size(MLT,1);
SYS = squeeze(oxd(hasMLT,blk,oxdMLT)); SYS = SYS./sum(SYS,2)*100;
OLV = zeros(length(hasMLT(hasOLV)),length(oxdMLT)); OLV(:,oxdOLV) = squeeze(oxd(hasOLV,olv,oxdOLV)); OLV = OLV./sum(OLV,2)*100; nOLV = size(OLV,1);
OPX = zeros(length(hasMLT(hasOPX)),length(oxdMLT)); OPX(:,oxdOPX) = squeeze(oxd(hasOPX,opx,oxdOPX)); OPX = OPX./sum(OPX,2)*100; nOPX = size(OPX,1);
CPX = zeros(length(hasMLT(hasCPX)),length(oxdMLT)); CPX(:,oxdCPX) = squeeze(oxd(hasCPX,cpx,oxdCPX)); CPX = CPX./sum(CPX,2)*100; nCPX = size(CPX,1);
FSP = zeros(length(hasMLT(hasFSP)),length(oxdMLT)); FSP(:,oxdFSP) = squeeze(oxd(hasFSP,fsp,oxdFSP)); FSP = FSP./sum(FSP,2)*100; nFSP = size(FSP,1);
SPN = zeros(length(hasMLT(hasSPN)),length(oxdMLT)); SPN(:,oxdSPN) = squeeze(oxd(hasSPN,spn,oxdSPN)); SPN = SPN./sum(SPN,2)*100; nSPN = size(SPN,1);
QTZ = zeros(length(hasMLT(hasQTZ)),length(oxdMLT)); QTZ(:,oxdQTZ) = squeeze(oxd(hasQTZ,qtz,oxdQTZ)); QTZ = QTZ./sum(QTZ,2)*100; nQTZ = size(QTZ,1);
T   = Tmp(hasMLT);
P   = Prs(hasMLT)*1e9;
H2O = 0*T;

PHS    = [phs(hasMLT,mlt)*100, ...
          phs(hasMLT,olv)*100./(1-phs(hasMLT,mlt)), ...
          phs(hasMLT,opx)*100./(1-phs(hasMLT,mlt)), ...
          phs(hasMLT,cpx)*100./(1-phs(hasMLT,mlt)), ...
          phs(hasMLT,fsp)*100./(1-phs(hasMLT,mlt)), ...
          phs(hasMLT,spn)*100./(1-phs(hasMLT,mlt)), ...
          phs(hasMLT,qtz)*100./(1-phs(hasMLT,mlt))];


%% olivine system
DATA.PRJCT  = 'LUNA';
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
cal_LUNA; % load melt model calibration

figure(102); clf;
subplot(1,3,1);
scatter(OLV (:,cal.Si),OLV (:,cal.Fe),25,T(hasOLV(hasMLT))); colormap('copper'); hold on
scatter(OLVp(:,cal.Si),OLVp(:,cal.Fe),25,T(hasOLV(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(1,3,2);
scatter(OLV (:,cal.Si),OLV (:,cal.Mg),25,T(hasOLV(hasMLT))); colormap('copper'); hold on
scatter(OLVp(:,cal.Si),OLVp(:,cal.Mg),25,T(hasOLV(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(1,3,3);
scatter(OLV (:,cal.Si),OLV (:,cal.Ca),25,T(hasOLV(hasMLT))); colormap('copper'); hold on
scatter(OLVp(:,cal.Si),OLVp(:,cal.Ca),25,T(hasOLV(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

sgtitle('OLV PCA',FS{:},TX{:})
drawnow


%% orthopyroxene system
DATA.PRJCT  = 'LUNA';
DATA.VNAMES = cal.oxdStr(oxdOPX);
DATA.SNAMES = {};
DATA.X = OPX(:,oxdOPX);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
OPX_PCA = DGN;
Xp = Xp(1:nOPX,:);

OPXp           = zeros(size(OPX));
OPXp(:,oxdOPX) = max(0,Xp)./sum(max(0,Xp),2)*100;
OPX_PCA.OPXp   = OPXp;

OPX_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
OPX_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%% plot orthopyroxene system
cal_LUNA; % load melt model calibration

figure(103); clf;
subplot(2,2,1);
scatter(OPX (:,cal.Si),OPX (:,cal.Al),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Al),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(OPX (:,cal.Si),OPX (:,cal.Fe),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Fe),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,2,3);
scatter(OPX (:,cal.Si),OPX (:,cal.Mg),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Mg),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,2,4);
scatter(OPX (:,cal.Si),OPX (:,cal.Ca),25,T(hasOPX(hasMLT))); colormap('copper'); hold on
scatter(OPXp(:,cal.Si),OPXp(:,cal.Ca),25,T(hasOPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

sgtitle('OPX PCA',FS{:},TX{:})
drawnow


%% clinopyroxene system
DATA.PRJCT  = 'LUNA';
DATA.VNAMES = cal.oxdStr(oxdCPX);
DATA.SNAMES = {};
DATA.X = CPX(:,oxdCPX);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
CPX_PCA = DGN;
Xp = Xp(1:nCPX,:);

CPXp           = zeros(size(CPX));
CPXp(:,oxdCPX) = max(0,Xp)./sum(max(0,Xp),2)*100;
CPX_PCA.CPXp   = CPXp;

CPX_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
CPX_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%%
cal_LUNA; % load melt model calibration

figure(104); clf;
subplot(2,3,1);
scatter(CPX (:,cal.Si),CPX (:,cal.Ti),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Ti),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tdp,cal.Si),cal.mem_oxd(cal.tdp,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(2,3,2);
scatter(CPX (:,cal.Si),CPX (:,cal.Al),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Al),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tdp,cal.Si),cal.mem_oxd(cal.tdp,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,3);
scatter(CPX (:,cal.Si),CPX (:,cal.Fe),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Fe),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tdp,cal.Si),cal.mem_oxd(cal.tdp,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Fe),200,'kh','filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,4);
scatter(CPX (:,cal.Si),CPX (:,cal.Mg),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Mg),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tdp,cal.Si),cal.mem_oxd(cal.tdp,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,5);
scatter(CPX (:,cal.Si),CPX (:,cal.Ca),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Ca),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tdp,cal.Si),cal.mem_oxd(cal.tdp,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Ca),200,'kh','filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,6);
scatter(CPX (:,cal.Si),CPX (:,cal.Na),25,T(hasCPX(hasMLT))); colormap('copper'); hold on
scatter(CPXp(:,cal.Si),CPXp(:,cal.Na),25,T(hasCPX(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.tdp,cal.Si),cal.mem_oxd(cal.tdp,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.dps,cal.Si),cal.mem_oxd(cal.dps,cal.Na),200,'kh','filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
% colorbar;
sgtitle('CPX PCA',FS{:},TX{:})
drawnow


%% feldspar system
DATA.PRJCT  = 'LUNA';
DATA.VNAMES = cal.oxdStr(oxdFSP);
DATA.SNAMES = {};
DATA.X = FSP(:,oxdFSP);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
FSP_PCA = DGN;
Xp = Xp(1:nFSP,:);

FSPp           = zeros(size(FSP));
FSPp(:,oxdFSP) = max(0,Xp)./sum(max(0,Xp),2)*100;
FSP_PCA.FSPp   = FSPp;

FSP_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
FSP_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);


%% plot feldspar system
cal_LUNA; % load melt model calibration

figure(106); clf;
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


%% oxides system
DATA.PRJCT  = 'LUNA';
DATA.VNAMES = cal.oxdStr(oxdSPN);
DATA.SNAMES = {};
DATA.X = SPN(:,oxdSPN);

DATA.X = DATA.X./sum(DATA.X,2);

unmix
SPN_PCA = DGN;
Xp = Xp(1:nSPN,:);

% Fe             = mean(DATA.X);
% Fi             = Fe;
% Xp             = repmat(Fe,nSPN,1);

SPNp           = zeros(size(SPN));
SPNp(:,oxdSPN) = max(0,Xp)./sum(max(0,Xp),2)*100;
SPN_PCA.SPNp   = SPNp;

SPN_PCA.EMExt  = round(max(0,Fe)./sum(max(0,Fe),2)*100,2);
SPN_PCA.EMInt  = round(max(0,Fi)./sum(max(0,Fi),2)*100,2);

%%
cal_LUNA; % load melt model calibration

figure(105); clf;
subplot(1,3,1);
scatter(SPN (:,cal.Fe),SPN (:,cal.Ti),25,T(hasSPN(hasMLT))); colormap('copper'); hold on
scatter(SPNp(:,cal.Fe),SPNp(:,cal.Ti),25,T(hasSPN(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ams,cal.Fe),cal.mem_oxd(cal.ams,cal.Ti),200,'kh','filled');
scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(1,3,2);
scatter(SPN (:,cal.Fe),SPN (:,cal.Al),25,T(hasSPN(hasMLT))); colormap('copper'); hold on
scatter(SPNp(:,cal.Fe),SPNp(:,cal.Al),25,T(hasSPN(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ams,cal.Fe),cal.mem_oxd(cal.ams,cal.Al),200,'kh','filled');
scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(1,3,3);
scatter(SPN (:,cal.Fe),SPN (:,cal.Mg),25,T(hasSPN(hasMLT))); colormap('copper'); hold on
scatter(SPNp(:,cal.Fe),SPNp(:,cal.Mg),25,T(hasSPN(hasMLT)),'filled');
scatter(cal.mem_oxd(cal.ams,cal.Fe),cal.mem_oxd(cal.ams,cal.Mg),200,'kh','filled');
scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
% colorbar;
sgtitle('OXS PCA',FS{:},TX{:})
drawnow


%% add up projected mineral compositions to solid composition
wt   = zeros(size(T)) + 1e-16;
SOLp = zeros(size(SOL));

SOLp(hasOLV(hasMLT),:) = SOLp(hasOLV(hasMLT),:) + OLVp.*PHS(hasOLV(hasMLT),olv+1);
wt(hasOLV(hasMLT))     = wt(hasOLV(hasMLT)) + PHS(hasOLV(hasMLT),olv+1);

SOLp(hasOPX(hasMLT),:) = SOLp(hasOPX(hasMLT),:) + OPXp.*PHS(hasOPX(hasMLT),opx+1);
wt(hasOPX(hasMLT))     = wt(hasOPX(hasMLT)) + PHS(hasOPX(hasMLT),opx+1);

SOLp(hasCPX(hasMLT),:) = SOLp(hasCPX(hasMLT),:) + CPXp.*PHS(hasCPX(hasMLT),cpx+1);
wt(hasCPX(hasMLT))     = wt(hasCPX(hasMLT)) + PHS(hasCPX(hasMLT),cpx+1);

SOLp(hasFSP(hasMLT),:) = SOLp(hasFSP(hasMLT),:) + FSPp.*PHS(hasFSP(hasMLT),fsp+1);
wt(hasFSP(hasMLT))     = wt(hasFSP(hasMLT)) + PHS(hasFSP(hasMLT),fsp+1);

SOLp(hasSPN(hasMLT),:) = SOLp(hasSPN(hasMLT),:) + SPNp.*PHS(hasSPN(hasMLT),spn+1);
wt(hasSPN(hasMLT))     = wt(hasSPN(hasMLT)) + PHS(hasSPN(hasMLT),spn+1);

SOLp(hasQTZ(hasMLT),:) = SOLp(hasQTZ(hasMLT),:) + QTZ.*PHS(hasQTZ(hasMLT),qtz+1);
wt(hasQTZ(hasMLT))     = wt(hasQTZ(hasMLT)) + PHS(hasQTZ(hasMLT),qtz+1);

SOLp = SOLp./wt;


% collate oxide compositions into global data array
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr(oxdMLT(1:H));
DATA.SNAMES = {};
DATA.X      = [MLT(:,1:end-1);SOLp(:,1:end-1)];
DATA.X      = DATA.X./sum(DATA.X,2);

unmix;
MLT_PCA = DGN;

MLTp = MLT;
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
SYSp = SOLp.*(1-PHS(hasMLT,1)/100) + MLTp.*PHS(hasMLT,1)/100;

close all;
save('MAGEMin_processed');


%% liquid, solid, mixture compositions
cal_LUNA;  % load melt model calibration
figure(107); clf;
subplot(2,3,1);
scatter(MLT(:,Si),MLT(:,Ti),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Ti),25,T,'o','filled');
scatter(SOL(:,Si),SOL(:,Ti),25,T,'s');
scatter(SOLp(:,Si),SOLp(:,Ti),25,T,'s','filled');
scatter(SYS(:,Si),SYS(:,Ti),25,T,'d');
scatter(SYSp(:,Si),SYSp(:,Ti),25,T,'d','filled');
scatter(EMInt(:,Si),EMInt(:,Ti),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Ti),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(2,3,2);
scatter(MLT(:,Si),MLT(:,Al),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Al),25,T,'o','filled');
scatter(SOL(:,Si),SOL(:,Al),25,T,'s');
scatter(SOLp(:,Si),SOLp(:,Al),25,T,'s','filled');
scatter(SYS(:,Si),SYS(:,Al),25,T,'d');
scatter(SYSp(:,Si),SYSp(:,Al),25,T,'d','filled');
scatter(EMInt(:,Si),EMInt(:,Al),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Al),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,3);
scatter(MLT(:,Si),MLT(:,FeO),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,FeO),25,T,'o','filled');
scatter(SOL(:,Si),SOL(:,FeO),25,T,'s');
scatter(SOLp(:,Si),SOLp(:,FeO),25,T,'s','filled');
scatter(SYS(:,Si),SYS(:,FeO),25,T,'d');
scatter(SYSp(:,Si),SYSp(:,FeO),25,T,'d','filled');
scatter(EMInt(:,Si),EMInt(:,FeO),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,FeO),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,4);
scatter(MLT(:,Si),MLT(:,Mg),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Mg),25,T,'o','filled');
scatter(SOL(:,Si),SOL(:,Mg),25,T,'s');
scatter(SOLp(:,Si),SOLp(:,Mg),25,T,'s','filled');
scatter(SYS(:,Si),SYS(:,Mg),25,T,'d');
scatter(SYSp(:,Si),SYSp(:,Mg),25,T,'d','filled');
scatter(EMInt(:,Si),EMInt(:,Mg),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Mg),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,5);
scatter(MLT(:,Si),MLT(:,Ca),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Ca),25,T,'o','filled');
scatter(SOL(:,Si),SOL(:,Ca),25,T,'s');
scatter(SOLp(:,Si),SOLp(:,Ca),25,T,'s','filled');
scatter(SYS(:,Si),SYS(:,Ca),25,T,'d');
scatter(SYSp(:,Si),SYSp(:,Ca),25,T,'d','filled');
scatter(EMInt(:,Si),EMInt(:,Ca),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Ca),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,6);
scatter(MLT(:,Si),MLT(:,Na),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(MLTp(:,Si),MLTp(:,Na),25,T,'o','filled');
scatter(SOL(:,Si),SOL(:,Na),25,T,'s');
scatter(SOLp(:,Si),SOLp(:,Na),25,T,'s','filled');
scatter(SYS(:,Si),SYS(:,Na),25,T,'d');
scatter(SYSp(:,Si),SYSp(:,Na),25,T,'d','filled');
scatter(EMInt(:,Si),EMInt(:,Na),200,'kh','filled');
scatter(EMExt(:,Si),EMExt(:,Na),200,'kh');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
% colorbar;

sgtitle('MLT \& SOL PCA',FS{:},TX{:})
drawnow


%% load projected data and prepare for fitting routines
load('MAGEMin_processed');

cal_LUNA;  % load melt model calibration
                % for fay ens hyp ads dps pig ant alb ams ulv qtz wat
indmem  = logical([1   0   0   0   0   0   0   0   0   0   0   0   0
                   1   1   1   0   0   0   0   0   0   0   0   0   0
                   0   1   1   1   1   0   0   0   0   0   0   0   0
                   0   0   0   1   1   1   0   1   0   0   0   0   0
                   0   0   0   0   0   1   1   1   1   1   0   0   0
                   0   0   0   0   0   0   1   0   1   0   1   1   0
                   0   0   0   0   0   0   0   0   0   0   0   0   1]);


% convert factor analysis end-member to mineral end-member proportions
cmp_oxd = 1.0*MLT_PCA.EMInt+0.0*MLT_PCA.EMExt;
cmp_oxd = cmp_oxd./sum(cmp_oxd,2)*100;

Xp = zeros(size(cmp_oxd,1),cal.nmem);
for ip = 1:size(cmp_oxd,1)
    Xp(ip,indmem(ip,:)) = lsqnonneg(cal.mem_oxd(indmem(ip,:),oxdSYS(1:end-1)).',cmp_oxd(ip,:).');
end
cmp_mem = Xp./sum(Xp,2)*100;
cmp_mem(2:end,:) = max(2*indmem(2:end-1,:),min(75,cmp_mem(2:end,:)));
cmp_mem = cmp_mem./sum(cmp_mem,2)*100;
while max(cmp_mem(2:end,:),[],'all')>=80 || min(cmp_mem(indmem(1:end-1,:)),[],'all')<=1
    cmp_mem(2:end,:) = max(2*indmem(2:end-1,:),min(75,cmp_mem(2:end,:)));
    cmp_mem = cmp_mem./sum(cmp_mem,2)*100;
end

cmp_mem_FINT = zeros(cal.ncmp,cal.nmem);
cmp_mem_FINT(1:end-1,:) = cmp_mem;
cmp_mem_FINT(end,cal.wat) = 100;
cmp_oxd_FINT = cmp_mem_FINT*cal.mem_oxd/100;

cmp_mem_MAP = cmp_mem_FINT;

T0_MAP = [1890   1540  1480  1200  1075   1015];
A_MAP  = (cal.T0+273.15)./350;
B_MAP  = [ 8.7   6.5   5.2   2.5   2.0   1.8];
r_MAP  = [25.0  55.0  17.5  11.0   6.5   9.5];


%% Monte Carlo Parameter Fitting
cal_LUNA;  % load melt model calibration

data   = [MLTp(:);SOLp(:);PHS(:);Tsol(:);Tliq(:)];

m0     = [T0_MAP.';A_MAP.';B_MAP.';r_MAP.';cmp_mem_MAP(:).*indmem(:);];
m0_lw  = m0 - [max(10,0.05*T0_MAP.');0*max(0.1,0.01*A_MAP.');0.5*max(1,0.3*B_MAP.');max(1,0.3*r_MAP.');max(100,1*cmp_mem_MAP(:)).*indmem(:)];
m0_up  = m0 + [max(10,0.05*T0_MAP.');0*max(0.1,0.01*A_MAP.');0.5*max(1,0.3*B_MAP.');max(1,0.3*r_MAP.');max(100,1*cmp_mem_MAP(:)).*indmem(:)];
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(   cal.ncmp-1 +(1:cal.ncmp-1)       ,:) = max(2,        mbnds(   cal.ncmp-1 +(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:),min(80,mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));
mbnds(m0==100 ,:) = 100;
mbnds(m0==T0_MAP(1),:) = T0_MAP(1);

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
sigma_TsolTliq  =  0.5 * ones(size([Tsol(:);Tliq(:)]));
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
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.1,3,model,cal);

% run MCMC algorithm
Niter = 2e5;

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.002 * diff(mbnds,1,2);
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

% retrieve distributions
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% liquid, solid, mixture compositions
% cal_LUNA; % load melt model calibration
 
figure(108); clf;

subplot(2,3,1);
scatter(MLTp(:,Si),MLTp(:,Ti),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Ti),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Ti),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ti),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ti),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ti),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})

subplot(2,3,2);
scatter(MLTp(:,Si),MLTp(:,Al),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Al),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Al),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Al),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Al),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Al),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})

subplot(2,3,3);
scatter(MLTp(:,Si),MLTp(:,FeO),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,FeO),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,FeO),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,FeO),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,FeO),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,FeO),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})

subplot(2,3,4);
scatter(MLTp(:,Si),MLTp(:,Mg),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Mg),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Mg),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Mg),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Mg),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Mg),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})

subplot(2,3,5);
scatter(MLTp(:,Si),MLTp(:,Ca),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Ca),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Ca),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ca),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ca),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ca),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

subplot(2,3,6);
scatter(MLTp(:,Si),MLTp(:,Na),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Na),25,T,'s');
scatter(SYSp(:,Si),SYSp(:,Na),25,T,'d');
scatter(MLTfit(:,Si),MLTfit(:,Na),25,T,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Na),25,T,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Na),25,T,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})

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
drawnow

figure(111); clf
plot(Tm,Psl,'LineWidth',1); axis ij tight; hold on
plot(Tsolfit,Psl,'b-','LineWidth',2);
plot(Tliqfit,Psl,'r-','LineWidth',2);
plot(Tsol,Psl,'k:','LineWidth',1.5);
plot(Tliq,Psl,'k:','LineWidth',1.5);

legend([cal.cmpStr(1:end-1),'Tsol','Tliq'],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Pressure [GPa]','Interpreter','latex','FontSize',15)
drawnow


%% save and display final state of calibration
save('LUNA_calibration_olv3');

% values to enter into cal file
cmp_mem = round(cmp_mem_MAP,2)
cmp_mem = round(cmp_oxd_MAP,2)
T0      = round(T0_MAP,0)
B       = round(B_MAP,2)
r       = round(r_MAP,1)
c0      = round(cmpSYS(1,1:end-1)./sum(cmpSYS(1,1:end-1),2),2)