% calibrate phase diagram
clear all; close all;

addpath('../../cal')
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
TINY = 1e-16;
FS = {'FontSize',15};
TX = {'Interpreter','latex'};

% run('../../usr/par_default');  % load default parameters

% calibration run options
runID     = 'cal_MORB';           % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set phase diagram parameters
cal_MORB;  % load melt model calibration

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

if ~holdfig; close all; end

%% load MAGEMin results
nc   = [1]; % number of compositions modelled
frct = [4];
MAG  = [];
ioxd = [1 8 2 5 4 3 7 6]; % oxide indices from MAGEMin to standard
Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; H = 8;
for ic = nc
    filename = ['MORB_fract',int2str(frct(ic)),'_H05_noK_out.mat'];
    load(filename);

    % lump in free O to FeO, Cr2O3 to Al2O3, normalise to anhydrous unit sum
    phs = fieldnames(OUT.PhaseProps);
    phs = [phs(:)',{'SYS'},{'sol'}];
    for iph = 1:length(phs)
        OUT.OxideFract.(phs{iph}) = zeros(size(OUT.OxideFractions.(phs{iph})));
        OUT.OxideFractions.(phs{iph})(:,7) = OUT.OxideFractions.(phs{iph})(:,7) + OUT.OxideFractions.(phs{iph})(:,6);  OUT.OxideFractions.(phs{iph})(:,6)  = 0;
        OUT.OxideFractions.(phs{iph})(:,5) = OUT.OxideFractions.(phs{iph})(:,5) + OUT.OxideFractions.(phs{iph})(:,9);  OUT.OxideFractions.(phs{iph})(:,9)  = 0;
        OUT.OxideFractions.(phs{iph})(:,2) = OUT.OxideFractions.(phs{iph})(:,2) + OUT.OxideFractions.(phs{iph})(:,10); OUT.OxideFractions.(phs{iph})(:,10) = 0;
        OUT.OxideFract.(phs{iph}) = OUT.OxideFractions.(phs{iph})(:,[ioxd 11]);
        OUT.OxideFract.(phs{iph}) = OUT.OxideFract.(phs{iph})./sum(OUT.OxideFract.(phs{iph})+1e-16,2);
    end

    % combine all feldspar instances
    if isfield(OUT.PhaseProps,'pl4T2')
        OUT.OxideFract.pl4T = (OUT.OxideFract.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.OxideFract.pl4T2.*OUT.PhaseProps.pl4T2(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1)+1e-16);
        OUT.EMFractions.pl4T = (OUT.EMFractions.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.EMFractions.pl4T2.*OUT.PhaseProps.pl4T2(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1)+1e-16);
        OUT.PhaseProps.pl4T = OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T2(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'pl4T2');
        OUT.EMFractions = rmfield(OUT.EMFractions,'pl4T2');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'pl4T2');
    end
    if isfield(OUT.PhaseProps,'pl4T3')
        OUT.OxideFract.pl4T = (OUT.OxideFract.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.OxideFract.pl4T3.*OUT.PhaseProps.pl4T3(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T3(:,1)+1e-16);
        OUT.EMFractions.pl4T = (OUT.EMFractions.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.EMFractions.pl4T3.*OUT.PhaseProps.pl4T3(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T3(:,1)+1e-16);
        OUT.PhaseProps.pl4T = OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.pl4T3(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'pl4T3');
        OUT.EMFractions = rmfield(OUT.EMFractions,'pl4T3');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'pl4T3');
    end

    % combine all orthopyroxene instances
    if isfield(OUT.PhaseProps,'opx2')
        OUT.OxideFract.opx = (OUT.OxideFract.opx.*OUT.PhaseProps.opx(:,1) + OUT.OxideFract.opx2.*OUT.PhaseProps.opx2(:,1)) ./ (OUT.PhaseProps.opx(:,1)+OUT.PhaseProps.opx2(:,1)+1e-16);
        OUT.EMFractions.opx = (OUT.EMFractions.opx.*OUT.PhaseProps.opx(:,1) + OUT.EMFractions.opx2.*OUT.PhaseProps.opx2(:,1)) ./ (OUT.PhaseProps.opx(:,1)+OUT.PhaseProps.opx2(:,1)+1e-16);
        OUT.PhaseProps.opx = OUT.PhaseProps.opx(:,1)+OUT.PhaseProps.opx2(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'opx2');
        OUT.EMFractions = rmfield(OUT.EMFractions,'opx2');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'opx2');
    end

    % combine all clinopyroxene instances
    if isfield(OUT.PhaseProps,'cpx2')
        OUT.OxideFract.cpx = (OUT.OxideFract.cpx.*OUT.PhaseProps.cpx(:,1) + OUT.OxideFract.cpx2.*OUT.PhaseProps.cpx2(:,1)) ./ (OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx2(:,1)+1e-16);
        OUT.EMFractions.cpx = (OUT.EMFractions.cpx.*OUT.PhaseProps.cpx(:,1) + OUT.EMFractions.cpx2.*OUT.PhaseProps.cpx2(:,1)) ./ (OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx2(:,1)+1e-16);
        OUT.PhaseProps.cpx = OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx2(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'cpx2');
        OUT.EMFractions = rmfield(OUT.EMFractions,'cpx2');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'cpx2');
    end
    if isfield(OUT.PhaseProps,'cpx3')
        OUT.OxideFract.cpx = (OUT.OxideFract.cpx.*OUT.PhaseProps.cpx(:,1) + OUT.OxideFract.cpx3.*OUT.PhaseProps.cpx3(:,1)) ./ (OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx3(:,1)+1e-16);
        OUT.EMFractions.cpx = (OUT.EMFractions.cpx.*OUT.PhaseProps.cpx(:,1) + OUT.EMFractions.cpx3.*OUT.PhaseProps.cpx3(:,1)) ./ (OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx3(:,1)+1e-16);
        OUT.PhaseProps.cpx = OUT.PhaseProps.cpx(:,1)+OUT.PhaseProps.cpx3(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'cpx3');
        OUT.EMFractions = rmfield(OUT.EMFractions,'cpx3');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'cpx3');
    end

    % combine all spinel instances
    if isfield(OUT.PhaseProps,'spn2')
        OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseProps.spn(:,1) + OUT.OxideFract.spn2.*OUT.PhaseProps.spn2(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.spn2(:,1)+1e-16);
        OUT.EMFractions.spn = (OUT.EMFractions.spn.*OUT.PhaseProps.spn(:,1) + OUT.EMFractions.spn2.*OUT.PhaseProps.spn2(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.spn2(:,1)+1e-16);
        OUT.PhaseProps.spn = OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.spn2(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'spn2');
        OUT.EMFractions = rmfield(OUT.EMFractions,'spn2');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'spn2');
    end

    % lump in ilmenite with spinel
    if isfield(OUT.PhaseProps,'ilm')
        OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseProps.spn(:,1) + OUT.OxideFract.ilm.*OUT.PhaseProps.ilm(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.ilm(:,1)+1e-16);
        OUT.PhaseProps.spn = OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.ilm(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'ilm');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'ilm');
    end

    MAG(ic).OUT = OUT;
end

% calibrate mineral end-members

%% olivine system
cal_MORB;  Fe=4;  % load melt model calibration
DATA.PRJCT  = 'BSE';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Fe,cal.Mg]);
DATA.SNAMES = {};
DATA.X      = [];
DATA.T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'ol')
        hasolv = MAG(ic).OUT.PhaseProps.ol(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
        DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.ol(hasolv,[Si,Fe,Mg]).*100];
        DATA.T = [DATA.T;MAG(ic).OUT.T(hasolv)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix

%% plot olivine system
figure(1); clf;
subplot(1,2,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Fe),100,[0.1 0.6 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),100,[0.1 0.6 0.1]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Fe),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Fe),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.olv,cal.Fe),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(1,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Mg),100,[0.1 0.6 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),100,[0.1 0.6 0.1]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Mg),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Mg),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.olv,cal.Mg),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
colorbar;
sgtitle('olivine system',FS{:},TX{:})

Fe=4;
for ic = nc
    MAG(ic).OUT.OxideFract.olp = zeros(size(MAG(ic).OUT.OxideFract.ol));
    MAG(ic).OUT.OxideFract.olp(hasolv,[Si,Fe,Mg]) = Xp;
end

%% orthopyroxene system
cal_MORB; % load melt model calibration
DATA.PRJCT  = 'BSE';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Al,cal.Fe,cal.Mg,cal.Ca]);
DATA.SNAMES = {};
DATA.X      = [];
DATA.T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'opx')
    hasopx = MAG(ic).OUT.PhaseProps.opx(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.opx(hasopx,[Si,Al,Fe,Mg,Ca]).*100];
    DATA.T = [DATA.T;MAG(ic).OUT.T(hasopx)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix

%% plot orthopyroxene system
figure(3); clf;
subplot(2,2,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),100,[0.4 0.2 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Al),100,[0.4 0.2 0.1]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Al),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Al),100,cal.T0(cal.bas),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),100,[0.4 0.2 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Fe),100,[0.4 0.2 0.1]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Fe),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Fe),100,cal.T0(cal.bas),'s','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,2,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),100,[0.4 0.2 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Mg),100,[0.4 0.2 0.1]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Mg),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Mg),100,cal.T0(cal.bas),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,2,4);
scatter(X(:,1)*100,X(:,5)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,5)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),100,[0.4 0.2 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Ca),100,[0.4 0.2 0.1]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Ca),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Ca),100,cal.T0(cal.bas),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
colorbar;
sgtitle('orthopyroxene system',FS{:},TX{:})

Fe=4;
for ic = nc
    MAG(ic).OUT.OxideFract.opxp = zeros(size(MAG(ic).OUT.OxideFract.opx));
    MAG(ic).OUT.OxideFract.opxp(hasopx,[Si,Al,Fe,Mg,Ca]) = Xp;
end

%% clinopyroxene system
cal_MORB; % load melt model calibration
DATA.PRJCT  = 'BSE';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Ti,cal.Al,cal.Fe,cal.Mg,cal.Ca,cal.Na]);
DATA.SNAMES = {};
DATA.X      = [];
DATA.T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hascpx = MAG(ic).OUT.PhaseProps.cpx(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.cpx(hascpx,[Si,Ti,Al,Fe,Mg,Ca,Na]).*100];
    DATA.T = [DATA.T;MAG(ic).OUT.T(hascpx)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix

%%
figure(4); clf;
subplot(2,3,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Ti),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ti),100,[0.6 0.0 0.3]*1.0,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Ti),100,[0.6 0.0 0.3]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Ti),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Ti),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Ti),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(2,3,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Al),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Al),100,[0.6 0.0 0.3]*1.0,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Al),100, [0.6 0.0 0.3]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Al),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Al),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Al),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Fe),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Fe),100,[0.6 0.0 0.3]*1.0,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Fe),100,[0.6 0.0 0.3]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Fe),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Fe),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Fe),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,4);
scatter(X(:,1)*100,X(:,5)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,5)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Mg),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Mg),100,[0.6 0.0 0.3]*1.0,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Mg),100,[0.6 0.0 0.3]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Mg),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Mg),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Mg),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,5);
scatter(X(:,1)*100,X(:,6)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,6)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Ca),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ca),100,[0.6 0.0 0.3]*1.0,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Ca),100,[0.6 0.0 0.3]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Ca),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Ca),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Ca),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,6);
scatter(X(:,1)*100,X(:,7)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,7)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Na),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Na),100,[0.6 0.0 0.3]*1.0,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Na),100,[0.6 0.0 0.3]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Na),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Na),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Na),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
colorbar;
sgtitle('clinopyroxene system',FS{:},TX{:})

Fe=4;
for ic = nc
    MAG(ic).OUT.OxideFract.cpxp = zeros(size(MAG(ic).OUT.OxideFract.cpx));
    MAG(ic).OUT.OxideFract.cpxp(hascpx,[Si,Ti,Al,Fe,Mg,Ca,Na]) = Xp;
end

%% oxides system
cal_MORB; Fe = 4; % load melt model calibration
DATA.PRJCT  = 'andesSVZ';
DATA.VNAMES = cal.oxdStr([cal.Fe,cal.Ti,cal.Al,cal.Mg]);
DATA.SNAMES = {};
DATA.X      = [];
DATA.T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'spn')
        hasoxs = MAG(ic).OUT.PhaseProps.spn(:,1)>0.0001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
        DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.spn(hasoxs,[Fe,Ti,Al,Mg]).*100];
        DATA.T = [DATA.T;MAG(ic).OUT.T(hasoxs)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix

%%
figure(2); clf;
subplot(1,3,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Ti),100,[0.4 0.4 0.4]*1.3,'filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Ti),100,[0.4 0.4 0.4]*1.0,'filled');
scatter(cal.mem_oxd(cal.tim,cal.Fe),cal.mem_oxd(cal.tim,cal.Ti),100,[0.4 0.4 0.4]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.gbr,cal.oxs,cal.Ti),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.bas,cal.oxs,cal.Ti),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.rhy,cal.oxs,cal.Ti),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
subplot(1,3,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Al),100,[0.4 0.4 0.4]*1.3,'filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Al),100,[0.4 0.4 0.4]*1.0,'filled');
scatter(cal.mem_oxd(cal.tim,cal.Fe),cal.mem_oxd(cal.tim,cal.Al),100,[0.4 0.4 0.4]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.gbr,cal.oxs,cal.Al),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.bas,cal.oxs,cal.Al),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.rhy,cal.oxs,cal.Al),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(1,3,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Mg),100,[0.4 0.4 0.4]*1.3,'filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Mg),100,[0.4 0.4 0.4]*1.0,'filled');
scatter(cal.mem_oxd(cal.tim,cal.Fe),cal.mem_oxd(cal.tim,cal.Mg),100,[0.4 0.4 0.4]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.gbr,cal.oxs,cal.Mg),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.bas,cal.oxs,cal.Mg),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.oxs,cal.Fe),cal.cmp_msy_oxd(cal.rhy,cal.oxs,cal.Mg),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
colorbar;
sgtitle('oxides system',FS{:},TX{:})

Fe=4;
for ic = nc
    MAG(ic).OUT.OxideFract.spnp = zeros(size(MAG(ic).OUT.OxideFract.spn));
    MAG(ic).OUT.OxideFract.spnp(hasoxs,[Fe,Ti,Al,Mg]) = Xp;
end

%% feldspar system
cal_MORB; % load melt model calibration
DATA.PRJCT  = 'BSE';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Al,cal.Ca,cal.Na]);
DATA.SNAMES = {};
DATA.X      = [];
DATA.T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'pl4T')
    hasfsp = MAG(ic).OUT.PhaseProps.pl4T(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.pl4T(hasfsp,[Si,Al,Ca,Na]).*100];
    DATA.T = [DATA.T;MAG(ic).OUT.T(hasfsp)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix

%% plot feldspar system
figure(5); clf;
subplot(2,2,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),100,[0.2 0.1 0.5]*1.3,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),100,[0.2 0.1 0.5]*1.0,'filled');
% scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Al),100,[0.2 0.1 0.5]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Al),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Al),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Al),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),100,[0.2 0.1 0.5]*1.3,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),100,[0.2 0.1 0.5]*1.0,'filled');
% scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Ca),100,[0.2 0.1 0.5]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Ca),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Ca),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Ca),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,2,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T,'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),100,[0.2 0.1 0.5]*1.3,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),100,[0.2 0.1 0.5]*1.0,'filled');
% scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Na),100,[0.2 0.1 0.5]*0.7,'filled');
scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Na),100,cal.T0(cal.gbr),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Na),100,cal.T0(cal.bas),'s','filled');
scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Na),100,cal.T0(cal.rhy),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
colorbar;
sgtitle('felspar system',FS{:},TX{:})

Fe=4;
for ic = nc
    MAG(ic).OUT.OxideFract.pl4Tp = zeros(size(MAG(ic).OUT.OxideFract.pl4T));
    MAG(ic).OUT.OxideFract.pl4Tp(hasfsp,[Si,Al,Ca,Na]) = Xp;
end

%% add up projected mineral compositions to solid composition
for ic = nc
    wt = zeros(size(MAG(ic).OUT.T)) + 1e-16;
    MAG(ic).OUT.OxideFract.solp = zeros(size(MAG(ic).OUT.OxideFract.sol));

    MAG(ic).OUT.OxideFract.solp(hasolv,:) = MAG(ic).OUT.OxideFract.solp(hasolv,:) + MAG(ic).OUT.OxideFract.olp(hasolv,:).*MAG(ic).OUT.PhaseProps.ol(hasolv,1);
    wt(hasolv) = wt(hasolv) + MAG(ic).OUT.PhaseProps.ol(hasolv,1);

    MAG(ic).OUT.OxideFract.solp(hasopx,:) = MAG(ic).OUT.OxideFract.solp(hasopx,:) + MAG(ic).OUT.OxideFract.opxp(hasopx,:).*MAG(ic).OUT.PhaseProps.opx(hasopx,1);
    wt(hasopx) = wt(hasopx) + MAG(ic).OUT.PhaseProps.opx(hasopx,1);

    MAG(ic).OUT.OxideFract.solp(hascpx,:) = MAG(ic).OUT.OxideFract.solp(hascpx,:) + MAG(ic).OUT.OxideFract.cpxp(hascpx,:).*MAG(ic).OUT.PhaseProps.cpx(hascpx,1);
    wt(hascpx) = wt(hascpx) + MAG(ic).OUT.PhaseProps.cpx(hascpx,1);

    MAG(ic).OUT.OxideFract.solp(hasfsp,:) = MAG(ic).OUT.OxideFract.solp(hasfsp,:) + MAG(ic).OUT.OxideFract.pl4Tp(hasfsp,:).*MAG(ic).OUT.PhaseProps.pl4T(hasfsp,1);
    wt(hasfsp) = wt(hasfsp) + MAG(ic).OUT.PhaseProps.pl4T(hasfsp,1);

    MAG(ic).OUT.OxideFract.solp(hasoxs,:) = MAG(ic).OUT.OxideFract.solp(hasoxs,:) + MAG(ic).OUT.OxideFract.spnp(hasoxs,:).*MAG(ic).OUT.PhaseProps.spn(hasoxs,1);
    wt(hasoxs) = wt(hasoxs) + MAG(ic).OUT.PhaseProps.spn(hasoxs,1);

    MAG(ic).OUT.OxideFract.solp = MAG(ic).OUT.OxideFract.solp./wt;
end

%% project melt composition into space of defined mineral end-members
for ic = nc
    hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001;
    oxdliq = MAG(ic).OUT.OxideFract.liq(hasmlt,:)*100;
    oxdliq(:,end-1) = 0;
    oxdliq = oxdliq./sum(oxdliq,2)*100;

    Xp = zeros(length(oxdliq),cal.nmem);
    for ip = 1:length(oxdliq)
        Xp(ip,:) = lsqnonneg(cal.mem_oxd.',oxdliq(ip,:).');
    end

    MAG(ic).OUT.OxideFract.liqp = zeros(size(MAG(ic).OUT.OxideFract.liq));
    MAG(ic).OUT.OxideFract.liqp(hasmlt,:) = Xp*cal.mem_oxd./100;
end

%% reconstruct bulk composition based on projected solid and liquid compositions
for ic = nc
    MAG(ic).OUT.OxideFract.SYSp           = zeros(size(MAG(ic).OUT.OxideFract.SYS));
    MAG(ic).OUT.OxideFract.SYSp(hasmlt,:) = MAG(ic).OUT.OxideFract.liqp(hasmlt,:).*MAG(ic).OUT.PhaseFractions.liq_wt(hasmlt,:) ...
                                          + MAG(ic).OUT.OxideFract.solp(hasmlt,:).*MAG(ic).OUT.PhaseFractions.sol_wt(hasmlt,:);

end

%% liquid, solid, mixture compositions
cal_MORB; Fe = 4; % load melt model calibration
figure(6); clf;
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001;
    subplot(2,3,1);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Ti).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Ti).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Ti).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Ti).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Ti).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Ti).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    scatter(cal.cmp_oxd(cal.ano,cal.Si),cal.cmp_oxd(cal.ano,cal.Ti),140,cal.T0(cal.ano),'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Ti),140,cal.T0(cal.gbr),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Ti),140,cal.T0(3),'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Ti),140,cal.T0(4),'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
    subplot(2,3,2);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    scatter(cal.cmp_oxd(cal.ano,cal.Si),cal.cmp_oxd(cal.ano,cal.Al),140,cal.T0(cal.ano),'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Al),140,cal.T0(cal.gbr),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Al),140,cal.T0(cal.bas),'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Al),140,cal.T0(cal.rhy),'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,3,3);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    scatter(cal.cmp_oxd(cal.ano,cal.Si),cal.cmp_oxd(cal.ano,cal.Fe),140,cal.T0(cal.ano),'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Fe),140,cal.T0(cal.gbr),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Fe),140,cal.T0(cal.bas),'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Fe),140,cal.T0(cal.rhy),'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(2,3,4);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    scatter(cal.cmp_oxd(cal.ano,cal.Si),cal.cmp_oxd(cal.ano,cal.Mg),140,cal.T0(cal.ano),'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Mg),140,cal.T0(cal.gbr),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Mg),140,cal.T0(cal.bas),'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Mg),140,cal.T0(cal.rhy),'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    subplot(2,3,5);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    scatter(cal.cmp_oxd(cal.ano,cal.Si),cal.cmp_oxd(cal.ano,cal.Ca),140,cal.T0(cal.ano),'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Ca),140,cal.T0(cal.gbr),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Ca),140,cal.T0(cal.bas),'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Ca),140,cal.T0(cal.rhy),'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    subplot(2,3,6);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    scatter(cal.cmp_oxd(cal.ano,cal.Si),cal.cmp_oxd(cal.ano,cal.Na),150,cal.T0(cal.ano),'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Na),150,cal.T0(cal.gbr),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Na),150,cal.T0(cal.bas),'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Na),150,cal.T0(cal.rhy),'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('melt, solid, mixture',FS{:},TX{:})


%% best fit cmp_mem
cal_MORB; Fe = 4; % load melt model calibration
oxdSYS = [];
oxdSOL = [];
oxdLIQ = [];
for ic = nc
    hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001;
    oxdSYS = [oxdSYS;MAG(ic).OUT.OxideFract.SYSp(hasmlt,:)*100];
    oxdLIQ = [oxdLIQ;MAG(ic).OUT.OxideFract.liqp(hasmlt,:)*100];
    oxdSOL = [oxdSOL;MAG(ic).OUT.OxideFract.solp(hasmlt,:)*100];
end
oxd0 = [oxdSYS];%;oxdLIQ;oxdSOL];

% DATA.PRJCT  = 'BSE';
% DATA.VNAMES = cal.oxdStr([cal.Si,cal.Ti,cal.Al,cal.Fe,cal.Mg,cal.Ca,cal.Na]);
% DATA.SNAMES = {};
% DATA.X      = oxd0(:,1:end-2);
% DATA.X      = DATA.X./sum(DATA.X,2);
% 
% unmix

cmp_mem = cal.cmp_mem;
cmp_mem_best = cmp_mem;
cmp_oxd = cal.cmp_oxd;

Xp = zeros(length(oxd0),cal.ncmp);
for ip = 1:length(oxd0)
    Xp(ip,:) = lsqnonneg(cmp_oxd.',oxd0(ip,:).');
end
Xp = Xp./sum(Xp,2);
oxdfit = Xp*cmp_oxd;

misfit  = norm((oxdfit-oxd0)./(oxd0+0.01))/sqrt(length(oxd0(:)));
bestfit = misfit*1.001;
dfit    = misfit*0.001;

figure(100); clf;

ifit = [2,3,4];
tol = 1e-2; it = 1;
while misfit>tol && dfit>1e-6 && it<1e5

    cmp_mem(ifit,1:end-1) = max(1e-6,cmp_mem_best(ifit,1:end-1) .* (1 + randn(size(cmp_mem(ifit,1:end-1))).*min(1e-1,bestfit^0.1/500)));
    cmp_mem(ifit,:) = cmp_mem(ifit,:)./sum(cmp_mem(ifit,:),2)*100;
    cmp_oxd = cmp_mem*cal.mem_oxd./100;

    Xp = zeros(length(oxd0),cal.ncmp);
    for ip = 1:length(oxd0)
        Xp(ip,:) = lsqnonneg(cmp_oxd.',oxd0(ip,:).');
    end
    Xp = Xp./sum(Xp,2);
    oxdfit = Xp*cmp_oxd;

    misfit = norm((oxdfit-oxd0)./(oxd0+0.01))/sqrt(length(oxd0(:)));

    if misfit<1.0*bestfit
        dfit = abs(bestfit-misfit)*0.1 + dfit*0.9;
        bestfit = misfit;
        cmp_mem_best = cmp_mem;
        cmp_oxd_best = cmp_oxd;
        cfit_best    = cfit;
        oxdfit_best  = oxdfit;
        subplot(2,1,1); loglog(it,(bestfit),'k.'); hold on; axis tight; box on;
        subplot(2,1,2); loglog(it,(dfit   ),'k.'); hold on; axis tight; box on;
        drawnow;
    end

    it = it+1;
end



%% best fit melting temperatures
% clear cal 
clear cal var x m f cx cm
cal_MORB;  % load melt model calibration
DATA.c_oxd  = [];
DATA.cm_oxd = [];
DATA.cx_oxd = [];
DATA.T      = [];
DATA.P      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'pl4T')
        hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001;
        DATA.c_oxd  = [DATA.c_oxd ;MAG(ic).OUT.OxideFract.SYSp(hasmlt,:).*100];
        DATA.cm_oxd = [DATA.cm_oxd;MAG(ic).OUT.OxideFract.liqp(hasmlt,:).*100];
        DATA.cx_oxd = [DATA.cx_oxd;MAG(ic).OUT.OxideFract.solp(hasmlt,:).*100];
        DATA.T = [DATA.T;MAG(ic).OUT.T(hasmlt)];
        DATA.P = [DATA.P;MAG(ic).OUT.P(hasmlt)*1e8];
    end
end
Xp = zeros(length(DATA.c_oxd),cal.ncmp);
for ip = 1:length(DATA.c_oxd)
    Xp(ip,:) = lsqnonneg(cal.cmp_oxd.',DATA.c_oxd(ip,:).');
end
c = Xp./sum(Xp,2);
T = DATA.T;
P = DATA.P;

% equilibrium phase fractions and compositions
Nz = length(T); Nx = 1;

figure(100); clf;

cal_best = cal;
misfit = 15; bestfit = 16; dfit = 0.01; tol = 1e-3; it = 1;
while misfit>tol && dfit>1e-5 && it<1e6

    % set pure component melting points T_m^i at P=0
    cal.T0(2:end) = max(750,cal_best.T0(2:end) .* (1 + randn(size(cal.T0(2:end))).*min(1e-2,bestfit^0.1/3e3)));
    cal.r         = max(  5,cal_best.r         .* (1 + randn(size(cal.r        )).*min(1e-2,bestfit^0.1/3e3)));
    cal.A         = (cal.T0+273.15)./300;

    var.m = 0.5; var.x = 0; var.f = 0; cm = 0.*c; cx = 0.*c; var.cm = c(1,:); var.cx = c(1,:); cm_oxd = c(1,:)*cal.cmp_oxd;
    for i=1:length(T)
        % update local phase equilibrium
        var.c      = c(i,:);        % component fractions [wt]
        var.T      = T(i);          % temperature [C]
        var.P      = P(i)/1e9;      % pressure [GPa]
        var.H2O    = c(i,end);      % water concentration [wt]
        var.SiO2m  = cm_oxd(:,1)./sum(cm_oxd(:,1:end-1)); % melt silica concentration [wt]
        cal.H2Osat = fluidsat(var.T,var.P*1e9,var.SiO2m,cal);
        [var,cal]  = meltmodel(var,cal,'E');

        m(i) = var.m;
        f(i) = var.f;
        x(i) = var.x;

        cx(i,:) = var.cx;
        cm(i,:) = var.cm;
        cm_oxd  = var.cm*cal.cmp_oxd;
    end

    c_oxd  = reshape(reshape( c,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
    cm_oxd = reshape(reshape(cm,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);
    cx_oxd = reshape(reshape(cx,Nz*Nx,cal.ncmp)*cal.cmp_oxd,Nz,Nx,cal.noxd);

    misfit = norm((squeeze(cm_oxd(:,1:end-1))-DATA.cm_oxd(:,1:end-1))./(squeeze(cm_oxd(:,1:end-1))+0.01)) + norm((squeeze(cx_oxd(:,1:end-1))-DATA.cx_oxd(:,1:end-1))./(squeeze(cx_oxd(:,1:end-1))+0.01))/10;

    if misfit<1.0*bestfit
        dfit = abs(bestfit-misfit)*0.1 + dfit*0.9;
        bestfit = misfit;
        cal_best = cal;
        cm_oxd_best = cm_oxd;
        cx_oxd_best = cx_oxd;
        subplot(2,1,1); loglog(it,(bestfit),'k.'); hold on; axis tight; box on;
        subplot(2,1,2); loglog(it,(dfit   ),'k.'); hold on; axis tight; box on;
        drawnow;
    end

    it = it+1;
end



%% plot phase diagram
figure(7); clf;
subplot(3,3,1)
plot(DATA.cm_oxd(:,cal.Si)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Si)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Si)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Si)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Si)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Si)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Si},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,2)
plot(DATA.cm_oxd(:,cal.Ti)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Ti)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Ti)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Ti)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Ti)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Ti)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Ti},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,3)
plot(DATA.cm_oxd(:,cal.Al)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Al)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Al)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Al)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Al)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Al)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Al},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,4)
plot(DATA.cm_oxd(:,cal.Fe)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Fe)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Fe)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Fe)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Fe)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Fe)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Fe},' [wt]'],'Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

subplot(3,3,5)
plot(DATA.cm_oxd(:,cal.Mg)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Mg)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Mg)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Mg)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Mg)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Mg)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Mg},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,6)
plot(DATA.cm_oxd(:,cal.Ca)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Ca)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Ca)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Ca)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Ca)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Ca)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Ca},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,7)
plot(DATA.cm_oxd(:,cal.Na)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Na)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Na)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Na)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.Na)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.Na)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Na},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,8)
plot(DATA.cm_oxd(:,cal.K)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.K)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.K)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.K)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.K)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.K)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.K},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,9)
plot(DATA.cm_oxd(:,cal.H)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.H)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.H)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.H)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_best(:,cal.H)./sum(cx_oxd_best(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_best(:,cal.H)./sum(cm_oxd_best(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.H},' [wt]'],'Interpreter','latex','FontSize',15)

% sgtitle('Phase Diagram','Interpreter','latex','FontSize',18)
% xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
% ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)


%% update material closures
Nz = length(T); Nx = 1; Ptop = min(P); Pt = P; etareg = 1; calibrt = 1; T = T+273.15;
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
