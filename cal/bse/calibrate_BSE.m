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

run('../../usr/par_default');  % load default parameters

% calibration run options
runID     = 'cal_BSE';           % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set phase diagram parameters
cal_BSE;  % load melt model calibration

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

if ~holdfig; close all; end

%% load MAGEMin results
nc  = [1]; % number of compositions modelled
MAG = [];
ioxd = [1 8 2 5 4 3 7 6]; % oxide indices from MAGEMin to standard
Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; K = 8;
for ic = nc
    filename = ['BSE_fract5_anh_out.mat'];
    load(filename);

    % lump in free O to FeO, Cr2O3 to Al2O3, normalise to anhydrous unit sum
    phs = fieldnames(OUT.PhaseProps);
    phs = [phs(:)',{'SYS'},{'sol'}];
    for iph = 1:length(phs)
        OUT.OxideFract.(phs{iph}) = zeros(size(OUT.OxideFractions.(phs{iph})));
        OUT.OxideFractions.(phs{iph})(:,5) = OUT.OxideFractions.(phs{iph})(:,5) + OUT.OxideFractions.(phs{iph})(:,9);
        OUT.OxideFractions.(phs{iph})(:,2) = OUT.OxideFractions.(phs{iph})(:,2) + OUT.OxideFractions.(phs{iph})(:,10);
        OUT.OxideFract.(phs{iph}) = OUT.OxideFractions.(phs{iph})(:,[ioxd 11]);
        OUT.OxideFract.(phs{iph})(:,1:cal.noxd) = OUT.OxideFract.(phs{iph})(:,1:cal.noxd)./sum(OUT.OxideFract.(phs{iph})(:,1:cal.noxd)+1e-16,2);
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


    % combine trd with qtz
    if isfield(OUT.PhaseProps,'trd')
        OUT.OxideFract.q = (OUT.OxideFract.q.*OUT.PhaseProps.q(:,1) + OUT.OxideFract.trd.*OUT.PhaseProps.trd(:,1)) ./ (OUT.PhaseProps.q(:,1)+OUT.PhaseProps.trd(:,1)+1e-16);
        OUT.PhaseProps.q = OUT.PhaseProps.q(:,1)+OUT.PhaseProps.trd(:,1);
        OUT.OxideFract = rmfield(OUT.OxideFract,'trd');
        OUT.PhaseProps = rmfield(OUT.PhaseProps,'trd');
    end

    MAG(ic).OUT = OUT;
end

% calibrate mineral end-members

%% olivine system
cal_BSE;  % load melt model calibration
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
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.dun,cal.Si),cal.mem_oxd(cal.dun,cal.Fe),100,1890,'filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Fe),100,cal.T0(cal.dun),'s','filled');
scatter(cal.cmp_msy_oxd(cal.fay,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.fay,cal.olv,cal.Fe),100,cal.T0(cal.fay),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(1,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.dun,cal.Si),cal.mem_oxd(cal.dun,cal.Mg),100,1890,'filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Mg),100,cal.T0(cal.dun),'s','filled');
scatter(cal.cmp_msy_oxd(cal.fay,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.fay,cal.olv,cal.Mg),100,cal.T0(cal.fay),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
colorbar;
sgtitle('olivine system',FS{:},TX{:})


%% orthopyroxene system
cal_BSE;  Fe=4; % load melt model calibration
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
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Al),100,1250,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),100,1030,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Al),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Al),100,cal.T0(2),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Al),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Al),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Fe),100,1250,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),100,1030,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Fe),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Fe),100,cal.T0(2),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Fe),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Fe),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,2,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Mg),100,1250,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),100,1030,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Mg),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Mg),100,cal.T0(2),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Mg),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Mg),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,2,4);
scatter(X(:,1)*100,X(:,5)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,5)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Ca),100,1250,'filled');
scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),100,1030,'filled');
scatter(cal.mem_oxd(cal.fsl,cal.Si),cal.mem_oxd(cal.fsl,cal.Ca),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.pxn,cal.opx,cal.Ca),100,cal.T0(2),'s','filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Ca),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.opx,cal.Ca),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
colorbar;
sgtitle('orthopyroxene system',FS{:},TX{:})


%% clinopyroxene system
cal_BSE; Fe = 4; % load melt model calibration
DATA.PRJCT  = 'BSE';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Al,cal.Fe,cal.Mg,cal.Ca,cal.Na]);
DATA.SNAMES = {};
DATA.X      = [];
DATA.T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hascpx = MAG(ic).OUT.PhaseProps.cpx(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.cpx(hascpx,[Si,Al,Fe,Mg,Ca,Na]).*100];
    DATA.T = [DATA.T;MAG(ic).OUT.T(hascpx)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);
unmix

%%
figure(4); clf;
subplot(2,3,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Al),100,1090,'filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Al),100, 950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Al),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Al),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Fe),100,1090,'filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Fe),100, 950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Fe),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Fe),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Mg),100,1090,'filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Mg),100, 950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Mg),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Mg),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,4);
scatter(X(:,1)*100,X(:,5)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,5)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ca),100,1090,'filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Ca),100, 950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Ca),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Ca),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,5);
scatter(X(:,1)*100,X(:,6)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,6)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Na),100,1090,'filled');
scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Na),100, 950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Na),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.cpx,cal.Na),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
colorbar;
sgtitle('clinopyroxene system',FS{:},TX{:})


%% feldspar system
cal_BSE; Fe = 4; % load melt model calibration
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
subplot(1,3,1);
scatter(X(:,1)*100,X(:,2)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),100,1150,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Al),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.fsp,cal.Al),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(1,3,2);
scatter(X(:,1)*100,X(:,3)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),100,1150,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Ca),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.fsp,cal.Ca),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(1,3,3);
scatter(X(:,1)*100,X(:,4)*100,25,DATA.T,'filled'); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,DATA.T);
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),100,1150,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),100,950,'filled');
scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Na),100,cal.T0(3),'s','filled');
scatter(cal.cmp_msy_oxd(cal.eut,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.eut,cal.fsp,cal.Na),100,cal.T0(4),'s','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
colorbar;
sgtitle('felspar system',FS{:},TX{:})


%% liquid, solid, mixture compositions
cal_BSE; Fe = 4; % load melt model calibration
figure(6); clf;
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001 & MAG(ic).OUT.T>=970;
    subplot(2,3,1);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'d'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Al),140,cal.T0(1),'filled','o');
    scatter(cal.cmp_oxd(cal.pxn,cal.Si),cal.cmp_oxd(cal.pxn,cal.Al),140,cal.T0(2),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Al),140,cal.T0(3),'filled','o');
    scatter(cal.cmp_oxd(cal.eut,cal.Si),cal.cmp_oxd(cal.eut,cal.Al),140,cal.T0(4),'filled','o');
%     scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Al),100,[0.7 0.2 0.1]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),100,[0.7 0.2 0.1]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Al),100,[0.6 0.1 0.4]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Al),100,[0.6 0.1 0.4]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),100,[0.1 0.2 0.7]*0.8,'filled','s');
%     scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),100,[0.1 0.2 0.7]*1.0,'filled','s');
%     scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Al),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,3,2);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Fe).*100,25,MAG(ic).OUT.T(hasmlt),'d'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Fe),140,cal.T0(1),'filled','o');
    scatter(cal.cmp_oxd(cal.pxn,cal.Si),cal.cmp_oxd(cal.pxn,cal.Fe),140,cal.T0(2),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Fe),140,cal.T0(3),'filled','o');
    scatter(cal.cmp_oxd(cal.eut,cal.Si),cal.cmp_oxd(cal.eut,cal.Fe),140,cal.T0(4),'filled','o');
%     scatter(cal.mem_oxd(cal.dun,cal.Si),cal.mem_oxd(cal.dun,cal.Fe),100,[0.3 0.6 0.1]*0.8,'filled','^');
%     scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),100,[0.3 0.6 0.1]*1.4,'filled','^');
%     scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Fe),100,[0.7 0.2 0.1]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),100,[0.7 0.2 0.1]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Fe),100,[0.6 0.1 0.4]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Fe),100,[0.6 0.1 0.4]*1.4,'filled','d');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(2,3,3);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'d'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Mg),140,cal.T0(1),'filled','o');
    scatter(cal.cmp_oxd(cal.pxn,cal.Si),cal.cmp_oxd(cal.pxn,cal.Mg),140,cal.T0(2),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Mg),140,cal.T0(3),'filled','o');
    scatter(cal.cmp_oxd(cal.eut,cal.Si),cal.cmp_oxd(cal.eut,cal.Mg),140,cal.T0(4),'filled','o');
%     scatter(cal.mem_oxd(cal.dun,cal.Si),cal.mem_oxd(cal.dun,cal.Mg),100,[0.3 0.6 0.1]*0.8,'filled','^');
%     scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),100,[0.3 0.6 0.1]*1.4,'filled','^');
%     scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Mg),100,[0.7 0.2 0.1]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),100,[0.7 0.2 0.1]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Mg),100,[0.6 0.1 0.4]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Mg),100,[0.6 0.1 0.4]*1.4,'filled','d');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    subplot(2,3,4);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'d'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Ca),140,cal.T0(1),'filled','o');
    scatter(cal.cmp_oxd(cal.pxn,cal.Si),cal.cmp_oxd(cal.pxn,cal.Ca),140,cal.T0(2),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Ca),140,cal.T0(3),'filled','o');
    scatter(cal.cmp_oxd(cal.eut,cal.Si),cal.cmp_oxd(cal.eut,cal.Ca),140,cal.T0(4),'filled','o');
%     scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Ca),100,[0.7 0.2 0.1]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),100,[0.7 0.2 0.1]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ca),100,[0.6 0.1 0.4]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Ca),100,[0.6 0.1 0.4]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),100,[0.1 0.2 0.7]*0.8,'filled','s');
%     scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),100,[0.1 0.2 0.7]*1.0,'filled','s');
%     scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Ca),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    subplot(2,3,5);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'d'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Na),150,cal.T0(1),'filled','o');
    scatter(cal.cmp_oxd(cal.pxn,cal.Si),cal.cmp_oxd(cal.pxn,cal.Na),150,cal.T0(2),'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Na),150,cal.T0(3),'filled','o');
    scatter(cal.cmp_oxd(cal.eut,cal.Si),cal.cmp_oxd(cal.eut,cal.Na),150,cal.T0(4),'filled','o');
%     scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Na),100,[0.6 0.1 0.4]*0.8,'filled','d');
%     scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Na),100,[0.6 0.1 0.4]*1.4,'filled','d');
%     scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),100,[0.1 0.2 0.7]*0.8,'filled','s');
%     scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),100,[0.1 0.2 0.7]*1.0,'filled','s');
%     scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Na),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('melt, solid, mixture',FS{:},TX{:})

%% best fit cmp_mem
cal_BSE; Fe = 4; % load melt model calibration
oxdSYS      = MAG(1).OUT.OxideFract.SYS(hasmlt,[Si Al Fe Mg Ca Na])./sum(MAG(1).OUT.OxideFract.SYS(hasmlt,[Si Al Fe Mg Ca Na]),2)*100;
oxdLIQ      = MAG(1).OUT.OxideFract.liq(hasmlt,[Si Al Fe Mg Ca Na])./sum(MAG(1).OUT.OxideFract.liq(hasmlt,[Si Al Fe Mg Ca Na]),2)*100;
oxdSOL      = MAG(1).OUT.OxideFract.sol(hasmlt,[Si Al Fe Mg Ca Na])./sum(MAG(1).OUT.OxideFract.sol(hasmlt,[Si Al Fe Mg Ca Na]),2)*100;
oxd0 = [oxdSYS;oxdLIQ];%;oxdSOL];

DATA.PRJCT  = 'BSE';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Al,cal.Fe,cal.Mg,cal.Ca,cal.Na]);
DATA.SNAMES = {};
DATA.X      = oxd0;
DATA.X      = DATA.X./sum(DATA.X,2);
unmix

oxd0(:,7) = 0.0;
cmp_mem = cal.cmp_mem;
cmp_mem_best = cmp_mem;
cmp_oxd = cal.cmp_oxd;

figure(100); clf;

ifit = [2,3,4];
misfit = 3; bestfit = 3.1; dfit = 0.01; tol = 1e-3; it = 1;
while misfit>tol && dfit>1e-5 && it<1e6

    cmp_mem(ifit,1:end-1) = max(1e-9,cmp_mem_best(ifit,1:end-1) .* (1 + randn(size(cmp_mem(ifit,1:end-1))).*min(1e-1,bestfit^0.25/20)));
    cmp_mem(ifit,:) = cmp_mem(ifit,:)./sum(cmp_mem(ifit,:),2)*100;
    cmp_oxd = cmp_mem*cal.mem_oxd./100;

    cfit   = max(0.0,cmp_oxd.'\oxd0.').'; cfit(:,end) = 0; cfit = cfit./sum(cfit,2);
    oxdfit = cfit*cmp_oxd;

    misfit = norm((oxdfit-oxd0)./(oxd0+0.1))/sqrt(length(oxd0(:)));

    if misfit<1.001*bestfit
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
cal_BSE;  % load melt model calibration
DATA.c_oxd  = [];
DATA.cm_oxd = [];
DATA.cx_oxd = [];
DATA.T      = [];
DATA.P      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'pl4T')
        hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001 & MAG(ic).OUT.T>=970;
        DATA.c_oxd  = [DATA.c_oxd ;MAG(ic).OUT.OxideFract.SYS(hasmlt,[Si Al Fe Mg Ca Na end]).*100];
        DATA.cm_oxd = [DATA.cm_oxd;MAG(ic).OUT.OxideFract.liq(hasmlt,[Si Al Fe Mg Ca Na end]).*100];
        DATA.cx_oxd = [DATA.cx_oxd;MAG(ic).OUT.OxideFract.sol(hasmlt,[Si Al Fe Mg Ca Na end]).*100];
        DATA.T = [DATA.T;MAG(ic).OUT.T(hasmlt)];
        DATA.P = [DATA.P;MAG(ic).OUT.P(hasmlt)*1e8];
    end
end
c      = max(1e-3,cal.cmp_oxd.'\DATA.c_oxd.').'; c(:,end) = 0; c = c./sum(c,2);
oxdfit = c*cal.cmp_oxd;
T = DATA.T;
P = DATA.P;

% equilibrium phase fractions and compositions
Nz = length(T); Nx = 1;

figure(100); clf;

% % set pure component melting points T_m^i at P=0
% cal.T0(cal.mdu) =  1553;
% cal.T0(cal.bas) =  1190;
% cal.T0(cal.bas) =  980;
% cal.T0(cal.eut) =  850;
% 
% % set coeff. for T-dependence of partition coefficients K^i [1/K]
% cal.r(cal.mdu)  =  18.0;
% cal.r(cal.bas)  =  16.0;
% cal.r(cal.bas)  =  14.0;
% cal.r(cal.eut)  =  12.0;
cal_best = cal;
misfit = 15; bestfit = 16; dfit = 0.01; tol = 1e-3; it = 1;
while misfit>tol && dfit>1e-5 && it<1e6

    % set pure component melting points T_m^i at P=0
    cal.T0(2:end) = max(800,cal_best.T0(2:end) .* (1 + randn(size(cal.T0(2:end))).*min(1e-2,bestfit^0.25/1e3)));
    cal.r         = max(  5,cal_best.r         .* (1 + randn(size(cal.r        )).*min(1e-2,bestfit^0.25/1e3)));

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

    misfit = norm((squeeze(cm_oxd(:,1:end-1))-DATA.cm_oxd(:,1:end-1))./(squeeze(cm_oxd(:,1:end-1))+1)) + norm((squeeze(cx_oxd)-DATA.cx_oxd)./(squeeze(cx_oxd)+10))/10;

    if misfit<1.001*bestfit
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
figure(7);
subplot(3,3,1)
plot(DATA.cm_oxd(:,cal.Si)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Si)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Si)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Si)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.Si)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.Si)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Si},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,3)
plot(DATA.cm_oxd(:,cal.Al)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Al)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Al)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Al)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.Al)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.Al)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Al},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,4)
plot(DATA.cm_oxd(:,cal.Fe)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Fe)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Fe)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Fe)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.Fe)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.Fe)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Fe},' [wt]'],'Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

subplot(3,3,5)
plot(DATA.cm_oxd(:,cal.Mg)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Mg)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Mg)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Mg)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.Mg)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.Mg)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Mg},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,6)
plot(DATA.cm_oxd(:,cal.Ca)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Ca)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Ca)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Ca)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.Ca)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.Ca)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Ca},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,7)
plot(DATA.cm_oxd(:,cal.Na)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.Na)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.Na)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.Na)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.Na)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.Na)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Na},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,9)
plot(DATA.cm_oxd(:,cal.H)./sum(DATA.cm_oxd(:,1:end-1),2),DATA.T,'ro'); axis tight; hold on; box on;
plot(DATA.cx_oxd(:,cal.H)./sum(DATA.cx_oxd(:,1:end-1),2),DATA.T,'bs');
plot(DATA.c_oxd (:,cal.H)./sum(DATA.c_oxd (:,1:end-1),2),DATA.T,'kd');

plot( c_oxd(:,cal.H)./sum( c_oxd(:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd(:,cal.H)./sum(cx_oxd(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd(:,cal.H)./sum(cm_oxd(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
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
