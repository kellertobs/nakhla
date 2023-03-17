% calibrate phase diagram
clear all;

addpath('../../cal')
addpath('../../src')
load ocean
TINY = 1e-16;
FS = {'FontSize',15};
TX = {'Interpreter','latex'};

run('../../usr/par_default');  % load default parameters

% calibration run options
runID     = 'cal_andes';         % run ID for output files; [system name_wt.% SiO2_wt.% H2O] 
holdfig   = 0;                   % set to 1 to hold figures, to 0 for new figures
linestyle = '-';                 % set line style for plots
save_plot = 0;                   % turn on (1) to save output file in /out directory

% set phase diagram parameters
cal_andes;  % load melt model calibration

% set model buoyancy parameters
dx       =  1e-3;                % crystal size [m]
df       =  1e-3;                % bubble size [m]
g0       =  10.;                 % gravity [m/s2]

if ~holdfig; close all; end

%% load MAGEMin results
nc  = [1 2 3]; % number of cluster compositions modelled
MAG = [];
ioxd = [1 8 2 5 4 3 7 6]; % oxide indices from MAGEMin to standard
for ic = nc
    filename = ['andes_Fc',int2str(ic),'_fract_anh_out.mat'];
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

%     % lump in ilmenite with spinel
%     if isfield(OUT.PhaseProps,'ilm')
%         OUT.OxideFract.spn = (OUT.OxideFract.spn.*OUT.PhaseProps.spn(:,1) + OUT.OxideFract.ilm.*OUT.PhaseProps.ilm(:,1)) ./ (OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.ilm(:,1)+1e-16);
%         OUT.PhaseProps.spn = OUT.PhaseProps.spn(:,1)+OUT.PhaseProps.ilm(:,1);
%         OUT.OxideFract = rmfield(OUT.OxideFract,'ilm');
%         OUT.PhaseProps = rmfield(OUT.PhaseProps,'ilm');
%     end
% 
%     % lump in rutile with quartz
%     if isfield(OUT.PhaseProps,'ru') && isfield(OUT.PhaseProps,'q')
%         OUT.OxideFract.q = (OUT.OxideFract.q.*OUT.PhaseProps.q(:,1) + OUT.OxideFract.ru.*OUT.PhaseProps.ru(:,1)) ./ (OUT.PhaseProps.q(:,1)+OUT.PhaseProps.ru(:,1)+1e-16);
%         OUT.PhaseProps.q = OUT.PhaseProps.q(:,1)+OUT.PhaseProps.q2(:,1);
%         OUT.OxideFract = rmfield(OUT.OxideFract,'ru');
%         OUT.PhaseProps = rmfield(OUT.PhaseProps,'ru');
%     end
% 
%     % lump in sillimanite with feldspar
%     if isfield(OUT.PhaseProps,'sill')
%         OUT.OxideFract.pl4T = (OUT.OxideFract.pl4T.*OUT.PhaseProps.pl4T(:,1) + OUT.OxideFract.sill.*OUT.PhaseProps.sill(:,1)) ./ (OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.sill(:,1)+1e-16);
%         OUT.PhaseProps.pl4T = OUT.PhaseProps.pl4T(:,1)+OUT.PhaseProps.sill(:,1);
%         OUT.OxideFract = rmfield(OUT.OxideFract,'sill');
%         OUT.PhaseProps = rmfield(OUT.PhaseProps,'sill');
%     end

    MAG(ic).OUT = OUT;
end

% calibrate mineral end-members

%% olivine system
cal_andes;  % load melt model calibration
figure(1); clf;
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'ol')
    hasolv = MAG(ic).OUT.PhaseProps.ol(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    subplot(1,2,1);
    scatter(MAG(ic).OUT.OxideFract.ol(hasolv,cal.Si).*100,MAG(ic).OUT.OxideFract.ol(hasolv,cal.Fe).*100,25,MAG(ic).OUT.T(hasolv)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Fe),100,1890,'filled');
    scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Fe),100,1850,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Fe),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Fe),100,1050,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(1,2,2);
    scatter(MAG(ic).OUT.OxideFract.ol(hasolv,cal.Si).*100,MAG(ic).OUT.OxideFract.ol(hasolv,cal.Mg).*100,25,MAG(ic).OUT.T(hasolv)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Mg),100,1890,'filled');
    scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.dun,cal.olv,cal.Mg),100,1850,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.olv,cal.Mg),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.olv,cal.Mg),100,1050,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('olivine system',FS{:},TX{:})

%% spinel system
cal_andes;  % load melt model calibration
figure(2); clf;
for ic = nc(2:end)
    if isfield(MAG(ic).OUT.PhaseProps,'spn')
    hasspn = MAG(ic).OUT.PhaseProps.spn(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    subplot(1,3,1);
    scatter(MAG(ic).OUT.OxideFract.spn(hasspn,cal.Fe).*100,MAG(ic).OUT.OxideFract.spn(hasspn,cal.Ti).*100,25,MAG(ic).OUT.T(hasspn)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Ti),100,860,'filled');
    scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Ti),100,1160,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.gbr,cal.spn,cal.Ti),100,1160,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.bas,cal.spn,cal.Ti),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.rhy,cal.spn,cal.Ti),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})
    subplot(1,3,2);
    scatter(MAG(ic).OUT.OxideFract.spn(hasspn,cal.Fe).*100,MAG(ic).OUT.OxideFract.spn(hasspn,cal.Al).*100,25,MAG(ic).OUT.T(hasspn)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Al),100,860,'filled');
    scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Al),100,1160,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.gbr,cal.spn,cal.Al),100,1160,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.bas,cal.spn,cal.Al),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.rhy,cal.spn,cal.Al),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(1,3,3);
    scatter(MAG(ic).OUT.OxideFract.spn(hasspn,cal.Fe).*100,MAG(ic).OUT.OxideFract.spn(hasspn,cal.Mg).*100,25,MAG(ic).OUT.T(hasspn)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Mg),100,860,'filled');
    scatter(cal.mem_oxd(cal.ulv,cal.Fe),cal.mem_oxd(cal.ulv,cal.Mg),100,1160,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.gbr,cal.spn,cal.Mg),100,1160,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.bas,cal.spn,cal.Mg),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.spn,cal.Fe),cal.cmp_msy_oxd(cal.rhy,cal.spn,cal.Mg),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('spinel system',FS{:},TX{:})


%% orthopyroxene system
cal_andes;  % load melt model calibration
figure(3); clf;
for ic = nc(2:end)
    if isfield(MAG(ic).OUT.PhaseProps,'opx')
    hasopx = MAG(ic).OUT.PhaseProps.opx(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    subplot(2,2,1);
    scatter(MAG(ic).OUT.OxideFract.opx(hasopx,cal.Si).*100,MAG(ic).OUT.OxideFract.opx(hasopx,cal.Al).*100,25,MAG(ic).OUT.T(hasopx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Al),100,1250,'filled');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Al),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Al),100,1050,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,2,2);
    scatter(MAG(ic).OUT.OxideFract.opx(hasopx,cal.Si).*100,MAG(ic).OUT.OxideFract.opx(hasopx,cal.Fe).*100,25,MAG(ic).OUT.T(hasopx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Fe),100,1250,'filled');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Fe),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Fe),100,1050,'s','filled');
    xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(2,2,3);
    scatter(MAG(ic).OUT.OxideFract.opx(hasopx,cal.Si).*100,MAG(ic).OUT.OxideFract.opx(hasopx,cal.Mg).*100,25,MAG(ic).OUT.T(hasopx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Mg),100,1250,'filled');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Mg),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Mg),100,1050,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    subplot(2,2,4);
    scatter(MAG(ic).OUT.OxideFract.opx(hasopx,cal.Si).*100,MAG(ic).OUT.OxideFract.opx(hasopx,cal.Ca).*100,25,MAG(ic).OUT.T(hasopx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Ca),100,1250,'filled');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.opx,cal.Ca),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.opx,cal.Ca),100,1050,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('orthopyroxene system',FS{:},TX{:})


%% clinopyroxene system
cal_andes;  % load melt model calibration
figure(4); clf;
for ic = nc(2:end)
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hascpx = MAG(ic).OUT.PhaseProps.cpx(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    subplot(2,3,1);
    scatter(MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Si).*100,MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Al).*100,25,MAG(ic).OUT.T(hascpx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Al),100,1200,'filled');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Al),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Al),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Al),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Al),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,3,2);
    scatter(MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Si).*100,MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Fe).*100,25,MAG(ic).OUT.T(hascpx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Fe),100,1200,'filled');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Fe),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Fe),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Fe),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Fe),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(2,3,3);
    scatter(MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Si).*100,MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Mg).*100,25,MAG(ic).OUT.T(hascpx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Mg),100,1200,'filled');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Mg),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Mg),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Mg),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Mg),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    subplot(2,3,4);
    scatter(MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Si).*100,MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Ca).*100,25,MAG(ic).OUT.T(hascpx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ca),100,1200,'filled');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Ca),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Ca),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Ca),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Ca),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    subplot(2,3,5);
    scatter(MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Si).*100,MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Na).*100,25,MAG(ic).OUT.T(hascpx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Na),100,1200,'filled');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Na),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Na),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Na),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Na),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
    subplot(2,3,6);
    scatter(MAG(ic).OUT.OxideFract.cpx(hascpx,cal.Si).*100,MAG(ic).OUT.OxideFract.cpx(hascpx,cal.K).*100,25,MAG(ic).OUT.T(hascpx)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.K),100,1200,'filled');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.K),100,900,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.cpx,cal.K),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.cpx,cal.K),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.cpx,cal.K),100, 850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.K),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('clinopyroxene system',FS{:},TX{:})


%% feldspar system
cal_andes;  % load melt model calibration
figure(5); clf;
for ic = nc(2:end)
    if isfield(MAG(ic).OUT.PhaseProps,'pl4T')
    hasfsp = MAG(ic).OUT.PhaseProps.pl4T(:,1)>0.001 & MAG(ic).OUT.PhaseFractions.liq_wt>=0.001;
    subplot(2,2,1);
    scatter(MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Si).*100,MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Al).*100,25,MAG(ic).OUT.T(hasfsp)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),100,1250,'filled');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),100,1050,'filled');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Al),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Al),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Al),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Al),100,850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,2,2);
    scatter(MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Si).*100,MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Ca).*100,25,MAG(ic).OUT.T(hasfsp)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),100,1250,'filled');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),100,1050,'filled');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Ca),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Ca),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Ca),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Ca),100,850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    subplot(2,2,3);
    scatter(MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Si).*100,MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Na).*100,25,MAG(ic).OUT.T(hasfsp)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),100,1250,'filled');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),100,1050,'filled');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Na),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Na),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Na),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Na),100,850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
    subplot(2,2,4);
    scatter(MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.Si).*100,MAG(ic).OUT.OxideFract.pl4T(hasfsp,cal.K).*100,25,MAG(ic).OUT.T(hasfsp)); colormap('copper'); hold on
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.K),100,1250,'filled');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.K),100,1050,'filled');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.K),100,850,'filled');
    scatter(cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.gbr,cal.fsp,cal.K),100,1250,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.bas,cal.fsp,cal.K),100,1050,'s','filled');
    scatter(cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.Si),cal.cmp_msy_oxd(cal.rhy,cal.fsp,cal.K),100,850,'s','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.K),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('felspar system',FS{:},TX{:})


%% liquid, solid, mixture compositions
cal_andes;  % load melt model calibration
figure(6); clf;
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001;
    subplot(2,4,1);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Ti).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Ti).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Ti),140,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Ti),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Ti),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Ti),140,cal.Tphs0,'filled','o');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,4,2);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Al).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Al).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Al),140,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Al),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Al),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Al),140,cal.Tphs0,'filled','o');
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Al),100,[0.7 0.2 0.1]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Al),100,[0.7 0.2 0.1]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Al),100,[0.6 0.1 0.4]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Al),100,[0.6 0.1 0.4]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),100,[0.1 0.2 0.7]*0.8,'filled','s');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),100,[0.1 0.2 0.7]*1.0,'filled','s');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Al),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,4,3);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Fe).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Fe).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Fe),140,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Fe),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Fe),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Fe),140,cal.Tphs0,'filled','o');
%     scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Fe),100,[0.3 0.6 0.1]*0.8,'filled','^');
    scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),100,[0.3 0.6 0.1]*1.4,'filled','^');
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Fe),100,[0.7 0.2 0.1]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Fe),100,[0.7 0.2 0.1]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Fe),100,[0.6 0.1 0.4]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Fe),100,[0.6 0.1 0.4]*1.4,'filled','d');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(2,4,4);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Mg).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Mg).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Mg),140,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Mg),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Mg),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Mg),140,cal.Tphs0,'filled','o');
    scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Mg),100,[0.3 0.6 0.1]*0.8,'filled','^');
%     scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),100,[0.3 0.6 0.1]*1.4,'filled','^');
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Mg),100,[0.7 0.2 0.1]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Mg),100,[0.7 0.2 0.1]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Mg),100,[0.6 0.1 0.4]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Mg),100,[0.6 0.1 0.4]*1.4,'filled','d');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    subplot(2,4,5);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Ca).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Ca).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Ca),140,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Ca),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Ca),140,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Ca),140,cal.Tphs0,'filled','o');
    scatter(cal.mem_oxd(cal.ens,cal.Si),cal.mem_oxd(cal.ens,cal.Ca),100,[0.7 0.2 0.1]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.hyp,cal.Si),cal.mem_oxd(cal.hyp,cal.Ca),100,[0.7 0.2 0.1]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Ca),100,[0.6 0.1 0.4]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Ca),100,[0.6 0.1 0.4]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),100,[0.1 0.2 0.7]*0.8,'filled','s');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),100,[0.1 0.2 0.7]*1.0,'filled','s');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Ca),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    subplot(2,4,6);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Na).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Na).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.Na),150,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.Na),150,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.Na),150,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.Na),150,cal.Tphs0,'filled','o');
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.Na),100,[0.6 0.1 0.4]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.Na),100,[0.6 0.1 0.4]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),100,[0.1 0.2 0.7]*0.8,'filled','s');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),100,[0.1 0.2 0.7]*1.0,'filled','s');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.Na),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
    subplot(2,4,7);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,cal.K).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,cal.Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,cal.K).*100,25,MAG(ic).OUT.T(hasmlt),'s'); colormap('copper');
    scatter(cal.cmp_oxd(cal.dun,cal.Si),cal.cmp_oxd(cal.dun,cal.K),150,cal.Tphs1,'filled','o');
    scatter(cal.cmp_oxd(cal.gbr,cal.Si),cal.cmp_oxd(cal.gbr,cal.K),150,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.bas,cal.Si),cal.cmp_oxd(cal.bas,cal.K),150,cal.perT ,'filled','o');
    scatter(cal.cmp_oxd(cal.rhy,cal.Si),cal.cmp_oxd(cal.rhy,cal.K),150,cal.Tphs0,'filled','o');
    scatter(cal.mem_oxd(cal.aug,cal.Si),cal.mem_oxd(cal.aug,cal.K),100,[0.6 0.1 0.4]*0.8,'filled','d');
    scatter(cal.mem_oxd(cal.pig,cal.Si),cal.mem_oxd(cal.pig,cal.K),100,[0.6 0.1 0.4]*1.4,'filled','d');
    scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.K),100,[0.1 0.2 0.7]*0.8,'filled','s');
    scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.K),100,[0.1 0.2 0.7]*1.0,'filled','s');
    scatter(cal.mem_oxd(cal.san,cal.Si),cal.mem_oxd(cal.san,cal.K),100,[0.1 0.2 0.7]*1.4,'filled','s');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.K),FS{:},TX{:})
    colorbar;
    end
end
sgtitle('melt, solid, mixture',FS{:},TX{:})


%% set ranges for control variables T, c, v, P
T = linspace(1400,700,400).';    % temperature range [degC]
c = linspace(0.51,0.51,400).';   % major component range [wt SiO2]
v = linspace(0.00,0.00,400).';   % volatile component range [wt H2O]
P = linspace(125,125,400).'*1e6; % pressure range [Pa]

% equilibrium phase fractions and compositions
cal_andes;  % load melt model calibration
c0 = c; res = 1; x = zeros(size(T)); f = zeros(size(T));
while res>1e-13
    ci = c;
    [x,cx,cm,f,vf,vm]  =  equilibrium(x,f,T,c,v,P,cal,TINY);
    c = c0.*(1-f);
    m = 1-f-x;
    res = norm(c-ci,'fro')./norm(c,'fro');
end

Nz = length(T); Nx = 1; Ptop = min(P); Pt = P; etareg = 1; calibrt = 1; T = T+273.15;
update;
T = T-273.15;

wm =  mu.*Ksgr_m .* (rhom-rho)*g0; % melt segregation speed
wx = chi.*Ksgr_x .* (rhox-rho)*g0; % crystal segregation speed
wf = phi.*Ksgr_f .* (rhof-rho)*g0; % fluid segregation speed


% block out values below solidus and above liquidus
ind = x<1e-9 | m<1e-9;
c_oxd  = squeeze(c_oxd); c_oxd(ind,:) = [];
cm_oxd = squeeze(cm_oxd); cm_oxd(ind,:) = [];
cx_oxd = squeeze(cx_oxd); cx_oxd(ind,:) = [];
c_cmp  = squeeze(c_cmp); c_cmp(ind,:) = [];
cm_cmp = squeeze(cm_cmp); cm_cmp(ind,:) = [];
cx_cmp = squeeze(cx_cmp); cx_cmp(ind,:) = [];
c_mem  = squeeze(c_mem); c_mem(ind,:) = [];
cm_mem = squeeze(cm_mem); cm_mem(ind,:) = [];
cx_mem = squeeze(cx_mem); cx_mem(ind,:) = [];
cx (ind) = [];
cm (ind) = [];
vf (ind) = [];
vm (ind) = [];
c   (ind) = [];
v   (ind) = [];
T   (ind) = [];
P   (ind) = [];
rhox(ind) = [];
rhom(ind) = [];
rhof(ind) = [];
rho (ind) = [];
etam(ind) = [];
eta (ind) = [];
wm  (ind) = [];
wx  (ind) = [];
wf  (ind) = [];
f  (ind) = [];
x  (ind) = [];
m  (ind) = [];
phi(ind) = [];
chi(ind) = [];
mu (ind) = [];
Ptop = min(P); Pt = P;

% plot phase diagram
figure(7); if ~holdfig; clf; end
for ic = nc
    ind = MAG(ic).OUT.PhaseFractions.liq_wt>=0.001 & MAG(ic).OUT.PhaseFractions.sol_wt>=0.001;
    plot(MAG(ic).OUT.OxideFract.liq(ind,1),MAG(ic).OUT.T(ind),'o','Color',[0.6,0.2,0.2]); axis tight; hold on; box on;
    plot(MAG(ic).OUT.OxideFract.sol(ind,1),MAG(ic).OUT.T(ind),'s','Color',[0.2,0.2,0.6]); axis tight; hold on; box on;
end


TT = [linspace(cal.Tphs0+Ptop*cal.clap,cal.perT+Ptop*cal.clap,200),linspace(cal.perT+Ptop*cal.clap,cal.Tphs1+Ptop*cal.clap,200)];
cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2,200),linspace((cal.perCx+cal.perCm)/2,cal.cphs0,200)];
[~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,0*TT,Ptop*ones(size(TT)),cal,TINY);
plot(CCx,TT,'k-','LineWidth',2); 
plot(CCm,TT,'k-','LineWidth',2);

% perTs  = cal.perT;
% Tphs0s = cal.Tphs0;
% Tphs1s = cal.Tphs1;
% vv = 0.10*ones(size(TT));
% xx = 0.50*ones(size(TT));
% ff = 0.05*ones(size(TT));
% for i = 1:5
%     TT = [linspace(Tphs0s+Ptop*cal.clap,perTs+Ptop*cal.clap,200),linspace(perTs+Ptop*cal.clap,Tphs1s+Ptop*cal.clap,200)];
%     cc = [linspace(cal.cphs1,(cal.perCx+cal.perCm)/2,200),linspace((cal.perCx+cal.perCm)/2,cal.cphs0,200)];
%     vmq_c0 = (4.7773e-7.*Ptop.^0.6 + 1e-11.*Ptop) .* exp(2565*(1./(TT+273.15)-1./(cal.perT+273.15))); % Katz et al., 2003; Moore et al., 1998
%     vmq_c1 = (3.5494e-3.*Ptop.^0.5 + 9.623e-8.*Ptop - 1.5223e-11.*Ptop.^1.5)./(TT+273.15) + 1.2436e-14.*Ptop.^1.5; % Liu et al., 2015
%     vmq0   = (1-cc).*vmq_c0 + cc.*vmq_c1;
%     perTs  = cal.perT -cal.dTH2O(2).*vmq0(round((perTs-Tphs1s)./(Tphs0s-Tphs1s)*200)).^0.75;
%     Tphs0s = cal.Tphs0-cal.dTH2O(1).*vmq0(1  ).^0.75;
%     Tphs1s = cal.Tphs1-cal.dTH2O(3).*vmq0(end).^0.75;
% end
% [~,CCx,CCm,~,~,~] = equilibrium(0*TT,0*TT,TT,cc,vmq0,Ptop*ones(size(TT)),cal,TINY);
% plot(CCx,TT,'k-','LineWidth',2); axis tight; hold on; box on;
% plot(CCm,TT,'k-','LineWidth',2);

Tplt = T - (Pt-Ptop)*cal.clap;
cplt = c./(1-f);
plot(cplt,Tplt,'k.','LineWidth',2,'MarkerSize',15);
plot(cx  ,Tplt,'b.','LineWidth',2,'MarkerSize',15);
plot(cm  ,Tplt,'r.','LineWidth',2,'MarkerSize',15);

% axis([0.4,0.8,700,1890]);

title('Phase Diagram','Interpreter','latex','FontSize',18)
xlabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

%% plot phase fractions
figure(2); if ~holdfig; clf; end
plot(T,x.*100,'k',T,m.*100,'r',T,f.*1000,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid $\times10$','Interpreter','latex','FontSize',15,'box','off','location','east')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Melting model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Phase fractions [wt\%]','Interpreter','latex','FontSize',15)

% plot major phase compositions
figure(3); if ~holdfig; clf; end
plot(T,cx.*100,'b',T,cm.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Major component [wt\% SiO$_2$]','Interpreter','latex','FontSize',15)

% plot volatile phase compositions
figure(4); if ~holdfig; clf; end
plot(T,vf/10.*100,'b',T,vm.*100,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('fluid $/10$','melt','Interpreter','latex','FontSize',15,'box','off','location','northeast')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase compositions','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Volatile component [wt\% H$_2$O]','Interpreter','latex','FontSize',15)

% plot phase densities
figure(5); if ~holdfig; clf; end
plot(T,rhox,'k',T,rhom,'r',T,rhof,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
plot(T,rho ,'Color',[0.5 0.5 0.5],'LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','mixture','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Density model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Density [kg/m$^3$]','Interpreter','latex','FontSize',15)

% plot mixture rheology
figure(6); if ~holdfig; clf; end
semilogy(T,eta,'k',T,etam,'r','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('mixture','melt','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Viscosity model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Viscosity [log$_{10}$ Pas]','Interpreter','latex','FontSize',15)

% plot phase segregation speeds
figure(7); if ~holdfig; clf; end
semilogy(T,max(1e-18,abs(chi.*wx)).*3600,'k',T,max(1e-18,abs(mu.*wm)).*3600,'r',T,max(1e-18,abs(phi.*wf)).*3600,'b','LineStyle',linestyle,'LineWidth',2); hold on; box on; axis tight;
legend('crystals','melt','fluid','Interpreter','latex','FontSize',15,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
title('Phase segregation model','Interpreter','latex','FontSize',18)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Segregation flux [m/hr]','Interpreter','latex','FontSize',15)

% plot oxide compositions
figure(8); if ~holdfig; clf; end
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
figure(9); if ~holdfig; clf; end
subplot(2,1,1)
sgtitle('Phase Component Fractions','Interpreter','latex','FontSize',18)
for i=1:cal.ncmp
    plot(T,cm_cmp(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); hold on; box on; axis tight;
end
legend(cal.cmpStr,'Interpreter','latex','FontSize',13,'box','off','location','best')
set(gca,'TickLabelInterpreter','latex','FontSize',13)
ylabel('Melt composition [wt\%]','Interpreter','latex','FontSize',15)
subplot(2,1,2)
for i=1:cal.ncmp
    plot(T,cx_cmp(:,i).*100,'LineStyle',linestyle,'LineWidth',2,'color',ocean(round((i-1)*213/cal.ncmp)+1,:)); hold on; box on; axis tight;
end
set(gca,'TickLabelInterpreter','latex','FontSize',13)
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Solid composition [wt\%]','Interpreter','latex','FontSize',15)

% plot simplified mineral assemblage
figure(10); clf;
patch([T;flipud(T)],[zeros(size(T));flipud(cx_mem(:,1))],0.9.*[0.6,0.8,0.5],'LineWidth',2); hold on; box on; axis tight;
patch([T;flipud(T)],[sum(cx_mem(:,1  ),2);flipud(sum(cx_mem(:,1:2),2))],1.1.*[0.6,0.8,0.5],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_mem(:,1:2),2);flipud(sum(cx_mem(:,1:4),2))],0.9.*[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_mem(:,1:4),2);flipud(sum(cx_mem(:,1:6),2))],1.1.*[0.6,0.6,0.6],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_mem(:,1:6),2);flipud(sum(cx_mem(:,1:8),2))],[0.9,0.9,0.9],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_mem(:,1:8),2);flipud(sum(cx_mem(:,1:9),2))],[0.4,0.4,0.4],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_mem(:,1:9),2);flipud(sum(cx_mem(:,1:10),2))],[1.0,0.8,0.7],'LineWidth',2);
patch([T;flipud(T)],[sum(cx_mem(:,1:10),2);flipud(sum(cx_mem(:,1:11),2))],[0.9,0.7,1.0],'LineWidth',2);
legend('~forsterite','~fayalite','~orthopyroxene','~clinopyroxene','~plagioclase','~spinel','~k-feldspar','~quartz','Interpreter','latex','FontSize',15,'box','off','location','best')
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
