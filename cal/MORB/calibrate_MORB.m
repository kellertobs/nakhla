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


%% load MAGEMin results
nc   = [1]; % number of compositions modelled
frct = [5];
MAG  = [];
%       Si Al Fe Mg Ca Na
ioxd = [ 1  2  5  4  3  7]; % oxide indices from MAGEMin to standard
Si = 1; Al = 2; FeO = 3; Mg = 4; Ca = 5; Na = 6; H = 7;
for ic = nc
    filename = ['MORB_fract',int2str(frct(ic)),'_H05_out.mat'];
    load(filename);

    % lump in free O to FeO, Cr2O3 to Al2O3, normalise to anhydrous unit sum
    phs = fieldnames(OUT.PhaseProps);
    phs = [phs(:)',{'SYS'},{'sol'}];
    for iph = 1:length(phs)
        OUT.OxideFract.(phs{iph}) = zeros(size(OUT.OxideFractions.(phs{iph})));
        OUT.OxideFractions.(phs{iph})(:, 8) = 0; % no Ti %OUT.OxideFractions.(phs{iph})(:,1) + OUT.OxideFractions.(phs{iph})(:,8);  OUT.OxideFractions.(phs{iph})(:,8)  = 0;
        OUT.OxideFractions.(phs{iph})(:, 6) = 0; % no K  %OUT.OxideFractions.(phs{iph})(:,7) + OUT.OxideFractions.(phs{iph})(:,6);  OUT.OxideFractions.(phs{iph})(:,6)  = 0;
        OUT.OxideFractions.(phs{iph})(:,10) = 0; % no Cr %OUT.OxideFractions.(phs{iph})(:,2) + OUT.OxideFractions.(phs{iph})(:,10); OUT.OxideFractions.(phs{iph})(:,10) = 0;
        OUT.OxideFractions.(phs{iph})(:, 5) = OUT.OxideFractions.(phs{iph})(:,5) + OUT.OxideFractions.(phs{iph})(:,9);  OUT.OxideFractions.(phs{iph})(:,9)  = 0; % combine FeO + O
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

    MAG(ic).OUT = OUT;
end

% calibrate mineral end-members

%% olivine system
cal_MORB;% load melt model calibration
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Fe,cal.Mg]);
DATA.SNAMES = {};
DATA.X      = [];
T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'ol')
        hasolv = MAG(ic).OUT.PhaseProps.ol(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
        DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.ol(hasolv,[Si,FeO,Mg]).*100];
        T = [T;MAG(ic).OUT.T(hasolv)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix
Xp = max(0,Xp); Xp = Xp./sum(Xp,2);


%% plot olivine system
figure(100); clf;
subplot(1,2,1);
scatter(X(:,1)*100,X(:,2)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Fe),100,[0.1 0.6 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Fe),100,[0.1 0.6 0.1]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(1,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.for,cal.Si),cal.mem_oxd(cal.for,cal.Mg),100,[0.1 0.6 0.1]*1.3,'filled');
scatter(cal.mem_oxd(cal.fay,cal.Si),cal.mem_oxd(cal.fay,cal.Mg),100,[0.1 0.6 0.1]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
colorbar;
sgtitle('olivine system',FS{:},TX{:})
drawnow

ndata = zeros(nc(end)+1,1);
for ic = nc
    hasolv = MAG(ic).OUT.PhaseProps.ol(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
    ndata(ic+1) = ndata(ic)+length(find(hasolv>0));
    MAG(ic).OUT.OxideFract.olp = zeros(size(MAG(ic).OUT.OxideFract.ol));
    MAG(ic).OUT.OxideFract.olp(hasolv,[Si,FeO,Mg]) = Xp(ndata(ic)+(1:length(find(hasolv>0))),:);
end


%% clinopyroxene system
cal_MORB; % load melt model calibration
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Al,cal.Fe,cal.Mg,cal.Ca,cal.Na]);
DATA.SNAMES = {};
DATA.X      = [];
T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'cpx')
    hascpx = MAG(ic).OUT.PhaseProps.cpx(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
    DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.cpx(hascpx,[Si,Al,FeO,Mg,Ca,Na]).*100];
    T = [T;MAG(ic).OUT.T(hascpx)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix
Xp = max(0,Xp); Xp = Xp./sum(Xp,2);

EM_cpx = max(0,Fe)*100;
EM_cpx = round(EM_cpx./sum(EM_cpx,2)*100,2)

%%
figure(102); clf;
subplot(2,3,1);
scatter(X(:,1)*100,X(:,2)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Al),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Al),100, [0.6 0.0 0.3]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,2);
scatter(X(:,1)*100,X(:,3)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Fe),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Fe),100,[0.6 0.0 0.3]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,3);
scatter(X(:,1)*100,X(:,4)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Mg),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Mg),100,[0.6 0.0 0.3]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,4);
scatter(X(:,1)*100,X(:,5)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,5)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Ca),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Ca),100,[0.6 0.0 0.3]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,5);
scatter(X(:,1)*100,X(:,6)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,6)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mau,cal.Si),cal.mem_oxd(cal.mau,cal.Na),100,[0.6 0.0 0.3]*1.3,'filled');
scatter(cal.mem_oxd(cal.fau,cal.Si),cal.mem_oxd(cal.fau,cal.Na),100,[0.6 0.0 0.3]*0.7,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
colorbar;
sgtitle('clinopyroxene system',FS{:},TX{:})
drawnow

ndata = zeros(nc(end)+1,1);
for ic = nc
    hascpx = MAG(ic).OUT.PhaseProps.cpx(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
    ndata(ic+1) = ndata(ic)+length(find(hascpx>0));
    MAG(ic).OUT.OxideFract.cpxp = zeros(size(MAG(ic).OUT.OxideFract.cpx));
    MAG(ic).OUT.OxideFract.cpxp(hascpx,[Si,Al,FeO,Mg,Ca,Na]) = Xp(ndata(ic)+(1:length(find(hascpx>0))),:);
end

%% oxides system
cal_MORB;  % load melt model calibration
DATA.PRJCT  = 'andesSVZ';
DATA.VNAMES = cal.oxdStr([cal.Al,cal.Fe,cal.Mg]);
DATA.SNAMES = {};
DATA.X      = [];
T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'spn')
        hasoxs = MAG(ic).OUT.PhaseProps.spn(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
        DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.spn(hasoxs,[Al,FeO,Mg]).*100];
        T = [T;MAG(ic).OUT.T(hasoxs)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

% unmix
X  = DATA.X;
EM_oxs = round(mean(X,1)*100,2)
Xp = repmat(EM_oxs/100,size(X,1),1);

%%
figure(103); clf;
subplot(1,2,1);
scatter(X(:,2)*100,X(:,1)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,2)*100,Xp(:,1)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Al),100,[0.4 0.4 0.4]*0.7,'filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(1,2,2);
scatter(X(:,2)*100,X(:,3)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,2)*100,Xp(:,3)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.mgt,cal.Fe),cal.mem_oxd(cal.mgt,cal.Mg),100,[0.4 0.4 0.4]*0.7,'filled');
xlabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
colorbar;
sgtitle('oxides system',FS{:},TX{:})
drawnow

ndata = zeros(nc(end)+1,1);
for ic = nc
    hasoxs = MAG(ic).OUT.PhaseProps.spn(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
    ndata(ic+1) = ndata(ic)+length(find(hasoxs>0));
    MAG(ic).OUT.OxideFract.spnp = zeros(size(MAG(ic).OUT.OxideFract.spn));
    MAG(ic).OUT.OxideFract.spnp(hasoxs,[Al,FeO,Mg]) = Xp(ndata(ic)+(1:length(find(hasoxs>0))),:);
end

%% feldspar system
cal_MORB; % load melt model calibration
DATA.PRJCT  = 'MORB';
DATA.VNAMES = cal.oxdStr([cal.Si,cal.Al,cal.Ca,cal.Na]);
DATA.SNAMES = {};
DATA.X      = [];
T      = [];
for ic = nc
    if isfield(MAG(ic).OUT.PhaseProps,'pl4T')
    hasfsp = MAG(ic).OUT.PhaseProps.pl4T(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
    DATA.X = [DATA.X;MAG(ic).OUT.OxideFract.pl4T(hasfsp,[Si,Al,Ca,Na]).*100];
    T = [T;MAG(ic).OUT.T(hasfsp)];
    end
end
DATA.X = DATA.X./sum(DATA.X,2);

unmix

%% plot feldspar system
figure(104); clf;
subplot(2,2,1);
scatter(X(:,1)*100,X(:,2)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,2)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Al),100,[0.2 0.1 0.5]*1.3,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Al),100,[0.2 0.1 0.5]*1.0,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,2,2);
scatter(X(:,1)*100,X(:,3)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,3)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Ca),100,[0.2 0.1 0.5]*1.3,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Ca),100,[0.2 0.1 0.5]*1.0,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,2,3);
scatter(X(:,1)*100,X(:,4)*100,25,T); colormap('copper'); hold on
scatter(Xp(:,1)*100,Xp(:,4)*100,25,T,'filled');
scatter(cal.mem_oxd(cal.ant,cal.Si),cal.mem_oxd(cal.ant,cal.Na),100,[0.2 0.1 0.5]*1.3,'filled');
scatter(cal.mem_oxd(cal.alb,cal.Si),cal.mem_oxd(cal.alb,cal.Na),100,[0.2 0.1 0.5]*1.0,'filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
colorbar;
sgtitle('felspar system',FS{:},TX{:})
drawnow


ndata = zeros(nc(end)+1,1);
for ic = nc
    hasfsp = MAG(ic).OUT.PhaseProps.pl4T(:,1)>=1e-4 & MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4;
    ndata(ic+1) = ndata(ic)+length(find(hasfsp>0));
    MAG(ic).OUT.OxideFract.pl4Tp = zeros(size(MAG(ic).OUT.OxideFract.pl4T));
    MAG(ic).OUT.OxideFract.pl4Tp(hasfsp,[Si,Al,Ca,Na]) = Xp(ndata(ic)+(1:length(find(hasfsp>0))),:);
end

%% add up projected mineral compositions to solid composition
for ic = nc
    wt = zeros(size(MAG(ic).OUT.T)) + 1e-16;
    MAG(ic).OUT.OxideFract.solp = zeros(size(MAG(ic).OUT.OxideFract.sol));

    MAG(ic).OUT.OxideFract.solp(hasolv,:) = MAG(ic).OUT.OxideFract.solp(hasolv,:) + MAG(ic).OUT.OxideFract.olp(hasolv,:).*MAG(ic).OUT.PhaseProps.ol(hasolv,1);
    wt(hasolv) = wt(hasolv) + MAG(ic).OUT.PhaseProps.ol(hasolv,1);

    MAG(ic).OUT.OxideFract.solp(hascpx,:) = MAG(ic).OUT.OxideFract.solp(hascpx,:) + MAG(ic).OUT.OxideFract.cpxp(hascpx,:).*MAG(ic).OUT.PhaseProps.cpx(hascpx,1);
    wt(hascpx) = wt(hascpx) + MAG(ic).OUT.PhaseProps.cpx(hascpx,1);

    MAG(ic).OUT.OxideFract.solp(hasoxs,:) = MAG(ic).OUT.OxideFract.solp(hasoxs,:) + MAG(ic).OUT.OxideFract.spnp(hasoxs,:).*MAG(ic).OUT.PhaseProps.spn(hasoxs,1);
    wt(hasoxs) = wt(hasoxs) + MAG(ic).OUT.PhaseProps.spn(hasoxs,1);

    MAG(ic).OUT.OxideFract.solp(hasfsp,:) = MAG(ic).OUT.OxideFract.solp(hasfsp,:) + MAG(ic).OUT.OxideFract.pl4Tp(hasfsp,:).*MAG(ic).OUT.PhaseProps.pl4T(hasfsp,1);
    wt(hasfsp) = wt(hasfsp) + MAG(ic).OUT.PhaseProps.pl4T(hasfsp,1);

    MAG(ic).OUT.OxideFract.solp = MAG(ic).OUT.OxideFract.solp./wt;
end

%% project melt composition into space of defined mineral end-members
% for ic = nc
%     hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4 & MAG(ic).OUT.PhaseFractions.sol_wt>=1e-4;
%     oxdliq = MAG(ic).OUT.OxideFract.liq(hasmlt,:)*100;
%     oxdliq = oxdliq./sum(oxdliq,2)*100;
% 
%     Xp = zeros(length(oxdliq),cal.nmem);
%     for ip = 1:length(oxdliq)
%         Xp(ip,:) = lsqnonneg(cal.mem_oxd.',oxdliq(ip,:).');
%     end
% 
%     MAG(ic).OUT.OxideFract.liqp = zeros(size(MAG(ic).OUT.OxideFract.liq));
%     MAG(ic).OUT.OxideFract.liqp(hasmlt,:) = Xp*cal.mem_oxd./100;
% end
hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4 & MAG(ic).OUT.PhaseFractions.sol_wt>=1e-4;
MAG(ic).OUT.OxideFract.liqp = zeros(size(MAG(ic).OUT.OxideFract.liq));
MAG(ic).OUT.OxideFract.liqp(hasmlt,:) = MAG(ic).OUT.OxideFract.liq(hasmlt,:);

%% reconstruct bulk composition based on projected solid and liquid compositions
for ic = nc
    MAG(ic).OUT.OxideFract.SYSp           = zeros(size(MAG(ic).OUT.OxideFract.SYS));
    MAG(ic).OUT.OxideFract.SYSp(hasmlt,:) = (MAG(ic).OUT.OxideFract.liqp(hasmlt,:).*MAG(ic).OUT.PhaseFractions.liq_wt(hasmlt,:)  ...
                                           + MAG(ic).OUT.OxideFract.solp(hasmlt,:).*MAG(ic).OUT.PhaseFractions.sol_wt(hasmlt,:)) ...
                                           ./(MAG(ic).OUT.PhaseFractions.liq_wt(hasmlt,:)+MAG(ic).OUT.PhaseFractions.sol_wt(hasmlt,:));

end

%% liquid, solid, mixture compositions
cal_MORB;  % load melt model calibration
figure(105); clf;
for ic = nc
    hasmlt = MAG(ic).OUT.PhaseFractions.liq_wt>=1e-4 & MAG(ic).OUT.PhaseFractions.sol_wt>=1e-4;% & MAG(ic).OUT.T>=800;
    subplot(2,3,1);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Al).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
    subplot(2,3,2);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,FeO).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,FeO).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,FeO).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,FeO).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,FeO).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,FeO).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
    subplot(2,3,3);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Mg).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
    subplot(2,3,4);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Ca).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
    subplot(2,3,5);
    scatter(MAG(ic).OUT.OxideFract.liq(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liq(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'o'); colormap('copper'); axis tight; hold on
    scatter(MAG(ic).OUT.OxideFract.liqp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.liqp(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'o','filled');
    scatter(MAG(ic).OUT.OxideFract.sol(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.sol(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'s');
    scatter(MAG(ic).OUT.OxideFract.solp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.solp(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'s','filled');
    scatter(MAG(ic).OUT.OxideFract.SYS(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYS(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'d');
    scatter(MAG(ic).OUT.OxideFract.SYSp(hasmlt,Si).*100,MAG(ic).OUT.OxideFract.SYSp(hasmlt,Na).*100,25,MAG(ic).OUT.T(hasmlt),'d','filled');
    xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
    ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
    colorbar;
end
sgtitle('melt, solid, mixture',FS{:},TX{:})
drawnow

save('MAGEMin_processed','MAG','hasmlt');

%% best fit cmp_mem by MCMC
load('MAGEMin_processed','MAG','hasmlt');
cal_MORB;  % load melt model calibration
oxdSYS = [];
oxdSOL = [];
oxdLIQ = [];
phs    = [];
T      = [];
P      = [];
for ic = nc
    oxdSYS = [oxdSYS;MAG(ic).OUT.OxideFract.SYSp(hasmlt,:)*100];
    oxdLIQ = [oxdLIQ;MAG(ic).OUT.OxideFract.liqp(hasmlt,:)*100];
    oxdSOL = [oxdSOL;MAG(ic).OUT.OxideFract.solp(hasmlt,:)*100];
    phs    = [phs;[MAG(ic).OUT.PhaseFractions.liq_wt(hasmlt)*100, ...
                   MAG(ic).OUT.PhaseProps.ol(  hasmlt,1)*100, ...
                   MAG(ic).OUT.PhaseProps.cpx( hasmlt,1)*100, ...
                   MAG(ic).OUT.PhaseProps.spn( hasmlt,1)*100, ...
                   MAG(ic).OUT.PhaseProps.pl4T(hasmlt,1)*100, ...
                   zeros(size(MAG(ic).OUT.PhaseProps.pl4T(hasmlt,1)))] ...
                   ./(MAG(ic).OUT.PhaseFractions.sol_wt(hasmlt)+MAG(ic).OUT.PhaseFractions.liq_wt(hasmlt))];
         T = [T;MAG(ic).OUT.T(hasmlt)];
         P = [P;MAG(ic).OUT.P(hasmlt)*1e8];
end

data = [oxdLIQ(:);oxdSOL(:);phs(:)];

indmem  = logical([0   0   0   0   0   1   0   0   0
                   1   1   1   1   0   1   1   0   0
                   1   1   1   1   1   1   1   0   0
                   1   1   1   1   1   1   1   0   0
                   1   1   0   1   1   0   1   1   0
                   0   0   0   0   0   0   0   0   1]);

%%
ioxd = [cal.Si,cal.Al,cal.Fe,cal.Mg,cal.Ca,cal.Na];
DATA.PRJCT  = 'ASVZ';
DATA.VNAMES = cal.oxdStr(ioxd);
DATA.SNAMES = {};
DATA.X      = [oxdLIQ(:,ioxd);oxdSOL(:,ioxd)];
DATA.X      = DATA.X./sum(DATA.X,2);

unmix;

cmp_oxd = max(0,Fi)*100;
cmp_oxd = cmp_oxd./sum(cmp_oxd,2)*100;
cmp_oxd(cmp_oxd(:,cal.Al)==max(cmp_oxd(:,cal.Al)),:) = [];
[~,iSi] = sort(cmp_oxd(:,Mg),'descend');
cmp_oxd_FINT = cmp_oxd(iMg,:);

Xp = zeros(size(cmp_oxd,1),cal.nmem);
for ip = 1:size(cmp_oxd,1)
    Xp(ip,indmem(ip+1,:)) = lsqnonneg(cal.mem_oxd(indmem(ip+1,:),ioxd).',cmp_oxd_FINT(ip,:).');
end
cmp_mem = Xp./sum(Xp,2)*100;

cmp_mem_FINT = zeros(cal.ncmp,cal.nmem);
cmp_mem_FINT(2:end-1,:) = cmp_mem;
cmp_mem_FINT(1,cal.ant) = 100;
cmp_mem_FINT(end,cal.wat) = 100;
cmp_oxd_FINT = cmp_mem_FINT*cal.mem_oxd/100;


%%

cal_MORB;  % load melt model calibration

m0     = cal.cmp_mem.*indmem;
m0_lw  = max(0,floor(m0 - max(1,0.10*m0)).*indmem);
m0_up  = max(0, ceil(m0 + max(1,0.10*m0)).*indmem);
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

sigma  = max(0.1,0.01.*data);
% sigma(length(data)/2+1:end) = sigma(length(data)/2+1:end)*2;

% function to calculate forward model
% m --> dhat
dhatFunc  = @(model) OxdFromCmpMem(model,data,oxdSYS,cal);

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
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,model,cal);

% function to calculate likelihood from model parameter values
% model --> dhat --> likelihood
LkMdFunc  = @(model) ProbFuncs('LikeFuncModel', dhatFunc, model, data, sigma);

% run MCMC algorithm
Niter = 1e5;

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.175 * diff(mbnds,1,2);
anneal.levels   = 3;
anneal.burnin   = Niter/20;
anneal.refine   = Niter/10;

% tic;
% [models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,m0,mbnds,anneal,Niter);
% RunTime(1) = toc;
% 
% % plot mcmc outputs
% xMAP = plotmcmc(models, prob, [], mbnds, anneal, mNames);

xMAP        = cal.cmp_mem(:);
cmp_mem_MAP = reshape(xMAP,cal.ncmp,cal.nmem);
cmp_oxd_MAP = cmp_mem_MAP*cal.mem_oxd/100;
dhat        = dhatFunc(xMAP);
oxdfit      = dhat(1:2*length(T)*cal.noxd);
phsfit      = reshape(dhat(2*length(T)*cal.noxd+1:end),[],cal.nmsy+1);
oxdLIQfit   = reshape(oxdfit(1:length(T)*cal.noxd,:),[],cal.noxd);
oxdSOLfit   = reshape(oxdfit(length(T)*cal.noxd+1:end,:),[],cal.noxd);
oxdSYSfit   = oxdLIQfit.*phsfit(:,1)/100 + oxdSOLfit.*(1-phsfit(:,1)/100);

OPTIONS.TolX = 1e-11;
Xp = zeros(length(oxdSYSfit),cal.ncmp);
for ip = 1:length(oxdSYSfit)
    Xp(ip,:) = lsqnonneg(cmp_oxd_MAP.',oxdSYSfit(ip,:).',OPTIONS);
end
c = Xp./sum(Xp,2);

% retrieve distributions
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% liquid, solid, mixture compositions
% cal_MORB; % load melt model calibration
 
figure(107); clf;
subplot(2,3,1);
scatter(oxdLIQ(:,Si),oxdLIQ(:,Al),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(oxdSOL(:,Si),oxdSOL(:,Al),25,T,'s');
scatter(oxdSYS(:,Si),oxdSYS(:,Al),25,T,'d');
scatter(oxdLIQfit(:,Si),oxdLIQfit(:,Al),25,T,'o','filled');
scatter(oxdSOLfit(:,Si),oxdSOLfit(:,Al),25,T,'s','filled');
scatter(oxdSYSfit(:,Si),oxdSYSfit(:,Al),25,T,'d','filled');
for i=1:cal.ncmp-1
    scatter(cmp_oxd_MAP(i,cal.Si),cmp_oxd_MAP(i,cal.Al),140,cal.T0(i),'filled','o');
end
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})
subplot(2,3,2);
scatter(oxdLIQ(:,Si),oxdLIQ(:,FeO),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(oxdSOL(:,Si),oxdSOL(:,FeO),25,T,'s');
scatter(oxdSYS(:,Si),oxdSYS(:,FeO),25,T,'d');
scatter(oxdLIQfit(:,Si),oxdLIQfit(:,FeO),25,T,'o','filled');
scatter(oxdSOLfit(:,Si),oxdSOLfit(:,FeO),25,T,'s','filled');
scatter(oxdSYSfit(:,Si),oxdSYSfit(:,FeO),25,T,'d','filled');
for i=1:cal.ncmp-1
    scatter(cmp_oxd_MAP(i,cal.Si),cmp_oxd_MAP(i,cal.Fe),140,cal.T0(i),'filled','o');
end
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})
subplot(2,3,3);
scatter(oxdLIQ(:,Si),oxdLIQ(:,Mg),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(oxdSOL(:,Si),oxdSOL(:,Mg),25,T,'s');
scatter(oxdSYS(:,Si),oxdSYS(:,Mg),25,T,'d');
scatter(oxdLIQfit(:,Si),oxdLIQfit(:,Mg),25,T,'o','filled');
scatter(oxdSOLfit(:,Si),oxdSOLfit(:,Mg),25,T,'s','filled');
scatter(oxdSYSfit(:,Si),oxdSYSfit(:,Mg),25,T,'d','filled');
for i=1:cal.ncmp-1
    scatter(cmp_oxd_MAP(i,cal.Si),cmp_oxd_MAP(i,cal.Mg),140,cal.T0(i),'filled','o');
end
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})
subplot(2,3,4);
scatter(oxdLIQ(:,Si),oxdLIQ(:,Ca),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(oxdSOL(:,Si),oxdSOL(:,Ca),25,T,'s');
scatter(oxdSYS(:,Si),oxdSYS(:,Ca),25,T,'d');
scatter(oxdLIQfit(:,Si),oxdLIQfit(:,Ca),25,T,'o','filled');
scatter(oxdSOLfit(:,Si),oxdSOLfit(:,Ca),25,T,'s','filled');
scatter(oxdSYSfit(:,Si),oxdSYSfit(:,Ca),25,T,'d','filled');
for i=1:cal.ncmp-1
    scatter(cmp_oxd_MAP(i,cal.Si),cmp_oxd_MAP(i,cal.Ca),140,cal.T0(i),'filled','o');
end
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})
subplot(2,3,5);
scatter(oxdLIQ(:,Si),oxdLIQ(:,Na),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(oxdSOL(:,Si),oxdSOL(:,Na),25,T,'s');
scatter(oxdSYS(:,Si),oxdSYS(:,Na),25,T,'d');
scatter(oxdLIQfit(:,Si),oxdLIQfit(:,Na),25,T,'o','filled');
scatter(oxdSOLfit(:,Si),oxdSOLfit(:,Na),25,T,'s','filled');
scatter(oxdSYSfit(:,Si),oxdSYSfit(:,Na),25,T,'d','filled');
for i=1:cal.ncmp-1
    scatter(cmp_oxd_MAP(i,cal.Si),cmp_oxd_MAP(i,cal.Na),140,cal.T0(i),'filled','o');
end
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})
subplot(2,3,6);
scatter(oxdLIQ(:,Si),oxdLIQ(:,H),25,T,'o'); colormap('copper'); axis tight; hold on
scatter(oxdSOL(:,Si),oxdSOL(:,H),25,T,'s');
scatter(oxdSYS(:,Si),oxdSYS(:,H),25,T,'d');
scatter(oxdLIQfit(:,Si),oxdLIQfit(:,H),25,T,'o','filled');
scatter(oxdSOLfit(:,Si),oxdSOLfit(:,H),25,T,'s','filled');
scatter(oxdSYSfit(:,Si),oxdSYSfit(:,H),25,T,'d','filled');
for i=1:cal.ncmp-1
    scatter(cmp_oxd_MAP(i,cal.Si),cmp_oxd_MAP(i,cal.H),140,cal.T0(i),'filled','o');
end
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.H),FS{:},TX{:})
colorbar;
sgtitle('melt, solid, mixture',FS{:},TX{:})

figure(108); clf; cmap = colororder;
plot(T,phs(:,1),'-',T,phsfit(:,1),'--','Color',cmap(2,:),'LineWidth',1.5); axis tight; hold on % melt
plot(T,phs(:,2),'-',T,phsfit(:,2),'--','Color',cmap(5,:),'LineWidth',1.5); % olv
plot(T,phs(:,3),'-',T,phsfit(:,3),'--','Color',cmap(3,:),'LineWidth',1.5); % cpx
plot(T,phs(:,4),'-',T,phsfit(:,4),'--','Color',cmap(4,:),'LineWidth',1.5); % spn
plot(T,phs(:,5),'-',T,phsfit(:,5),'--','Color',cmap(1,:),'LineWidth',1.5); % fsp
xlabel('Temperature [$^\circ$C]',FS{:},TX{:})
ylabel('Phase proportions [wt\%]',FS{:},TX{:})
sgtitle('Phase stability',FS{:},TX{:})
drawnow

%% best fit melting temperatures
clear cal var x m f cx cm
cal_MORB;  % load melt model calibration

OPTIONS.TolX = 1e-16;

Xp = zeros(length(oxdSYSfit),cal.ncmp);
for ip = 1:length(oxdSYSfit)
    Xp(ip,:) = lsqnonneg(cal.cmp_oxd.',oxdSYSfit(ip,:).',OPTIONS);
end
c = Xp./sum(Xp,2);

% equilibrium phase fractions and compositions

data   = [oxdLIQfit(:);oxdSOLfit(:)];

m0    = [cal.T0,cal.r].';
m0_lw  = max(5,floor(m0 - max(1,0.05*m0)));
m0_up  = max(5, ceil(m0 + max(1,0.05*m0)));
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

sigma  = max(0.1,0.01.*data);
% sigma(length(data)/2+1:end) = sigma(length(data)/2+1:end)*2;

% function to calculate forward model
% m --> dhat
dhatFunc  = @(model) OxdFromMeltTemp(model,T,P,c,cal);

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
Niter = 1e5;
Nbins = min(500,Niter/20);

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.03 * diff(mbnds,1,2);
anneal.levels   = 3;
anneal.burnin   = Niter/20;
anneal.refine   = Niter/10;

tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,m0,mbnds,anneal,Niter);
RunTime(1) = toc;

% plot mcmc outputs
xMAP = plotmcmc(models, prob, [], mbnds, anneal, mNames);

T0_MAP = xMAP(1:cal.ncmp-1);
r_MAP  = xMAP(cal.ncmp:end);
dhat   = dhatFunc(xMAP.');
cm_oxd_MAP = reshape(dhat(1:length(dhat)/2),[],cal.noxd);%*cmp_oxd_MAP;
cx_oxd_MAP = reshape(dhat(length(dhat)/2+1:end),[],cal.noxd);%*cmp_oxd_MAP;
c_oxd_MAP  = c*cmp_oxd_MAP;

% retrieve distributions
[ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% plot phase diagram
figure(6); clf;
subplot(3,3,1)
plot(oxdLIQ(:,cal.Si)./sum(oxdLIQ(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(oxdSOL(:,cal.Si)./sum(oxdSOL(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(oxdSYS(:,cal.Si)./sum(oxdSYS(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(oxdLIQfit(:,cal.Si)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.Si)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.Si)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Si)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.Si)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.Si)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Si},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,2)
plot(oxdLIQ(:,cal.Al)./sum(oxdLIQ(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(oxdSOL(:,cal.Al)./sum(oxdSOL(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(oxdSYS(:,cal.Al)./sum(oxdSYS(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(oxdLIQfit(:,cal.Al)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.Al)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.Al)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Al)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.Al)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.Al)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Al},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,3)
plot(oxdLIQ(:,cal.Fe)./sum(oxdLIQ(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(oxdSOL(:,cal.Fe)./sum(oxdSOL(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(oxdSYS(:,cal.Fe)./sum(oxdSYS(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(oxdLIQfit(:,cal.Fe)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.Fe)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.Fe)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Fe)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.Fe)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.Fe)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Fe},' [wt]'],'Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

subplot(3,3,4)
plot(oxdLIQ(:,cal.Mg)./sum(oxdLIQ(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(oxdSOL(:,cal.Mg)./sum(oxdSOL(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(oxdSYS(:,cal.Mg)./sum(oxdSYS(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(oxdLIQfit(:,cal.Mg)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.Mg)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.Mg)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Mg)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.Mg)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.Mg)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Mg},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,5)
plot(oxdLIQ(:,cal.Ca)./sum(oxdLIQ(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(oxdSOL(:,cal.Ca)./sum(oxdSOL(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(oxdSYS(:,cal.Ca)./sum(oxdSYS(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(oxdLIQfit(:,cal.Ca)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.Ca)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.Ca)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Ca)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.Ca)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.Ca)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Ca},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,6)
plot(oxdLIQ(:,cal.Na)./sum(oxdLIQ(:,1:end-1),2),T,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(oxdSOL(:,cal.Na)./sum(oxdSOL(:,1:end-1),2),T,'s','Color',[0.7 0.7 0.7]);
plot(oxdSYS(:,cal.Na)./sum(oxdSYS(:,1:end-1),2),T,'d','Color',[0.7 0.7 0.7]);

plot(oxdLIQfit(:,cal.Na)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.Na)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.Na)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.Na)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.Na)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.Na)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.Na},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(3,3,7)
plot(oxdLIQfit(:,cal.H)./sum(oxdLIQfit(:,1:end-1),2),T,'ro'); axis tight; hold on; box on;
plot(oxdSOLfit(:,cal.H)./sum(oxdSOLfit(:,1:end-1),2),T,'bs');
plot(oxdSYSfit (:,cal.H)./sum(oxdSYSfit (:,1:end-1),2),T,'kd');

plot(c_oxd_MAP (:,cal.H)./sum(c_oxd_MAP (:,1:end-1),2),T,'k.','LineWidth',2,'MarkerSize',15);
plot(cx_oxd_MAP(:,cal.H)./sum(cx_oxd_MAP(:,1:end-1),2),T,'b.','LineWidth',2,'MarkerSize',15);
plot(cm_oxd_MAP(:,cal.H)./sum(cm_oxd_MAP(:,1:end-1),2),T,'r.','LineWidth',2,'MarkerSize',15);
xlabel([cal.oxdStr{cal.H},' [wt]'],'Interpreter','latex','FontSize',15)

PlotPhaseDiagrams

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
