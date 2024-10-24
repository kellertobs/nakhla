%**************************************************************************
%*****  CALIBRATE MULTI-COMPONENT MELTING MODEL  **************************
%**************************************************************************

% prepare workspace
clear all; close all;

addpath(genpath('../'));
addpath('../../cal')
addpath('../../src')
addpath('../../../unmix')
addpath('../../../unmix/src')
load ocean
Fs = {'FontSize',12};
FS = {'FontSize',15};
FL = {'FontSize',18};
Ms = {'MarkerSize',8};
MS = {'MarkerSize',10};
ML = {'MarkerSize',12};
TX = {'Interpreter','latex'};
TL = {'TickLabelInterpreter','latex'};
LB = {'Location','best'};
LO = {'Location','bestoutside'};


%% *****  load calibration data  ******************************************

% load experimental data from Schmidt & Kraettli (2020), Table 3
Load_SKTable3;

PHS_frc = phs(:,[mlt,olv,opx,cpx,plg,spn,qtz]);
PHS_oxd = oxd(:,[mlt,olv,opx,cpx,plg,spn,qtz],:);

PHS_oxd(:,:,P ) = [];  % ignore P2O5
PHS_oxd(:,:,K ) = [];  % ignore K2O
PHS_oxd(:,:,Mn) = [];  % ignore MnO
PHS_oxd(:,:,Cr) = [];  % ignore Cr2O3

PHS_frc = PHS_frc./(sum(phs,2)+eps)*100;      % normalise remaining phases to 1
PHS_oxd = PHS_oxd./(sum(PHS_oxd,3)+eps)*100;  % normalise remaining oxides to 100%

phs = {'liq','olv','opx','cpx','fsp','spn','qtz'};

liq = 1; olv=2; opx=3; cpx=4; fsp=5; spn=6; qtz=7; nphs = 7;
Si=1; Ti=2; Al=3; Fe=4; Mg=5; Ca=6; Na=7; H=8; noxd = 8;

PHS_oxd(:,:,H) = 0;

nstg = size(PHS_oxd,1);

% lump in opx with cpx
PHS_oxd(:,cpx,:) = (PHS_frc(:,cpx).*PHS_oxd(:,cpx,:) + PHS_frc(:,opx).*PHS_oxd(:,opx,:)) ./ (PHS_frc(:,cpx) + PHS_frc(:,opx) + eps);
PHS_oxd(:,opx,:) = [];
PHS_frc(:,cpx)   =  PHS_frc(:,cpx) + PHS_frc(:,opx); 
PHS_frc(:,opx)   = [];
phs(opx)         = [];
nphs             = nphs-1;

liq = 1; olv=2; cpx=3; fsp=4; spn=5; qtz=6; nphs = 6;

% interpolate experimental results to finer grid to expand fitting constraints
Tmp     = interp1((1:nstg).',Tmp,(1:0.5:nstg).','linear','extrap');
Prs     = interp1((1:nstg).',Prs*10,(1:0.5:nstg).','linear','extrap');
PHS_frc = interp1((1:nstg).',PHS_frc,(1:0.5:nstg).','linear','extrap');
PHS_oxd = interp1((1:nstg).',PHS_oxd,(1:0.5:nstg).','linear','extrap');
PHS_oxd = PHS_oxd./max(1e-16,sum(PHS_oxd,3))*100;
npts    = length(Prs);

% load data for solidus and liquidus
Psl  = linspace(0,45,19)';
Tsol = [1128.3 1166.7 1198.3 1213.3 1235.0 1283.3 1326.7 1366.7 1391.7 1425.0 1453.3 1478.3 1505.0 1526.7 1548.3 1568.3 1585.0 1603.3 1618.3].';
Tliq = [1720.0 1731.7 1741.7 1753.3 1761.7 1771.7 1776.7 1786.7 1791.7 1798.3 1803.3 1806.7 1810.0 1815.0 1818.3 1821.7 1825.0 1838.3 1850.0].';

% detect which phases are present in which stages
hasphs = logical(PHS_frc);

% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

% remove minor oxides from phases (mean<0.50; max<1.0)
for iph = 1:nphs
    ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.50 & max(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),[],1)<1.00);
    hasoxd(iph,ilim) = false;
    PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./(sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)+eps)*100;
end
PHS_oxdp = PHS_oxd;

% extract melt oxide composition
MLT_oxd = squeeze(PHS_oxd(:,1,:));  % melt oxide composition

% extract solid phase oxide composition
SOL_oxd = zeros(npts,noxd);
wt  = zeros(size(SOL_oxd)) + eps;
for iph = 2:nphs
    SOL_oxd = SOL_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SOL_oxd = SOL_oxd./wt;  % solid oxide composition

% extract system oxide composition
SYS_oxd = zeros(npts,noxd);
wt  = zeros(size(SYS_oxd)) + eps;
for iph = 1:nphs
    SYS_oxd = SYS_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxd = SYS_oxd./wt;  % system oxide composition


%% *****  simplify mineral systems and extract end-member compositions  ***

% !!!  Run Section as is, follow unmix prompts on command line  !!!
cal_LUNA_upd;  % read cal.oxdStr from calibration file

% prep auxiliary parameters
DATA.PRJCT  = 'cal';
figno = 100;

% initialise lists
PHS_nmem = zeros(nphs,1);
MEM_oxd = [];

% loop through all solid phases
for iph=2:nphs

    % extract indices and number of oxides present in phase
    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);

    % load phase compositions into data array for analysis
    X = squeeze(PHS_oxd(hasphs(:,iph)==1,iph,hasoxd(iph,:)));
    X = X./sum(X,2);

    % if more than 2 oxides, 
    % use unmix tool to perform PCA, end-member extraction
    if nox>2 && size(X,1) >= size(X,2)
        DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
        DATA.SNAMES = {};
        DATA.X      = X;
        unmix
    % if 2 or less oxides use mean composition as pure-phase end-member
    else
        DGN.p = 1;
    end

    if DGN.p == 1
        FExt = mean(X);
    end

    % process external end-members for phase composition
    EMExt = zeros(DGN.p,noxd);
    EMExt(:,hasoxd(iph,:))  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
    EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

    % sort end-members from primitive to evolved
    [~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na)+EMExt(:,Ti)/2,'ascend');
    EMExt = EMExt(is,:);

    % add processed end-members to list
    MEM_oxd       = [MEM_oxd;EMExt];
    PHS_nmem(iph) = DGN.p;
end

% add water as last end-member
nmem    = sum(PHS_nmem);
MEM_oxd = [MEM_oxd;zeros(1,noxd-1),100.0];

% record final end-member count and display results
nmem    = sum(PHS_nmem)+1;
PHS_nmem(1) = nmem;
disp(MEM_oxd)

% !!!  set MEM_oxd => cal.mem_oxd in cal_LUNA_upd.m  !!!


%% *****  use end-members to project reduced solid, melt, system compositions

% !!! update calibration file name on following line, then Run Section  !!!
cal_LUNA_upd;  % read cal.mem_oxd from calibration file

% extract melt phase end-member composition and project back to 
% reduced oxide composition

PHS_mem = zeros(npts,nphs,nmem);
MLT_mem = zeros(npts,nmem);
SOL_mem = zeros(npts,nmem);
kmem = 1;
for iph = 1:nphs
    imem  = kmem:kmem+min(nmem-1,PHS_nmem(iph))-1;
    A     = cal.mem_oxd(imem,1:end-1).';
    b     = squeeze(PHS_oxd(:,iph,1:end-1));
    PHS_mem (:,iph,imem) = lsqregcmp(A,b,[0.01 0 1])*100;
    PHS_oxdp(:,iph,:   ) = squeeze(PHS_mem (:,iph,imem))*cal.mem_oxd(imem,:)/100 .* max(hasoxd);
    if iph==1
        MLT_mem  = squeeze(PHS_mem (:,1,:));
        MLT_oxdp = squeeze(PHS_oxdp(:,1,:));
    else
        SOL_mem(:,imem) = squeeze(PHS_mem (:,iph,imem));
        SOL_mem(:,imem) = SOL_mem(:,imem) .* PHS_frc(:,iph)./(100-PHS_frc(:,1)+eps);
        kmem = kmem+PHS_nmem(iph); 
    end
end
SOL_oxdp = SOL_mem*cal.mem_oxd/100;

% reconstitute projected system oxide composition
SYS_oxdp = zeros(npts,noxd);
wt  = zeros(size(SYS_oxdp)) + eps;
for iph = 1:nphs
    SYS_oxdp = SYS_oxdp + squeeze(PHS_oxdp(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxdp = SYS_oxdp./wt;  % projected system oxide composition


%% *****  reduce dimensionality by selecting number of pseudo-components **

% load phase compositions into data array for analysis
X = [MLT_oxdp(:,1:end-1);SOL_oxdp(:,1:end-1)];
X = X./sum(X,2);

% if more than 2 oxides,
% use unmix tool to perform PCA, end-member extraction
DATA.VNAMES = cal.oxdStr(1:end-1);
DATA.SNAMES = {};
DATA.X      = X;
unmix

ncmp = DGN.p+1;

% process external end-members for phase composition
EMInt = round(max(0,FInt)./sum(max(0,FInt),2)*100,2);
EMInt(EMInt==max(EMInt,[],2)) = EMInt(EMInt==max(EMInt,[],2)) + 100 - sum(EMInt,2);

% sort end-members from primitive to evolved
[~,is] = sort(-EMInt(:,cal.Mg)+EMInt(:,cal.Ti),'ascend');
EMInt  = EMInt(is,:);

% process external end-members for phase composition
EMExt = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

% sort end-members from primitive to evolved
[~,is] = sort(-EMExt(:,cal.Mg)+EMExt(:,cal.Ti),'ascend');
EMExt = EMExt(is,:);



%% *****  visualised calibrated end-member, phase compositions  ***********

% !!! update calibration file name on following line, then Run Section  !!!
cal_LUNA_upd;

% plot selected end-member and projected mineral system compositions
kmem = 1;
for iph=2:nphs-1

    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);

    figure(figno); clf; figno=figno+1;

    spz = ceil(sqrt(nox-1));
    spx = ceil((nox-1)/spz);

    kk = 2;
    for ix = 1:spx
        for iz = 1:spz
            if kk<=nox
                subplot(spz,spx,kk-1);
                scatter(squeeze(PHS_oxd (hasphs(:,iph)==1,iph,iox(1))),squeeze(PHS_oxd (hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1)); colormap('copper'); hold on
                scatter(squeeze(PHS_oxdp(hasphs(:,iph)==1,iph,iox(1))),squeeze(PHS_oxdp(hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1),'filled');
                for iem = kmem:kmem+sum(cal.msy_mem(iph-1,:))-1
                    scatter(MEM_oxd    (iem,iox(1)),MEM_oxd    (iem,iox(kk)),200,'kh');
                    scatter(cal.mem_oxd(iem,iox(1)),cal.mem_oxd(iem,iox(kk)),200,'kh','filled');
                end
                xlabel(cal.oxdStr(iox(1 )),FS{:},TX{:})
                ylabel(cal.oxdStr(iox(kk)),FS{:},TX{:})
                kk = kk+1;
            else 
                break;
            end
        end
    end
    sgtitle(phs(iph),FL{:},TX{:});
    kmem = kmem+sum(cal.msy_mem(iph-1,:));
    drawnow;
end

% plot fitted liquid, solid, mixture compositions
figure(figno); clf; figno=figno+1;

spz = ceil(sqrt(noxd-1));
spx = ceil((noxd-1)/spz);

kk = 2;
ioxd = [5 1 2 3 4 6 7 8];
for ix = 1:spx
    for iz = 1:spz
        if kk<=noxd
            subplot(spz,spx,kk-1);
            scatter(MLT_oxd (:,ioxd(1)),MLT_oxd (:,ioxd(kk)),25,Tmp,'o'); colormap('copper'); axis tight; hold on
            scatter(SOL_oxd (:,ioxd(1)),SOL_oxd (:,ioxd(kk)),25,Tmp,'s');
            scatter(SYS_oxd (:,ioxd(1)),SYS_oxd (:,ioxd(kk)),25,Tmp,'d');
            scatter(MLT_oxdp(:,ioxd(1)),MLT_oxdp(:,ioxd(kk)),25,Tmp,'o','filled');
            scatter(SOL_oxdp(:,ioxd(1)),SOL_oxdp(:,ioxd(kk)),25,Tmp,'s','filled');
            scatter(SYS_oxdp(:,ioxd(1)),SYS_oxdp(:,ioxd(kk)),25,Tmp,'d','filled');
            for icp = 1:ncmp-1
                if kk<noxd
                    scatter(EMInt(icp,5),EMInt(icp,ioxd(kk)),200,'kh','filled');
                end
            end
            if kk==noxd; legend([{'orig. mlt'},{'orig. sol'},{'orig. sys'},{'proj. sol'},{'proj. mlt'},{'proj. sys'},{'init cmp'}],Fs{:},TX{:},LO{:}); end
            xlabel(cal.oxdStr(ioxd( 1)),FS{:},TX{:})
            ylabel(cal.oxdStr(ioxd(kk)),FS{:},TX{:})
            set(gca,Fs{:},TL{:});
            kk = kk+1;
        else
            break;
        end
    end
end
sgtitle('MLT \& SOL MCMC fit',FL{:},TX{:})
drawnow

%% *****  save progress for later use  ************************************

% !!!  Run Section to save calibrated end-members and reduced compositions  !!!
close all;
save('DATA_processed');


%% *****  prepare for pseudo-component calibration  ***********************

% !!!  Run Section to load end-member calibration and prepare for pseudo-component calibration  !!!
load('DATA_processed');

cal_LUNA_upd;  % read calibration file

% !!!  Edit end-member appearances in pseudo-components  !!!
% - number and sequence of end-members must correspond to list in cal file
% - in sequence of appearance, add one new mineral
%   system to next pseudo-component
% - in sequence of appearance, add one new end-member of each mineral
%   system to next pseudo-component
% - phase out mineral systems and their end-members in accordance with
%   their fading or disappearance in PHS_frc

%                  for fay  ens dps pig  ant alb  ulv qtz wat
indmem  = logical([ 1   1    0   0   0    0   0    0   0   0
                    1   1    1   0   0    0   0    0   0   0
                    1   1    1   1   0    1   0    0   0   0
                    0   1    0   1   1    1   1    1   0   0
                    0   0    0   0   1    0   1    1   1   0
                    0   0    0   0   0    0   0    0   0   1]);

cmp_oxd = 1.0*EMInt + 0.0*EMExt;%*cal.mem_oxd(1:end-1,1:end-1)/100;
cmp_oxd = [cmp_oxd,zeros(ncmp-1,1)];
cmp_oxd = [cmp_oxd;zeros(1,noxd)];
cmp_oxd(end,end) =100;

cmp_mem = zeros(ncmp-1,cal.nmem);
for ic = 1:ncmp
    imem  = indmem(ic,:);
    A     = cal.mem_oxd(imem,:).';
    b     = squeeze(cmp_oxd(ic,:));
    cmp_mem(ic,imem) = lsqregcmp(A,b,[0.01 0 0 1e-3])*100;
end

cmp_mem_init = round(cmp_mem,1);
cmp_mem_init = cmp_mem_init./sum(cmp_mem_init,2)*100;
% indmem = logical(cmp_mem_init);
cmp_mem_best = cmp_mem_init;

cmp_oxd_init = cmp_mem_init*cal.mem_oxd/100;
cmp_oxd_best = cmp_oxd_init;

% set initial guess for melting point parameters
T0_init = [ 1800   1550   1250   1100   1000];  T0_best = T0_init;
A_init  = [ 7.40   4.80   3.40   2.90   2.60];   A_best =  A_init;
B_init  = [ 7.40   4.80   3.40   2.90   2.60];   B_best =  B_init;
r_init  = [22.00  21.00  18.00   7.00  12.00];   r_best =  r_init;
dT_init = 1400 * 1200./T0_init;  dT_best = dT_init;

% compose initial parameter guess
m0     = [T0_init.';A_init.';B_init.';r_init.';dT_init.';cmp_mem_init(:).*indmem(:);];

% set function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,Tmp,Prs,SYS_oxdp,PHS_frc,Psl,cal,[1,5,0.75,1e-3]);

% test fit function for initial guess
[~,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,Tsolfit,Tliqfit,~] = dhatFunc(m0);

% evaluate melting points as function of pressure
if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,max(Psl)*2,50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_init;
cal.A      = A_init;
cal.B      = B_init;
cal.r      = r_init;
cal.dTH2O  = dT_init;
var.c      = repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(PP),1);   % component fractions [wt]
var.P      = PP/10;
var.T      = 1000+PP*5e-8;
var.H2O    = zeros(size(PP));
cal.H2Osat = var.H2O+0.01;
[~,cal]    = meltmodel(var,cal,'T');
Tm         = cal.Tm;

% rescale mineral phase fractions
PHS_frc(:,2:end) = PHS_frc(:,2:end)./(100-PHS_frc(:,1)+eps)*100;

% plot basic information for initial fit
level = 1;
run('../MCMC/PlotFit.m')


%% *****  calibrate pseudo-components and melting point parameters  *******

cal_LUNA_upd;  % read calibration file


% uncomment following lines to run MCMC again with previous best fit as initial guess
T0_init = T0_best; A_init = A_best; B_init = B_best; r_init = r_best; cmp_mem_init = cmp_mem_best;
m0      = [T0_init.';A_init.';B_init.';r_init.';dT_init.';cmp_mem_init(:).*indmem(:)];

% !!!  set MCMC parameters then Run Section to execute MCMC routine  !!!
Niter           = 1e6;              % number of samples to take
anneal.initstep = 1e-4;             % adjust step size to get reasonable acceptance ratio 20-30%
anneal.levels   = 1;                % select number of annealing levels
anneal.burnin   = max(1,Niter/5 );  % set length of initial burn-in sequence
anneal.refine   = max(1,Niter/10);  % set length of final refinement sequence

% !!!  set data uncertainties to weight likelihood function  !!!
MLT_scl   = max(0.01,(MLT_oxdp(:)-min(MLT_oxdp(:)))./(max(MLT_oxdp(:))-min(MLT_oxdp(:))));
% SOL_scl   = max(0.01,(SOL_oxdp(:)-min(SOL_oxdp(:)))./(max(SOL_oxdp(:))-min(SOL_oxdp(:))));
MEM_scl   = max(0.01,(SOL_mem (:)-min(SOL_mem (:)))./(max(SOL_mem (:))-min(SOL_mem (:))));
PHS_scl   = max(0.01,(PHS_frc (:)-min(PHS_frc (:)))./(max(PHS_frc (:))-min(PHS_frc (:))));
sigma_MLT =  0.1  * MLT_scl.^0.25;       % uncertainty of melt oxide composition
% sigma_SOL =  1e6  * SOL_scl.^0.25;       % uncertainty of melt oxide composition
sigma_MEM =  0.1  * MEM_scl.^0.25;       % uncertainty of solid end-member composition
sigma_PHS =  0.1  * PHS_scl.^0.25;       % uncertainty of phase fractions
sigma_TSL =  0.3  * ones(size([Tsol(:);Tliq(:)])); % uncertainty of solidus/liquidus Temp
sigma = [sigma_MLT;sigma_MEM;sigma_PHS;sigma_TSL]; % combine all as in data vector

% load calibration data constraints into data vector
data   = [MLT_oxd(:);SOL_mem(:);PHS_frc(:);Tsol(:);Tliq(:)];

% construct lower and upper parameter bounds
dm    =[1*max(10.0,0.05*T0_init.'); ...
        1*max(0.25,0.25* A_init.'); ...
        1*max(0.25,0.25* B_init.'); ...
        1*max(0.50,0.25* r_init.'); ...
        0*max(10.0,0.25*dT_init.'); ...
        1*max(5.00,min(20.0,cmp_mem_init(:))).*indmem(:)];
m0_lw  = m0 - dm;
m0_up  = m0 + dm;
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(1*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(0.5,                  mbnds(1*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(2*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(1.0,                  mbnds(2*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(3.0,                  mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(100,                  mbnds(4*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/10,min(99.9,mbnds(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));
mbnds(m0==100 ,:) = 100;
% mbnds(m0==7.0 ,1) = 5.0;
mbnds(m0==T0_init(1),:) = T0_init(1);
mbnds(m0==T0_init(end),:) = T0_init(end);
anneal.initstep = anneal.initstep * diff(mbnds,1,2);  % resize step according to bounded bracket

% set parameter names according to info from calibration file
mNames = cell(cal.ncmp*cal.nmem+4*(cal.ncmp-1),1);
k = 1;
for j=1:5
    for i=1:cal.ncmp-1
        if j==1
            mNames{k} = ['T0:',cal.cmpStr{i}];
        elseif j==2
            mNames{k} = ['A:',cal.cmpStr{i}];
        elseif j==3
            mNames{k} = ['B:',cal.cmpStr{i}];
        elseif j==4
            mNames{k} = ['r:',cal.cmpStr{i}];
        elseif j==5
            mNames{k} = ['dT_H:',cal.cmpStr{i}];
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

% set function to apply further constraints to a proposed set of parameter values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% set function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% set function to calculate likelihood of forward model
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.1,1,max(Psl)*2,model,cal);

bestfit = m0;  % initialise bestfit from initial conditions

%*****  RUN MCMC PARAMETER FITTING ROUTINE  *******************************
tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),m0,mbnds,anneal,Niter);
RunTime(1) = toc;
%**************************************************************************

% uncomment following line to plot likelihood histograms for fitted parameters (slow!)
% plotmcmc(models, prob, [], mbnds, anneal, mNames); 

% extract best fit parameters
T0_best       = bestfit(               (1:cal.ncmp-1)).';
A_best        = bestfit(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
B_best        = bestfit(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
r_best        = bestfit(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
dT_best       = bestfit(4*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem_best  = reshape(bestfit(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd_best  = cmp_mem_best*cal.mem_oxd/100;

% evaluate forward model for best fit parameters
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmpfit,MLT_cmpfit,SYS_cmpfit,Tsolfit,Tliqfit,~] = dhatFunc(bestfit);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

% evaluate melting points as function of pressure
if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,max(Psl)*2,50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_best;
cal.A      = A_best;
cal.B      = B_best;
cal.r      = r_best;
cal.dTH2O  = dT_best;
var.c      = repmat(SYS_cmpfit(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmpfit(1,1:end-1),2),length(PP),1);   % component fractions [wt]
var.P      = PP/10;
var.T      = 1000+PP*5e-8;
var.H2O    = zeros(size(PP));
cal.H2Osat = var.H2O+0.001;
[~,cal]    = meltmodel(var,cal,'T');
Tm         = cal.Tm;

% uncomment following line to retrieve probability distributions (slow!)
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% *****  visualise best fit calibration  *********************************

%!!!  adjust desired level of detail to plot then run section  !!!
%     level = 1     only simple line plots
%     level = 2     add Harker diagrams and T-X pseudo-sections
%     level = 3     add mineral systems and display best fit parameters

level = 3;
run('../MCMC/PlotFit.m')

%% save and display calibration
save('LUNA_upd_calibration');
