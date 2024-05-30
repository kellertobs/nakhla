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
TINY = 1e-16;
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

% load MAGEMin results (window with table opens, using default selection
% click Import Selection => Import Data
filename = '/Users/tokeller/Documents/Research/nakhla/cal/MORB/MAR_mod_f25_H03_ig.csv';
uiopen(filename,1)

% Tsol, Tliq at selected P,X0 as additional constraint
Tsol = [1116.1; 1111.5; 1106.1; 1100.3; 1093.7; 1086.2; 1077.9;  964.7;  954.4;  943.6;  932.1;  919.8;  906.7;  892.7];   % solidus estimate from MAGEMin
Tliq = [1326.2; 1325.4; 1324.6; 1323.6; 1322.7; 1321.7; 1320.6; 1068.7; 1064.8; 1060.7; 1056.5; 1052.2; 1047.6; 1042.9];   % liquidus estimate from MAGEMin
Psl  = [   4.0;    3.5;    3.0;    2.5;    2.0;    1.5;    1.0;    4.0;    3.5;    3.0;    2.5;    2.0;    1.5;    1.0]; % P [Pa]


%% *****  unpack calibration data  ****************************************

% !!!  update table name on following line, then Run Section  !!!
DAT = MARmodf25H03ig;  % table name must correspond to table header above

% load phase names in order of appearance, liq first
phs = unique(string(DAT.phase),'stable');                                  % load phase list
phs(phs=='system') = [];                                                   % discard system
phs(phs=='qfm') = [];                                                      % discard fO2 buffer
phs(phs=='fl') = [];                                                       % discard fluid phase
nphs = length(phs);                                                        % record number of phases
iliq = find(strcmp(phs,'liq'));                                            % ensure liq comes first
iphs = 1:nphs; iphs(iliq) = []; iphs = [iliq,iphs];
phs  = phs(iphs);

liq = 1; olv = 2; fsp = 3; cpx = 4; spn = 5; opx = 6; ilm = 7; qtz = 8;    % set shortcut phase indices

% set oxide list in preferred sequence
oxd  = ["SiO2";"TiO2";"Al2O3";"FeO";"MgO";"CaO";"Na2O";"K2O";"H2O"];       % set major oxides
noxd = length(oxd);                                                        % record number of oxides

% extract calculation points
pts  = unique(DAT.point,'stable');                                         % point numbers
Tmp  = unique(DAT.TC,'stable');                                            % point temperatures
Prs  = unique(DAT.Pkbar,'stable');                                         % point pressures
npts = length(pts);                                                        % number of points

Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; K = 8; H = 9;      % set shortcut oxide indices

% detect which phases are stable on which points
hasphs = zeros(npts,nphs);
for iph = 1:nphs
    for ipt = 1:npts
        hasphs(ipt,iph) = any(table2array(DAT(DAT.point==ipt,'phase'))==phs(iph));
    end
end

% extract phase fractions in [wt%]
PHS_frc = zeros(npts,nphs);
for iph = 1:nphs
    PHS_frc(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'modewt'));
end
PHS_frc = PHS_frc./sum(PHS_frc,2)*100;

% extract phase densities in [kg/m3]
RHO = zeros(npts,nphs);
for iph = 1:nphs
    RHO(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'densitykgm3'));
end

% extract phase oxide compositions in [wt%]
PHS_oxd  = zeros(npts,nphs,noxd);
for iph = 1:nphs
    PHS_oxd(hasphs(:,iph)==1,iph,:) = table2array(DAT(DAT.phase==phs(iph),{'SiO2wt','TiO2wt','Al2O3wt','FeOwt','MgOwt','CaOwt','Na2Owt','K2Owt','H2Owt'}));
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)*100;
end

% lump in ilmenite with spinel
PHS_oxd(:,spn,:) = (PHS_frc(:,spn).*PHS_oxd(:,spn,:) + PHS_frc(:,ilm).*PHS_oxd(:,ilm,:)) ./ (PHS_frc(:,spn) + PHS_frc(:,ilm) + TINY);
PHS_oxd(:,ilm,:) = [];
RHO(:,spn)       = (PHS_frc(:,spn)+PHS_frc(:,ilm))./(PHS_frc(:,spn)./(RHO(:,spn)+TINY) + PHS_frc(:,ilm)./(RHO(:,ilm)+TINY) + TINY);
RHO(:,ilm)       = [];
PHS_frc(:,spn)   =  PHS_frc(:,spn) + PHS_frc(:,ilm); 
PHS_frc(:,ilm)   = [];
phs(ilm)         = [];
hasphs(:,spn)    = max(hasphs(:,spn),hasphs(:,ilm));
hasphs(:,ilm)    = [];
nphs             = nphs-1;

liq = 1; olv = 2; fsp = 3; cpx = 4; spn = 5; opx = 6; qtz = 7;             % update phase indices

% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

% remove minor oxides from phases (mean<0.18; max<1.8)
for iph = 1:nphs
    ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.18 | max(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),[],1)<1.8);
    hasoxd(iph,ilim) = false;
    PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)*100;
end
PHS_oxdp = PHS_oxd;


%% *****  simplify mineral systems and extract end-member compositions  ***

% !!!  Run Section as is, follow unmix prompts on command line  !!!
cal_MORB;  % read cal.oxdStr from calibration file

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
    if nox>2
        DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
        DATA.SNAMES = {};
        DATA.X      = X;
        unmix
    % if 2 or less oxides use mean composition as pure-phase end-member
    else
        DGN.p = 1;
        FExt = mean(X);
    end

    % process external end-members for phase composition
    EMExt = zeros(DGN.p,noxd);
    EMExt(:,hasoxd(iph,:))  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
    EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

    % sort end-members from primitive to evolved
    [~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na)+EMExt(:,K)+EMExt(:,Ti)/2,'ascend');
    EMExt = EMExt(is,:);

    % add processed end-members to list
    MEM_oxd       = [MEM_oxd;EMExt];
    PHS_nmem(iph) = DGN.p;
end

% add water as last end-member
nmem = sum(PHS_nmem);
MEM_oxd = [MEM_oxd;zeros(1,noxd-1),100.0];

% record final end-member count and display results
nmem = sum(PHS_nmem)+1;
disp(MEM_oxd)

% !!!  set MEM_oxd => cal.mem_oxd in cal_MORB.m  !!!


%% *****  use end-members to project reduced solid, melt, system compositions

% !!! update calibration file name on following line, then Run Section  !!!
cal_MORB;  % read cal.mem_oxd from calibration file

% extract solid phase oxide composition
SOL_oxd = zeros(npts,noxd);
wt  = zeros(size(SOL_oxd)) + TINY;
for iph = 2:nphs
    SOL_oxd = SOL_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SOL_oxd = SOL_oxd./wt;  % solid oxide composition

% extract solid phase end-member composition and project back to 
% reduced oxide composition
kmem   = 1;
SOL_mem = zeros(npts,cal.nmem);
for iph = 2:nphs
    for ip = 1:npts
        imem = kmem:kmem+PHS_nmem(iph)-1;
        SOL_mem(ip,imem) = lsqnonneg(cal.mem_oxd(imem,:).',squeeze(PHS_oxd(ip,iph,:)));
    end
    SOL_mem(:,imem) = SOL_mem(:,imem)./sum(SOL_mem(:,imem),2);
    SOL_mem(isnan(SOL_mem)) = 0;
    PHS_oxdp(:,iph,:) = SOL_mem(:,imem)*cal.mem_oxd(imem,:);

    SOL_mem(:,imem) = SOL_mem(:,imem) .* PHS_frc(:,iph)./(100-PHS_frc(:,1)) * 100;
    kmem = kmem+PHS_nmem(iph);
end
SOL_oxdp = SOL_mem*cal.mem_oxd/100;  % projected solid oxide composition

% extract melt oxide composition
iph = 1;
MLT_oxd  = squeeze(PHS_oxd(:,1,:));  % melt oxide composition

% extract melt end-member composition and project back to 
% reduced oxide composition
MLT_mem = zeros(npts,cal.nmem);
for ip = 1:npts
    MLT_mem(ip,:) = lsqnonneg(cal.mem_oxd.',MLT_oxd(ip,:).');
end
MLT_mem = MLT_mem./sum(MLT_mem,2) * 100;

MLT_oxdp = MLT_mem*cal.mem_oxd/100;  % projected melt oxide composition
PHS_oxdp(:,iph,:) = MLT_oxdp;

% reconstitute system oxide composition
SYS_oxd = zeros(npts,noxd);
wt  = zeros(size(SYS_oxd)) + TINY;
for iph = 1:nphs
    SYS_oxd = SYS_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxd = SYS_oxd./wt;  % system oxide composition

% reconstitute projected system oxide composition
SYS_oxdp = zeros(npts,noxd);
wt  = zeros(size(SYS_oxdp)) + TINY;
for iph = 1:nphs
    SYS_oxdp = SYS_oxdp + squeeze(PHS_oxdp(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxdp = SYS_oxdp./wt;  % projected system oxide composition


%% *****  extract end-members encompassing all melt/solid compositions  ***

% !!!  Run Section as is, follow unmix prompts on command line  !!!

% prep auxiliary parameters
iph = 1;
DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
DATA.SNAMES = {};

% extract indices and number of oxides present in phase
iox = find(hasoxd(iph,1:end-1)==1);
nox = length(iox);

% load phase compositions into data array for analysis
X = [MLT_oxdp(:,1:end-1);SOL_oxdp(:,1:end-1)];
X = X./sum(X,2);

% use unmix tool to perform PCA, end-member extraction
DATA.X = X;
unmix

% process internal end-members for melt/solid compositions
EMInt = zeros(DGN.p,noxd);
EMInt(:,1:end-1)  = round(max(0,FInt)./sum(max(0,FInt),2)*100,2);
EMInt(EMInt==max(EMInt,[],2)) = EMInt(EMInt==max(EMInt,[],2)) + 100 - sum(EMInt,2);

% sort end-members from primitive to evolved
[~,is] = sort(EMInt(:,Si)-EMInt(:,Mg)+EMInt(:,Na)+EMInt(:,K),'ascend');
EMInt = EMInt(is,:);

% process external end-members for melt/solid compositions
EMExt = zeros(DGN.p,noxd);
EMExt(:,1:end-1)  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

% sort end-members from primitive to evolved
[~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na)+EMExt(:,K),'ascend');
EMExt = EMExt(is,:);


%% *****  visualised calibrated end-member, phase compositions  ***********

% !!! update calibration file name on following line, then Run Section  !!!
cal_MORB;

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
                    scatter(MEM_oxd        (iem,iox(1)),MEM_oxd        (iem,iox(kk)),200,'kh');
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

% plot selected end-member and melt/solid/system compositions
figure(figno); clf; figno=figno+1;

iox = 1:noxd;
nox = noxd;

spz = ceil(sqrt(nox-1));
spx = ceil((nox-1)/spz);

kk = 2;
for ix = 1:spx
    for iz = 1:spz
        if kk<=nox
            subplot(spz,spx,kk-1);
            scatter(MLT_oxd (:,1),MLT_oxd (:,kk),25,Tmp,'o'); colormap('copper'); hold on
            scatter(MLT_oxdp(:,1),MLT_oxdp(:,kk),25,Tmp,'o','filled');
            scatter(SOL_oxd (:,1),SOL_oxd (:,kk),25,Tmp,'s'); colormap('copper'); hold on
            scatter(SOL_oxdp(:,1),SOL_oxdp(:,kk),25,Tmp,'s','filled');
            scatter(SYS_oxd (:,1),SYS_oxd (:,kk),25,Tmp,'d'); colormap('copper'); hold on
            scatter(SYS_oxdp(:,1),SYS_oxdp(:,kk),25,Tmp,'d','filled');
            for iem = 1:size(EMInt,1)
                scatter(EMInt(iem,1),EMInt(iem,kk),200,'kh','filled');
                scatter(EMExt(iem,1),EMExt(iem,kk),200,'kh');
            end
            xlabel(cal.oxdStr(1 ),FS{:},TX{:})
            ylabel(cal.oxdStr(kk),FS{:},TX{:})
            kk = kk+1;
        else
            break;
        end
    end
end
sgtitle('MLT \& SOL PCA',FL{:},TX{:})
drawnow


%% *****  save progress for later use  ************************************

% !!!  Run Section to save calibrated end-members and reduced compositions  !!!
close all;
save('MAGEMin_processed');


%% *****  prepare for pseudo-component calibration  ***********************

% !!!  Run Section to load end-member calibration and prepare for pseudo-component calibration  !!!
load('MAGEMin_processed');

cal_MORB;  % read calibration file

% !!!  Edit end-member appearances in pseudo-components  !!!
% - number and sequence of end-members must correspond to list in cal file
% - in sequence of appearance, add one new mineral
%   system to next pseudo-component
% - in sequence of appearance, add one new end-member of each mineral
%   system to next pseudo-component
% - phase out mineral systems and their end-members in accordance with
%   their fading or disappearance in PHS_frc

                % for fay  ant alb san  dps aug  ulv mgt ilm  hyp fsl  qtz wat
indmem  = logical([1   1    0   0   0    0   0    0   0   0    0   0    0   0    % dun
                   1   1    1   0   0    0   0    0   0   0    0   0    0   0    % tro
                   1   1    1   1   0    1   0    1   0   0    0   0    0   0    % ogb
                   0   1    1   1   1    1   1    1   1   0    1   0    0   0    % fbs
                   0   0    0   1   1    1   1    0   1   1    1   1    0   0    % tra
                   0   0    0   0   1    0   1    0   0   1    0   1    1   0    % rhy
                   0   0    0   0   0    0   0    0   0   0    0   0    0   1]); % vol

% set initial guess for pseudo-component compositions
cmp_oxd = 1.0*EMInt + 0.0*EMExt;
cmp_oxd = cmp_oxd./sum(cmp_oxd,2)*100;

Xp = zeros(cal.ncmp-1,cal.nmem);
for ic = 1:cal.ncmp-1
    Xp(ic,indmem(ic,:)) = lsqnonneg(cal.mem_oxd(indmem(ic,:),1:end-1).',cmp_oxd(ic,1:end-1).');
end
cmp_mem = Xp./sum(Xp,2)*100;
cmp_mem = max(1*indmem(1:end-1,:),min(99,cmp_mem));
cmp_mem = cmp_mem./sum(cmp_mem,2)*100;
while min(cmp_mem(indmem(1:end-1,:)),[],'all')<1
    cmp_mem = max(1*indmem(1:end-1,:),min(99,cmp_mem));
    cmp_mem = cmp_mem./sum(cmp_mem,2)*100;
end

cmp_mem_init = zeros(cal.ncmp,cal.nmem);
cmp_mem_init(1:end-1,:) = cmp_mem;
cmp_mem_init(end,cal.wat) = 100;
cmp_oxd_init = cmp_mem_init*cal.mem_oxd/100;

% set initial guess for melting point parameters
T0_init = [   1880    1265    1150    1090      980     800];
A_init  = [ 6.3000  4.2000  4.0000  3.9000   3.6000  3.1000];
B_init  = [ 8.7000  3.0000  2.6000  2.3000   1.9000  1.4000];
r_init  = [23.8000 12.0000  4.4000  9.1000  17.5000 14.0000];


%% *****  calibrate pseudo-components and melting point parameters  *******

cal_MORB;  % read calibration file

% !!!  set MCMC parameters then Run Section to execute MCMC routine  !!!
Niter           = 1e5;              % number of samples to take
anneal.initstep = 0.25e-3;           % adjust step size to get reasonable acceptance ratio 20-30%
anneal.levels   = 1;                % select number of annealing levels
anneal.burnin   = max(1,Niter/10);  % set length of initial burn-in sequence
anneal.refine   = max(1,Niter/10);  % set length of final refinement sequence

% !!!  set data uncertainties to weight likelihood function  !!!
sigma_MLT =  0.1  * ones(size([MLT_oxdp(:)]));     % uncertainty of melt oxide composition
sigma_SOL =  0.2  * ones(size([SOL_mem(:)]));      % uncertainty of solid end-member composition
sigma_PHS =  0.4  * ones(size([PHS_frc(:)]));      % uncertainty of phase fractions
sigma_TSL =  0.8  * ones(size([Tsol(:);Tliq(:)])); % uncertainty of solidus/liquidus Temp
sigma = [sigma_MLT;sigma_SOL;sigma_PHS;sigma_TSL]; % combine all as in data vector

% uncomment following line to run MCMC again with previous best fit as initial guess
% T0_init = T0_best; A_init = A_best; B_init = B_best; r_init = r_best; cmp_mem_init = cmp_mem_best;

% load calibration data constraints into data vector
data   = [MLT_oxdp(:);SOL_mem(:);PHS_frc(:);Tsol(:);Tliq(:)];

% compose initial parameter guess
m0     = [T0_init.';A_init.';B_init.';r_init.';cmp_mem_init(:).*indmem(:);];

% construct lower and upper parameter bounds
m0_lw  = m0 - [max(10,0.02*T0_init.');max(0.5,0.2*A_init.');max(0.5,0.2*B_init.');max(2,0.2*r_init.');max(100,1*cmp_mem_init(:)).*indmem(:)];
m0_up  = m0 + [max(10,0.02*T0_init.');max(0.5,0.2*A_init.');max(0.5,0.2*B_init.');max(2,0.2*r_init.');max(100,1*cmp_mem_init(:)).*indmem(:)];
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(3,                   mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/100,min(99,mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));
mbnds(m0==100 ,:) = 100;
mbnds(m0==T0_init(1),:) = T0_init(1);
anneal.initstep = anneal.initstep * diff(mbnds,1,2);  % resize step according to bounded bracket

% set parameter names according to info from calibration file
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

% set function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,Tmp,Prs,MLT_oxdp,SOL_mem,SYS_oxdp,PHS_frc,Psl,cal,[1e-4,1e-2]);

% set function to apply further constraints to a proposed set of parameter values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% set function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% set function to calculate likelihood of forward model
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0,1,model,cal);

bestfit = m0;  % initialise bestfit from initial conditions

%*****  RUN MCMC PARAMETER FITTING ROUTINE  *******************************
tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),m0,mbnds,anneal,Niter);
RunTime(1) = toc;
%**************************************************************************

% uncomment following line to plot likelihood histograms for fitted parameters (slow!)
plotmcmc(models, prob, [], mbnds, anneal, mNames); 

% extract best fit parameters
T0_best       = bestfit(               (1:cal.ncmp-1)).';
A_best        = bestfit(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
B_best        = bestfit(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
r_best        = bestfit(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem_best  = reshape(bestfit(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd_best  = cmp_mem_best*cal.mem_oxd/100;

% evaluate forward model for best fit parameters
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmp,MLT_cmp,SYS_cmp,Tsolfit,Tliqfit,~] = dhatFunc([T0_best.';A_best.';B_best.';r_best.';cmp_mem_best(:)]);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

% evaluate melting points as function of pressure
if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,2*max(Psl),50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_best;
cal.A      = A_best;
cal.B      = B_best;
cal.r      = r_best;
var.c      = repmat(SYS_cmp(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmp(1,1:end-1),2),length(PP),1);   % component fractions [wt]
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

% !!!  Run Section to visualise outcome of MCMC fitting routine !!!

% plot fitted mineral system compositions
kmem = 1;
figno = 102;
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
                scatter(squeeze(PHS_oxdp  (hasphs(:,iph)==1,iph,iox(1))),squeeze(PHS_oxdp  (hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1)); colormap('copper'); axis tight; hold on
                scatter(squeeze(PHS_oxdfit(hasphs(:,iph)==1,iph,iox(1))),squeeze(PHS_oxdfit(hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1),'filled');
                for iem = kmem:kmem+sum(cal.msy_mem(iph-1,:))-1
                    scatter(cal.mem_oxd(iem,iox(1)),cal.mem_oxd(iem,iox(kk)),200,'kh','filled');
                end
                if kk==nox; legend([{'proj.'},{'fit'},{'MEM'}],Fs{:},TX{:},LB{:}); end
                xlabel(cal.oxdStr(iox(1 )),FS{:},TX{:})
                ylabel(cal.oxdStr(iox(kk)),FS{:},TX{:})
                set(gca,Fs{:},TL{:});
                kk = kk+1;
            else 
                break;
            end
        end
    end
    sgtitle([char(phs(iph)),' MCMC fit'],FL{:},TX{:});
    kmem = kmem+sum(cal.msy_mem(iph-1,:));
    drawnow
end


% plot fitted liquid, solid, mixture compositions
figure(figno); clf; figno=figno+1;

spz = ceil(sqrt(noxd-1));
spx = ceil((noxd-1)/spz);

kk = 2;
for ix = 1:spx
    for iz = 1:spz
        if kk<=noxd
            subplot(spz,spx,kk-1);
            scatter(MLT_oxdp  (:,1),MLT_oxdp  (:,kk),25,Tmp,'o'); colormap('copper'); axis tight; hold on
            scatter(SOL_oxdp  (:,1),SOL_oxdp  (:,kk),25,Tmp,'s'); 
            scatter(SYS_oxdp  (:,1),SYS_oxdp  (:,kk),25,Tmp,'d'); 
            scatter(MLT_oxdfit(:,1),MLT_oxdfit(:,kk),25,Tmp,'o','filled');
            scatter(SOL_oxdfit(:,1),SOL_oxdfit(:,kk),25,Tmp,'s','filled');
            scatter(SYS_oxdfit(:,1),SYS_oxdfit(:,kk),25,Tmp,'d','filled');
            for iem = 1:cal.ncmp-1
                scatter(cmp_oxd_best(iem,1),cmp_oxd_best(iem,kk),200,'kh','filled');
                scatter(cmp_oxd_init(iem,1),cmp_oxd_init(iem,kk),200,'kh');
            end
            if kk==noxd; legend([{'proj. mlt'},{'proj. sol'},{'proj. sys'},{'fit sol'},{'fit mlt'},{'fit sys'},{'best cmp'},{'init cmp'}],Fs{:},TX{:},LO{:}); end
            xlabel(cal.oxdStr(1 ),FS{:},TX{:})
            ylabel(cal.oxdStr(kk),FS{:},TX{:})
            set(gca,Fs{:},TL{:});
            kk = kk+1;
        else
            break;
        end
    end
end
sgtitle('MLT \& SOL MCMC fit',FL{:},TX{:})
drawnow


% plot fitted T-X diagrams
figure(figno); clf; figno=figno+1;

spz = ceil(sqrt(noxd));
spx = ceil((noxd)/spz);

kk = 1;
for ix = 1:spx
    for iz = 1:spz
        if kk<=noxd
            subplot(spz,spx,kk);
            scatter(MLT_oxdp  (:,kk),Tmp,25,[0.7,0.7,0.7],'o'); axis tight; hold on
            scatter(SOL_oxdp  (:,kk),Tmp,25,[0.7,0.7,0.7],'s');
            scatter(SYS_oxdp  (:,kk),Tmp,25,[0.7,0.7,0.7],'d');
            scatter(MLT_oxdfit(:,kk),Tmp,25,[0.7,0.1,0.2],'o','filled');
            scatter(SOL_oxdfit(:,kk),Tmp,25,[0.2,0.1,0.7],'s','filled');
            scatter(SYS_oxdfit(:,kk),Tmp,25,[0.1,0.1,0.1],'d','filled');
            for iem = 1:size(EMInt,1)
                scatter(cmp_oxd_best(iem,kk),min(1400,T0_best(iem)),200,'kh','filled');
                scatter(cmp_oxd_init(iem,kk),min(1400,T0_init(iem)),200,'kh');
            end
            if kk==noxd; legend([{'proj. mlt'},{'proj. sol'},{'proj. sys'},{'fit mlt'},{'fit sol'},{'fit sys'},{'best cmp'},{'init cmp'}],Fs{:},TX{:},LO{:}); end
            xlabel(cal.oxdStr(kk),FS{:},TX{:})
            ylabel('Temperature [C]',FS{:},TX{:})
            set(gca,Fs{:},TL{:});
            kk = kk+1;
        else
            break;
        end
    end
end
sgtitle('MLT \& SOL MCMC fit',FL{:},TX{:})
drawnow


% plot fitted phase fractions
figure(figno); clf; figno=figno+1;
cmap = [colororder;[0 0 0]];

for iph=1:cal.nmsy+1
    plot(Tmp,PHS_frc   (:,iph),'-'  ,'Color',cmap(iph,:),'LineWidth',1.5); axis tight; hold on
end
for iph=1:cal.nmsy+1
    plot(Tmp,PHS_frcfit(:,iph),'--','Color',cmap(iph,:),'LineWidth',1.5); axis tight;
end
legend(['mlt',cal.msyStr],Fs{:},TX{:},LB{:})
xlabel('Temperature [$^\circ$C]',FS{:},TX{:})
ylabel('Phase proportions [wt\%]',FS{:},TX{:})
title('Phase stability MCMC fit',FL{:},TX{:})
set(gca,Fs{:},TL{:});
drawnow


% plot fitted melting points
figure(figno); clf; figno=figno+1;

plot(Tm,PP/10,'LineWidth',1); axis ij tight; hold on
plot(Tsol,Psl/10,'kd','LineWidth',1.5);
plot(Tliq,Psl/10,'ko','LineWidth',1.5);
plot(Tsolfit,Psl/10,'bd','LineWidth',2);
plot(Tliqfit,Psl/10,'ro','LineWidth',2);

legend([cal.cmpStr(1:end-1),'Tsol','Tliq','Tsol fit','Tliq fit'],Fs{:},TX{:},LB{:})
xlabel('Temperature [$^\circ$C]',TX{:},FS{:})
ylabel('Pressure [GPa]',TX{:},FS{:})
title('Melting points MCMC fit',FL{:},TX{:})
set(gca,Fs{:},TL{:});
drawnow


% system components fit
figure(figno); clf; figno=figno+1;

subplot(3,1,1)
plot(Tmp,MLT_cmp*100,'LineWidth',1.5); axis tight
ylabel('Melt comp. [wt\%]',TX{:},FS{:})
set(gca,Fs{:},TL{:});

subplot(3,1,2)
plot(Tmp,SOL_cmp*100,'LineWidth',1.5); axis tight
ylabel('Solid comp. [wt\%]',TX{:},FS{:})
set(gca,Fs{:},TL{:});

subplot(3,1,3)
plot(Tmp,SYS_cmp*100,'LineWidth',1.5); axis tight
legend(cal.cmpStr,Fs{:},TX{:},LB{:})
xlabel('Temperature [$^\circ$C]',TX{:},FS{:})
ylabel('System comp. [wt\%]',TX{:},FS{:})
set(gca,Fs{:},TL{:});

sgtitle('Pseudo-component evolution',FL{:},TX{:})
drawnow

% save and display calibration parameters
save('MORB_calibration');

% !!! enter the below values into cal file !!!
cmp_oxd = round(cmp_oxd_best,2)     
cmp_mem = round(cmp_mem_best,2)     % => cal.cmp_mem
T0      = round(T0_best,0)          % => cal.T0
A       = round(A_best,2)           % => cal.A
B       = round(B_best,2)           % => cal.B
r       = round(r_best,1)           % => cal.r
c0      = round(SYS_cmp(1,1:end-1)./sum(SYS_cmp(1,1:end-1),2),2) % => cal.c0
