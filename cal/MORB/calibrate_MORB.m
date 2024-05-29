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

%% load MAGEMin results

filename = '/Users/tokeller/Documents/Research/nakhla/cal/MORB/MORB_frac25_nobuff_H03_ig.csv';
uiopen(filename,1)

% record Tsol, Tliq at selected P,X0 as additional constraint
Tsol = [1119.1;1114.6;1109.6;1103.8;1097.4;1090.1;1081.9; 967.0; 956.8; 945.9; 934.4; 922.1; 909.0; 895.0];   % solidus estimate from MAGEMin
Tliq = [1417.3;1416.1;1415.0;1413.7;1412.5;1411.2;1409.8; 1067.4; 1063.4; 1059.3; 1055.0; 1050.5; 1045.9; 1041.0];   % liquidus estimate from MAGEMin
Psl  = [4.0;3.5;3.0;2.5;2.0;1.5;1.0;4.0;3.5;3.0;2.5;2.0;1.5;1.0]; % P [Pa]

%% unpack MAGEMin results
DAT = MORBfrac25nobuffH03ig;

phs = unique(string(DAT.phase),'stable');
phs(phs=='system') = [];
phs(phs=='qfm') = [];
phs(phs=='fl') = [];
nphs = length(phs);
oxd  = ["SiO2";"TiO2";"Al2O3";"FeO";"MgO";"CaO";"Na2O";"K2O";"H2O"];
noxd = length(oxd);

% extract calculation points
pts  = unique(DAT.point,'stable');
Tmp  = unique(DAT.TC,'stable');
Prs  = unique(DAT.Pkbar,'stable');
npts = length(pts);

% detect which phases are stable on which points
hasphs = zeros(npts,nphs);
for iph = 1:nphs
    for ipt = 1:npts
        hasphs(ipt,iph) = any(table2array(DAT(DAT.point==ipt,'phase'))==phs(iph));
    end
end

% extract phase mode in [wt%]
PHS_frc = zeros(npts,nphs);
for iph = 1:nphs
    PHS_frc(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'modewt'));
end
PHS_frc = PHS_frc./sum(PHS_frc,2)*100;

% get phase densities
RHO = zeros(npts,nphs);
for iph = 1:nphs
    RHO(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'densitykgm3'));
end

liq = 1; olv = 2; fsp = 3; cpx = 4; spn = 5; opx = 6; ilm = 7; qtz = 8;

% extract phase composition in [wt%]
PHS_oxd  = zeros(npts,nphs,noxd);
for iph = 1:nphs
    PHS_oxd(hasphs(:,iph)==1,iph,:) = table2array(DAT(DAT.phase==phs(iph),{'SiO2wt','TiO2wt','Al2O3wt','FeOwt','MgOwt','CaOwt','Na2Owt','K2Owt','H2Owt'}));
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)*100;
end

Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; K = 8; H = 9;

% lump in spinel with ilmenite
PHS_oxd(:,spn,:) = (PHS_frc(:,spn).*PHS_oxd(:,spn,:) + PHS_frc(:,ilm).*PHS_oxd(:,ilm,:)) ./ (PHS_frc(:,spn) + PHS_frc(:,ilm) + 1e-16);
PHS_oxd(:,ilm,:) = [];
PHS_frc(:,spn)   =  PHS_frc(:,spn) + PHS_frc(:,ilm); 
PHS_frc(:,ilm) = [];
phs(ilm) = [];
hasphs(:,spn) = max(hasphs(:,spn),hasphs(:,ilm));
hasphs(:,ilm) = [];
nphs = nphs-1;

liq = 1; olv = 2; fsp = 3; cpx = 4; spn = 5; opx = 6; qtz = 7;


% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(PHS_oxd,1)));

% remove most minor oxides from phases
for iph = 1:nphs
    ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.2);
    hasoxd(iph,ilim) = false;
    PHS_oxd(:,iph,~hasoxd(iph,:)) = 0;
    PHS_oxd(hasphs(:,iph)==1,iph,:) = PHS_oxd(hasphs(:,iph)==1,iph,:)./sum(PHS_oxd(hasphs(:,iph)==1,iph,:),3)*100;
end
PHS_oxdp = PHS_oxd;


%% simplify mineral systems and extract end-member compositions
DATA.PRJCT  = 'MORB';
figno = 100;

PHS_nmem = zeros(nphs,1);
MEM_oxd = [];

for iph=2:nphs

    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);

    X = squeeze(PHS_oxd(hasphs(:,iph)==1,iph,hasoxd(iph,:)));
    X = X./sum(X,2);

    if nox>2
        DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
        DATA.SNAMES = {};
        DATA.X      = X;
        unmix
    else
        DGN.p = 1;
        FExt = mean(X);
        Xp   = FExt.*ones(size(X));
    end

    PHS_oxd (hasphs(:,iph)==1,iph,hasoxd(iph,:)) = max(0,X )./sum(max(0,X ),2)*100;

    EMExt = zeros(DGN.p,noxd);
    EMExt(:,hasoxd(iph,:))  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
    EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

    [~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na)+EMExt(:,K)+EMExt(:,Ti)/2,'ascend');
    EMExt = EMExt(is,:);

    MEM_oxd       = [MEM_oxd;EMExt];
    PHS_nmem(iph) = DGN.p;
end
nmem = sum(PHS_nmem)+1;
disp(MEM_oxd)


%% add up projected mineral compositions to solid composition
cal_MORB;

SOL_oxd = zeros(npts,noxd);
wt  = zeros(size(SOL_oxd)) + 1e-16;
for iph = 2:nphs
    SOL_oxd = SOL_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SOL_oxd = SOL_oxd./wt;

SOL_mem = [];
kmem   = 1;
SOL_mem = zeros(npts,nmem);
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
SOL_oxdp = SOL_mem*cal.mem_oxd/100;


%% extract end-members encompassing all melt/solid compositions
cal_MORB;

iph = 1;
MLT_oxd  = squeeze(PHS_oxd(:,1,:));

MLT_mem = zeros(npts,nmem);
for ip = 1:npts
    MLT_mem(ip,:) = lsqnonneg(cal.mem_oxd.',MLT_oxd(ip,:).');
end
MLT_mem = MLT_mem./sum(MLT_mem,2) * 100;

MLT_oxdp = MLT_mem*cal.mem_oxd/100;
PHS_oxdp(:,iph,:) = MLT_oxdp;


DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
DATA.SNAMES = {};
iox = find(hasoxd(iph,1:end-1)==1);
nox = length(iox);

X = [MLT_oxdp(:,1:end-1);SOL_oxdp(:,1:end-1)];
X = X./sum(X,2);

DATA.X      = X;
unmix

EMInt = zeros(DGN.p,noxd);
EMInt(:,1:end-1)  = round(max(0,FInt)./sum(max(0,FInt),2)*100,2);
EMInt(EMInt==max(EMInt,[],2)) = EMInt(EMInt==max(EMInt,[],2)) + 100 - sum(EMInt,2);

[~,is] = sort(EMInt(:,Si)-EMInt(:,Mg)+EMInt(:,Na)+EMInt(:,K),'ascend');
EMInt = EMInt(is,:);

EMExt = zeros(DGN.p,noxd);
EMExt(:,1:end-1)  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

[~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na)+EMExt(:,K),'ascend');
EMExt = EMExt(is,:);


%% reconstruct bulk composition based on projected solid and liquid compositions
SYS_oxd = zeros(npts,noxd);
wt  = zeros(size(SYS_oxd)) + 1e-16;
for iph = 1:nphs
    SYS_oxd = SYS_oxd + squeeze(PHS_oxd(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxd = SYS_oxd./wt;

SYS_oxdp = zeros(npts,noxd);
wt  = zeros(size(SYS_oxdp)) + 1e-16;
for iph = 1:nphs
    SYS_oxdp = SYS_oxdp + squeeze(PHS_oxdp(:,iph,:)).*PHS_frc(:,iph);
    wt  = wt + PHS_frc(:,iph);
end
SYS_oxdp = SYS_oxdp./wt;


%% mineral endmember compositions
cal_MORB;

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
    sgtitle(phs(iph),FS{1},FS{2}+3,TX{:});
    kmem = kmem+sum(cal.msy_mem(iph-1,:));
    drawnow;
end


%% liquid, solid, mixture compositions
cal_MORB;

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

sgtitle('MLT \& SOL PCA',FS{:},TX{:})
drawnow


%% save progress for later use
close all;
save('MAGEMin_processed');


%% load projected data and prepare for fitting routines
load('MAGEMin_processed');
cal_MORB;  % load melt model calibration
                % for fay ant alb san dps aug ulv mgt hyp fsl qtz wat
indmem  = logical([1   0   0   0   0   0   0   0   0   0   0   0   0
                   1   1   1   0   0   0   0   0   0   0   0   0   0
                   1   1   1   1   0   1   0   0   0   0   0   0   0
                   0   1   1   1   1   1   1   1   0   1   0   0   0
                   0   0   0   1   1   1   1   1   1   1   1   0   0
                   0   0   0   0   1   0   1   0   1   0   1   1   0
                   0   0   0   0   0   0   0   0   0   0   0   0   1]);


% convert factor analysis end-member to mineral end-member proportions
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

T0_init = [1890  1411  1165  1087  983  780];
A_init  = (T0_init+273.15)./350;
B_init  = [9.0  3.2  2.7  2.3  1.9  1.4];
r_init  = [22.0  14.0  3.0  9.0  17.0  14.0];


%%
cal_MORB;  % load melt model calibration

% T0_init = T0_MAP; A_init = A_MAP; B_init = B_MAP; r_init = r_MAP; cmp_mem_init = cmp_mem_MAP;

data   = [MLT_oxdp(:);SOL_mem(:);PHS_frc(:);Tsol(:);Tliq(:)];

m0     = [T0_init.';A_init.';B_init.';r_init.';cmp_mem_init(:).*indmem(:);];
m0_lw  = m0 - [max(10,0.05*T0_init.');max(0.5,0.25*A_init.');max(0.5,0.25*B_init.');max(2,0.25*r_init.');max(100,1*cmp_mem_init(:)).*indmem(:)];
m0_up  = m0 + [max(10,0.05*T0_init.');max(0.5,0.25*A_init.');max(0.5,0.25*B_init.');max(2,0.25*r_init.');max(100,1*cmp_mem_init(:)).*indmem(:)];
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(3,                   mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/100,min(99,mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));
mbnds(m0==100 ,:) = 100;
mbnds(m0==T0_init(1),:) = T0_init(1);
% mbnds(m0==T0_init(2),:) = T0_init(2);

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
sigma_MLT =  0.1  * ones(size([MLT_oxdp(:)]));
sigma_SOL =  0.1  * ones(size([SOL_mem(:)]));
sigma_PHS =  0.1  * ones(size([PHS_frc(:)]));
% sigma_OXD =  0.5  * ones(size([PHS_oxd(:)]));
sigma_TSL =  0.5  * ones(size([Tsol(:);Tliq(:)]));
sigma = [sigma_MLT;sigma_SOL;sigma_PHS;sigma_TSL];

% function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,Tmp,Prs,MLT_oxdp,SOL_mem,SYS_oxdp,PHS_frc,Psl,cal,[0.001,0.2]);

% function to apply further constraints to a proposed set of model param values
% m --> m
ConstrFunc = @(model) ConstrFuncs('SumConstr', model, cal.ncmp, cal.nmem, 100);

% function to calculate prior probability given a set of model param values
% m --> prior prob
PriorFunc = @(model) ProbFuncs('PriorFunc', model, mbnds, 'uniform');

% function to calculate likelihood of dhat
% dhat --> likelihood 
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.1,1,model,cal);

% run MCMC algorithm
Niter = 1e4;

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.0005 * diff(mbnds,1,2);
anneal.levels   = 1;
anneal.burnin   = max(1,Niter/10);
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

[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmp,MLT_cmp,SYS_cmp,Tsolfit,Tliqfit,~] = dhatFunc([T0_MAP.';A_MAP.';B_MAP.';r_MAP.';cmp_mem_MAP(:)]);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,30,50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_MAP;
cal.A      = A_MAP;
cal.B      = B_MAP;
cal.r      = r_MAP;
var.c      = repmat(SYS_cmp(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmp(1,1:end-1),2),length(PP),1);   % component fractions [wt]
var.P      = PP/10;         % pressure [GPa]
var.T      = 1000+PP*5e-8;             % temperature [C]
var.H2O    = zeros(size(PP)); % water concentration [wt]
cal.H2Osat = var.H2O+0.001;
[~,cal]    = meltmodel(var,cal,'T');

Tm         = cal.Tm;

% retrieve distributions
% Nbins = min(500,Niter/20);
% [ppd_mcmc.m, ppd_mcmc.prob] = CalcPDF(mbnds, models(anneal.burnin:end,:), Nbins);


%% mineral endmember compositions
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
                    scatter(MEM_oxd(iem,iox(1)),MEM_oxd(iem,iox(kk)),200,'kh');
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
    sgtitle([char(phs(iph)),' MCMC fit'],FS{1},FS{2}+3,TX{:});
    kmem = kmem+sum(cal.msy_mem(iph-1,:));
    drawnow
end


% liquid, solid, mixture compositions
figure(figno); clf; figno=figno+1;

spz = ceil(sqrt(noxd-1));
spx = ceil((noxd-1)/spz);

kk = 2;
for ix = 1:spx
    for iz = 1:spz
        if kk<=noxd
            subplot(spz,spx,kk-1);
            scatter(MLT_oxdp  (:,1),MLT_oxdp  (:,kk),25,Tmp,'o'); colormap('copper'); axis tight; hold on
            scatter(MLT_oxdfit(:,1),MLT_oxdfit(:,kk),25,Tmp,'o','filled');
            scatter(SOL_oxdp  (:,1),SOL_oxdp  (:,kk),25,Tmp,'s'); colormap('copper'); hold on
            scatter(SOL_oxdfit(:,1),SOL_oxdfit(:,kk),25,Tmp,'s','filled');
            scatter(SYS_oxdp  (:,1),SYS_oxdp  (:,kk),25,Tmp,'d'); colormap('copper'); hold on
            scatter(SYS_oxdfit(:,1),SYS_oxdfit(:,kk),25,Tmp,'d','filled');

            for iem = 1:cal.ncmp-1
                scatter(cmp_oxd_MAP (iem,1),cmp_oxd_MAP (iem,kk),200,'kh','filled');
                scatter(cmp_oxd_init(iem,1),cmp_oxd_init(iem,kk),200,'kh');
            end
            xlabel(cal.oxdStr(1 ),FS{:},TX{:})
            ylabel(cal.oxdStr(kk),FS{:},TX{:})
            kk = kk+1;
        else
            break;
        end
    end
end

sgtitle('MLT \& SOL MCMC fit',FS{:},TX{:})
drawnow


% liquid, solid, mixture compositions
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
                scatter(cmp_oxd_MAP (iem,kk),T0_MAP (iem),200,'kh','filled');
                scatter(cmp_oxd_init(iem,kk),T0_init(iem),200,'kh');
            end
            xlabel(cal.oxdStr(kk),FS{:},TX{:})
            ylabel('Temperature [C]',FS{:},TX{:})
            kk = kk+1;
        else
            break;
        end
    end
end

sgtitle('MLT \& SOL MCMC fit',FS{:},TX{:})
drawnow


% phase stability fit
figure(figno); clf; figno=figno+1;
cmap = [colororder;[0 0 0]];

for iph=1:cal.nmsy+1
    plot(Tmp,PHS_frc   (:,iph),'-'  ,'Color',cmap(iph,:),'LineWidth',1.5); axis tight; hold on
end
for iph=1:cal.nmsy+1
    plot(Tmp,PHS_frcfit(:,iph),'--','Color',cmap(iph,:),'LineWidth',1.5); axis tight;
end

legend(['mlt',cal.msyStr],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]',FS{:},TX{:})
ylabel('Phase proportions [wt\%]',FS{:},TX{:})
title('Phase stability MCMC fit',FS{:},TX{:})
drawnow


% melting temperatures fit
figure(figno); clf; figno=figno+1;

plot(Tm,PP/10,'LineWidth',1); axis ij tight; hold on
plot(Tsolfit,Psl/10,'bd','LineWidth',2);
plot(Tliqfit,Psl/10,'ro','LineWidth',2);
plot(Tsol,Psl/10,'kd','LineWidth',1.5);
plot(Tliq,Psl/10,'ko','LineWidth',1.5);

legend([cal.cmpStr(1:end-1),'Tsol','Tliq'],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Pressure [GPa]','Interpreter','latex','FontSize',15)
title('Melting points MCMC fit',FS{:},TX{:})
drawnow


% system components fit
figure(figno); clf; figno=figno+1;

plot(SYS_cmp*100,'LineWidth',1.5); axis tight

legend(cal.cmpStr,FS{:},TX{:})
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Component fract. [wt\%]','Interpreter','latex','FontSize',15)
title('System components fit',FS{:},TX{:})
drawnow

% save and display final state of calibration
save('MORB_calibration');

% values to enter into cal file
cmp_mem = round(cmp_mem_MAP,2)
cmp_oxd = round(cmp_oxd_MAP,2)
T0      = round(T0_MAP,0)
B       = round(B_MAP,2)
r       = round(r_MAP,1)
c0      = round(SYS_cmp(1,1:end-1)./sum(SYS_cmp(1,1:end-1),2),2)
