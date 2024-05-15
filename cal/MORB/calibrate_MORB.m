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

filename = '/Users/tokeller/Documents/Research/nakhla/cal/MORB/MORB_frac20_nobuff_H03_ig.csv';
uiopen(filename,1)

% record Tsol, Tliq at selected P,X0 as additional constraint
Tsol = [1122.6;1118.6;1114.1;1108.8;1103.0;1096.4];   % solidus estimate from MAGEMin
Tliq = [1417.5;1416.3;1415.2;1414.0;1412.7;1411.4];   % liquidus estimate from MAGEMin
Psl  = [4.0;3.5;3.0;2.5;2.0;1.5]; % P [Pa]

%% unpack MAGEMin results
DAT = MORBfrac20nobuffH03ig;

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
PHS = zeros(npts,nphs);
for iph = 1:nphs
    PHS(hasphs(:,iph)==1,iph) = table2array(DAT(DAT.phase==phs(iph),'modewt'));
end
PHS = PHS./sum(PHS,2)*100;

% detect which phases are present on which points
hasphs = logical(PHS);

liq = 1; olv = 2; cpx = 3; fsp = 4; opx = 6; spn = 6; ilm = 7; qtz = 8;

% extract phase composition in [wt%]
OXD  = zeros(npts,nphs,noxd);
for iph = 1:nphs
    OXD(hasphs(:,iph)==1,iph,:) = table2array(DAT(DAT.phase==phs(iph),{'SiO2wt','TiO2wt','Al2O3wt','FeOwt','MgOwt','CaOwt','Na2Owt','K2Owt','H2Owt'}));
    OXD(hasphs(:,iph)==1,iph,:) = OXD(hasphs(:,iph)==1,iph,:)./sum(OXD(hasphs(:,iph)==1,iph,:),3)*100;
end

Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; K = 8; H = 9;

% lump in spinel with ilmenite
OXD(:,spn,:) = (PHS(:,spn).*OXD(:,spn,:) + PHS(:,ilm).*OXD(:,ilm,:)) ./ (PHS(:,spn) + PHS(:,ilm) + 1e-16);
OXD(:,ilm,:) = [];
PHS(:,spn)   =  PHS(:,spn) + PHS(:,ilm); 
PHS(:,ilm) = [];
phs(ilm) = [];
hasphs(:,spn) = max(hasphs(:,spn),hasphs(:,ilm));
hasphs(:,ilm) = [];
nphs = nphs-1;

% detect which oxides are present in which phases
hasoxd = logical(squeeze(sum(OXD,1)));

OXDp = OXD;

% manually remove oxides from phases as required for simplicity
% hasoxd(olv,Ca) = false;  % no Ca in olv
% hasoxd(opx,Ti) = false;  % no Ti in opx
% hasoxd(opx,Na) = false;  % no Na in opx


%% simplify mineral systems and extract end-member compositions
DATA.PRJCT  = 'MORB';
figno = 100;

kmem = 1;
MEM = [];

for iph=2:nphs

    iox = find(hasoxd(iph,:)==1);
    nox = length(iox);

    X = squeeze(OXD(hasphs(:,iph)==1,iph,hasoxd(iph,:)));
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

    OXD (hasphs(:,iph)==1,iph,hasoxd(iph,:)) = max(0,X )./sum(max(0,X ),2)*100;
    OXDp(hasphs(:,iph)==1,iph,hasoxd(iph,:)) = max(0,Xp)./sum(max(0,Xp),2)*100;

    EMExt = zeros(DGN.p,noxd);
    EMExt(:,hasoxd(iph,:))  = round(max(0,FExt)./sum(max(0,FExt),2)*100,2);
    EMExt(EMExt==max(EMExt,[],2)) = EMExt(EMExt==max(EMExt,[],2)) + 100 - sum(EMExt,2);

    [~,is] = sort(EMExt(:,Si)-EMExt(:,Mg)+EMExt(:,Na)+EMExt(:,K),'ascend');
    EMExt = EMExt(is,:);

    MEM    = [MEM;EMExt];

    figure(figno); clf; figno=figno+1;

    spz = ceil(sqrt(nox-1));
    spx = ceil((nox-1)/spz);

    kk = 2;
    for ix = 1:spx
        for iz = 1:spz
            if kk<=nox
                subplot(spz,spx,kk-1);
                scatter(squeeze(OXD (hasphs(:,iph)==1,iph,iox(1))),squeeze(OXD (hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1)); colormap('copper'); hold on
                scatter(squeeze(OXDp(hasphs(:,iph)==1,iph,iox(1))),squeeze(OXDp(hasphs(:,iph)==1,iph,iox(kk))),25,Tmp(hasphs(:,iph)==1),'filled');
                for iem = kmem:kmem+DGN.p-1
                    scatter(MEM(iem,iox(1)),MEM(iem,iox(kk)),200,'kh','filled');
                end
                xlabel(cal.oxdStr(iox(1 )),FS{:},TX{:})
                ylabel(cal.oxdStr(iox(kk)),FS{:},TX{:})
                kk = kk+1;
            else 
                break;
            end
        end
        sgtitle(phs(iph),FS{1},FS{2}+3,TX{:});
    end
    kmem = kmem+DGN.p;

end


%% add up projected mineral compositions to solid composition
SOL = zeros(npts,noxd);
wt  = zeros(size(SOL)) + 1e-16;
for iph = 2:nphs
    SOL = SOL + squeeze(OXD(:,iph,:)).*PHS(:,iph);
    wt  = wt + PHS(:,iph);
end
SOL = SOL./wt;

SOLp = zeros(npts,noxd);
wt  = zeros(size(SOLp)) + 1e-16;
for iph = 2:nphs
    SOLp = SOLp + squeeze(OXDp(:,iph,:)).*PHS(:,iph);
    wt  = wt + PHS(:,iph);
end
SOLp = SOLp./wt;


%% extract end-members encompassing all melt/solid compositions
iph = 1;
MLT  = squeeze(OXD(:,1,:));

DATA.VNAMES = cal.oxdStr(hasoxd(iph,:));
DATA.SNAMES = {};
iox = find(hasoxd(iph,1:end-1)==1);
nox = length(iox);

X = [MLT(:,1:end-1);SOLp(:,1:end-1)];
X = X./sum(X,2);

DATA.X      = X;
unmix

MLTp = MLT;
MLTp(:,1:end-1) = max(0,Xp(1:npts,:))./sum(max(0,Xp(1:npts,:)),2).*(100-OXD(:,iph,end));

OXDp(:,iph,1:end-1) = max(0,Xp(1:npts,:))./sum(max(0,Xp(1:npts,:)),2).*(100-OXD(:,iph,end));

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
SYS = zeros(npts,noxd);
wt  = zeros(size(SYS)) + 1e-16;
for iph = 1:nphs
    SYS = SYS + squeeze(OXD(:,iph,:)).*PHS(:,iph);
    wt  = wt + PHS(:,iph);
end
SYS = SYS./wt;

SYSp = zeros(npts,noxd);
wt  = zeros(size(SYSp)) + 1e-16;
for iph = 1:nphs
    SYSp = SYSp + squeeze(OXDp(:,iph,:)).*PHS(:,iph);
    wt  = wt + PHS(:,iph);
end
SYSp = SYSp./wt;


%% liquid, solid, mixture compositions
figure(figno); clf; figno=figno+1;

mrk = {'o','d','s','^','h','p','v','+'};
iox = 1:noxd;
nox = noxd;

spz = ceil(sqrt(nox-1));
spx = ceil((nox-1)/spz);

kk = 2;
for ix = 1:spx
    for iz = 1:spz
        if kk<=nox
            subplot(spz,spx,kk-1);
            scatter(squeeze(MLT (:,iox(1))),squeeze(MLT (:,iox(kk))),25,Tmp,'o'); colormap('copper'); hold on
            scatter(squeeze(MLTp(:,iox(1))),squeeze(MLTp(:,iox(kk))),25,Tmp,'o','filled');
            scatter(squeeze(SOL (:,iox(1))),squeeze(SOL (:,iox(kk))),25,Tmp,'s'); colormap('copper'); hold on
            scatter(squeeze(SOLp(:,iox(1))),squeeze(SOLp(:,iox(kk))),25,Tmp,'s','filled');
            scatter(squeeze(SYS (:,iox(1))),squeeze(SYS (:,iox(kk))),25,Tmp,'d'); colormap('copper'); hold on
            scatter(squeeze(SYSp(:,iox(1))),squeeze(SYSp(:,iox(kk))),25,Tmp,'d','filled');

            for iem = 1:size(EMInt,1)
                scatter(EMInt(iem,iox(1)),EMInt(iem,iox(kk)),200,'kh','filled');
                scatter(EMExt(iem,iox(1)),EMExt(iem,iox(kk)),200,'kh');
            end
            xlabel(cal.oxdStr(iox(1 )),FS{:},TX{:})
            ylabel(cal.oxdStr(iox(kk)),FS{:},TX{:})
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
                % for fay dps mau fau ant alb san ens fsl ulv mgt ilm qtz wat
indmem  = logical([1   1   0   0   0   0   0   0   0   0   0   0   0   0   0
                   1   1   1   0   0   1   0   0   0   0   0   0   0   0   0
                   1   1   1   1   0   1   1   0   1   0   1   0   0   0   0
                   1   1   1   1   1   1   1   1   1   1   1   1   0   0   0
                   0   1   0   1   1   0   1   1   1   1   0   1   1   0   0
                   0   0   0   0   1   0   0   1   0   1   0   0   1   1   0
                   0   0   0   0   0   0   0   0   0   0   0   0   0   0   1]);


% convert factor analysis end-member to mineral end-member proportions
cmp_oxd = 1.0*EMInt + 0.0*EMExt;
cmp_oxd = cmp_oxd./sum(cmp_oxd,2)*100;

Xp = zeros(size(cmp_oxd,1),cal.nmem);
for ip = 1:size(cmp_oxd,1)
    Xp(ip,indmem(ip,:)) = lsqnonneg(cal.mem_oxd(indmem(ip,:),1:end-1).',cmp_oxd(ip,1:end-1).');
end
cmp_mem = Xp./sum(Xp,2)*100;
cmp_mem(2:end,:) = max(2*indmem(2:end-1,:),min(75,cmp_mem(2:end,:)));
cmp_mem = cmp_mem./sum(cmp_mem,2)*100;
while max(cmp_mem(2:end,:),[],'all')>=80 || min(cmp_mem(indmem(1:end-1,:)),[],'all')<=1
    cmp_mem(2:end,:) = max(2*indmem(2:end-1,:),min(75,cmp_mem(2:end,:)));
    cmp_mem = cmp_mem./sum(cmp_mem,2)*100;
end

cmp_mem_MAP = zeros(cal.ncmp,cal.nmem);
cmp_mem_MAP(1:end-1,:) = cmp_mem;
cmp_mem_MAP(end,cal.wat) = 100;
cmp_oxd_MAP = cmp_mem_MAP*cal.mem_oxd/100;

T0_MAP = [1880  1195  1140  1060  905  795];
A_MAP  = (T0_MAP+273.15)./350;
B_MAP  = [8.7  5.0  4.2  3.4  3.0  2.5];
r_MAP  = [29.5  2.3  3.0  9.9  13.7  10.9];


%%
cal_MORB;  % load melt model calibration

PHSs = PHS;
data   = [MLTp(:);SOLp(:);PHS(:);Tsol(:);Tliq(:)];

m0     = [T0_MAP.';A_MAP.';B_MAP.';r_MAP.';cmp_mem_MAP(:).*indmem(:);];
m0_lw  = m0 - [max(10,0.05*T0_MAP.');0*max(0.1,0.01*A_MAP.');max(1,0.4*B_MAP.');max(1,0.4*r_MAP.');max(100,1*cmp_mem_MAP(:)).*indmem(:)];
m0_up  = m0 + [max(10,0.05*T0_MAP.');0*max(0.1,0.01*A_MAP.');max(1,0.4*B_MAP.');max(1,0.4*r_MAP.');max(100,1*cmp_mem_MAP(:)).*indmem(:)];
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(2,               mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/10,min(99,mbnds(4*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));
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
sigma_PHS       =  0.5 * ones(size([PHS(:)]));
sigma_TsolTliq  =  1.0 * ones(size([Tsol(:);Tliq(:)]));
sigma = [sigma_MLTSOL;sigma_PHS;sigma_TsolTliq];

% function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,Tmp,Prs,MLTp,SOLp,SYSp,PHS(:,1),Psl,cal);

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
Niter = 2e4;

% adjust step size to get reasonable acceptance ratio ~26%
anneal.initstep = 0.0003 * diff(mbnds,1,2);
anneal.levels   = 3;
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

[dhat,MLTfit,SOLfit,SYSfit,PHSfit,cmpSYS,Tsolfit,Tliqfit,Tm] = dhatFunc([T0_MAP.';A_MAP.';B_MAP.';r_MAP.';cmp_mem_MAP(:)]);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,30,50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_MAP;
cal.A      = (T0_MAP+273.15)./350;
cal.r      = r_MAP;
var.c      = repmat(mean(cmpSYS(1:5,:).*[ones(1,cal.ncmp-1),0]./sum(cmpSYS(1:5,1:end-1),2),1),length(PP),1);   % component fractions [wt]
var.P      = PP/10;         % pressure [GPa]
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
scatter(MLTp(:,Si),MLTp(:,Ti),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Ti),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,Ti),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ti),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ti),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ti),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ti),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ti),FS{:},TX{:})

subplot(2,4,2);
scatter(MLTp(:,Si),MLTp(:,Al),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Al),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,Al),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,Al),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Al),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Al),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Al),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Al),FS{:},TX{:})

subplot(2,4,3);
scatter(MLTp(:,Si),MLTp(:,Fe),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Fe),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,Fe),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,Fe),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Fe),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Fe),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Fe),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Fe),FS{:},TX{:})

subplot(2,4,4);
scatter(MLTp(:,Si),MLTp(:,Mg),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Mg),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,Mg),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,Mg),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Mg),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Mg),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Mg),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Mg),FS{:},TX{:})

subplot(2,4,5);
scatter(MLTp(:,Si),MLTp(:,Ca),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Ca),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,Ca),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,Ca),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Ca),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Ca),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Ca),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Ca),FS{:},TX{:})

subplot(2,4,6);
scatter(MLTp(:,Si),MLTp(:,Na),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,Na),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,Na),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,Na),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,Na),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,Na),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.Na),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.Na),FS{:},TX{:})

subplot(2,4,7);
scatter(MLTp(:,Si),MLTp(:,K ),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,K ),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,K ),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,K ),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,K ),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,K ),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.K ),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.K ),FS{:},TX{:})

subplot(2,4,8);
scatter(MLTp(:,Si),MLTp(:,H ),25,Tmp,'o'); colormap('copper'); axis tight; hold on
scatter(SOLp(:,Si),SOLp(:,H ),25,Tmp,'s');
scatter(SYSp(:,Si),SYSp(:,H ),25,Tmp,'d');
scatter(MLTfit(:,Si),MLTfit(:,H ),25,Tmp,'o','filled');
scatter(SOLfit(:,Si),SOLfit(:,H ),25,Tmp,'s','filled');
scatter(SYSfit(:,Si),SYSfit(:,H ),25,Tmp,'d','filled');
scatter(cmp_oxd_MAP(1:end-1,cal.Si),cmp_oxd_MAP(1:end-1,cal.H ),200,'kh','filled');
xlabel(cal.oxdStr(cal.Si),FS{:},TX{:})
ylabel(cal.oxdStr(cal.H ),FS{:},TX{:})

sgtitle('MCMC component fit',FS{:},TX{:})
drawnow

figure(109); clf; cmap = colororder;
plot(Tmp,PHS(:,1),'-','Color',cmap(2,:),'LineWidth',1.5); axis tight; hold on % melt
plot(Tmp,PHS(:,2),'-','Color',cmap(5,:),'LineWidth',1.5); % olv
plot(Tmp,PHS(:,3),'-','Color',cmap(6,:),'LineWidth',1.5); % opx
plot(Tmp,PHS(:,4),'-','Color',cmap(3,:),'LineWidth',1.5); % cpx
plot(Tmp,PHS(:,5),'-','Color',cmap(4,:),'LineWidth',1.5); % fsp
plot(Tmp,PHS(:,6),'-','Color',cmap(1,:),'LineWidth',1.5); % spn
plot(Tmp,PHS(:,7),'-','Color',cmap(7,:),'LineWidth',1.5); % qtz

plot(Tmp,PHSfit(:,1),'--','Color',cmap(2,:),'LineWidth',1.5); % mlt
plot(Tmp,PHSfit(:,2),'--','Color',cmap(5,:),'LineWidth',1.5); % olv
plot(Tmp,PHSfit(:,3),'--','Color',cmap(6,:),'LineWidth',1.5); % opx
plot(Tmp,PHSfit(:,4),'--','Color',cmap(3,:),'LineWidth',1.5); % cpx
plot(Tmp,PHSfit(:,5),'--','Color',cmap(4,:),'LineWidth',1.5); % fsp
plot(Tmp,PHSfit(:,6),'--','Color',cmap(1,:),'LineWidth',1.5); % spn
plot(Tmp,PHSfit(:,7),'--','Color',cmap(7,:),'LineWidth',1.5); % qtz
legend(['mlt',cal.msyStr],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]',FS{:},TX{:})
ylabel('Phase proportions [wt\%]',FS{:},TX{:})
sgtitle('Phase stability',FS{:},TX{:})
drawnow

% plot phase diagram
figure(110); clf;
subplot(2,4,1)
plot(MLTp(:,Si)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Si)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Si)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Si)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Si)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Si)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Si},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,2)
plot(MLTp(:,Ti)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Ti)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Ti)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Ti)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Ti)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Ti)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Ti},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,3)
plot(MLTp(:,Al)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Al)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Al)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Al)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Al)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Al)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Al},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,4)
plot(MLTp(:,Fe)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Fe)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Fe)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Fe)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Fe)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Fe)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Fe},' [wt]'],'Interpreter','latex','FontSize',15)
ylabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)

subplot(2,4,5)
plot(MLTp(:,Mg)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Mg)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Mg)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Mg)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Mg)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Mg)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Mg},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,6)
plot(MLTp(:,Ca)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Ca)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Ca)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Ca)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Ca)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Ca)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Ca},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,7)
plot(MLTp(:,Na)./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,Na)./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,Na)./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,Na)./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,Na)./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,Na)./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{Na},' [wt]'],'Interpreter','latex','FontSize',15)

subplot(2,4,8)
plot(MLTp(:,K )./sum(MLTp(:,1:end-1),2)*100,Tmp,'o','Color',[0.7 0.7 0.7]); axis tight; hold on; box on;
plot(SOLp(:,K )./sum(SOLp(:,1:end-1),2)*100,Tmp,'s','Color',[0.7 0.7 0.7]);
plot(SYSp(:,K )./sum(SYSp(:,1:end-1),2)*100,Tmp,'d','Color',[0.7 0.7 0.7]);

plot(MLTfit(:,K )./sum(MLTfit(:,1:end-1),2)*100,Tmp,'ro'); axis tight; hold on; box on;
plot(SOLfit(:,K )./sum(SOLfit(:,1:end-1),2)*100,Tmp,'bs');
plot(SYSfit(:,K )./sum(SYSfit(:,1:end-1),2)*100,Tmp,'kd');

xlabel([cal.oxdStr{K },' [wt]'],'Interpreter','latex','FontSize',15)
drawnow

figure(111); clf
plot(Tm,PP/10,'LineWidth',1); axis ij tight; hold on
plot(Tsolfit,Psl/10,'bo','LineWidth',2);
plot(Tliqfit,Psl/10,'ro','LineWidth',2);
plot(Tsol,Psl/10,'ko','LineWidth',1.5);
plot(Tliq,Psl/10,'ko','LineWidth',1.5);

legend([cal.cmpStr(1:end-1),'Tsol','Tliq'],FS{:},TX{:})
xlabel('Temperature [$^\circ$C]','Interpreter','latex','FontSize',15)
ylabel('Pressure [GPa]','Interpreter','latex','FontSize',15)
drawnow

%% save and display final state of calibration
save('MORB_calibration');

% values to enter into cal file
cmp_mem = round(cmp_mem_MAP,2)
cmp_oxd = round(cmp_oxd_MAP,2)
T0      = round(T0_MAP,0)
B       = round(B_MAP,2)
r       = round(r_MAP,1)
c0      = round(cmpSYS(1,1:end-1)./sum(cmpSYS(1,1:end-1),2),2)
