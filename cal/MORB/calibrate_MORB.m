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
pts  = unique(DAT.point,'stable'); offset = min(pts)-1; pts = pts-offset;  % point numbers
Tmp  = unique(DAT.TC,'stable');                                            % point temperatures
Prs  = unique(DAT.Pkbar,'stable');                                         % point pressures
npts = length(pts);                                                        % number of points

Si = 1; Ti = 2; Al = 3; Fe = 4; Mg = 5; Ca = 6; Na = 7; K = 8; H = 9;      % set shortcut oxide indices

% detect which phases are stable on which points
hasphs = zeros(npts,nphs);
for iph = 1:nphs
    for ipt = 1:npts
        hasphs(ipt,iph) = any(table2array(DAT(DAT.point==ipt+offset,'phase'))==phs(iph));
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
    ilim = find(mean(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),1)<0.5 & max(squeeze(PHS_oxd(hasphs(:,iph)==1,iph,:)),[],1)<1.5);
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

% extract melt oxide composition
iph = 1;
MLT_oxd  = squeeze(PHS_oxd(:,1,:));  % melt oxide composition

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


%% *****  visualised calibrated end-member, phase compositions  ***********

% !!! update calibration file name on following line, then Run Section  !!!
cal_MORB_6c;

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


%% *****  save progress for later use  ************************************

% !!!  Run Section to save calibrated end-members and reduced compositions  !!!
close all;
save('MAGEMin_processed');


%% *****  prepare for pseudo-component calibration  ***********************

% !!!  Run Section to load end-member calibration and prepare for pseudo-component calibration  !!!
load('MAGEMin_processed');

cal_MORB_6c;  % read calibration file

% !!!  Edit end-member appearances in pseudo-components  !!!
% - number and sequence of end-members must correspond to list in cal file
% - in sequence of appearance, add one new mineral
%   system to next pseudo-component
% - in sequence of appearance, add one new end-member of each mineral
%   system to next pseudo-component
% - phase out mineral systems and their end-members in accordance with
%   their fading or disappearance in PHS_frc

                % for fay  ant alb san  dps aug  ulv mgt ilm  hyp fsl  qtz wat
cmp_mem_init  = [100   0    0   0   0    0   0    0   0   0    0   0    0   0    % dun
                  34  22   42   0   0    2   0    0   0   0    0   0    0   0    % tro
                   5   4   32  10   0   37  10    2   0   0    0   0    0   0    % ogb
                   5  11   24   7   0   12  25    5   6   0    5   0    0   0    % fbs
                   0   2    0  65   5    2  11    0   2   1    3   9    0   0    % tra
                   0   0    0   0  42    0   4    0   0   1    0   2   51   0    % rhy
                   0   0    0   0   0    0   0    0   0   0    0   0    0 100];  % vol
indmem = logical(cmp_mem_init);

cmp_oxd_init = cmp_mem_init*cal.mem_oxd/100;

% set initial guess for melting point parameters
T0_init = [   1875    1180    1120    1040     965     795];
A_init  = [ 4.2000  2.7000  2.5000  2.3000  2.2000  1.7000];
B_init  = [ 8.0000  2.8000  2.6000  2.4000  2.2000  1.7000];
r_init  = [30.0000   6.000  4.0000  5.0000 16.0000 10.0000];
dT_init = [   1390    1610    1725    1750    1950    1975];

% compose initial parameter guess
m0     = [T0_init.';A_init.';B_init.';r_init.';dT_init.';cmp_mem_init(:).*indmem(:);];

% set function to calculate forward model
% m --> dhat
% dhatFunc  = @(model) OxdFromCmpMem(model,MLTp,SOLp,PHS(:,1),cal);
dhatFunc  = @(model) ModelFitP(model,Tmp,Prs,MLT_oxd,SOL_mem,SYS_oxdp,PHS_frc,Psl,cal,[1e-4,1e-2]);

% test fit function for initial guess
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmp,MLT_cmp,SYS_cmp,Tsolfit,Tliqfit,~] = dhatFunc(m0);

% evaluate melting points as function of pressure
if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,max(Psl)*3,50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_init;
cal.A      = A_init;
cal.B      = B_init;
cal.r      = r_init;
cal.dTH2O  = dT_init;
var.c      = repmat(SYS_cmp(1,:).*[ones(1,cal.ncmp-1),0]./sum(SYS_cmp(1,1:end-1),2),length(PP),1);   % component fractions [wt]
var.P      = PP/10;
var.T      = 1000+PP*5e-8;
var.H2O    = zeros(size(PP));
cal.H2Osat = var.H2O+0.001;
[~,cal]    = meltmodel(var,cal,'T');
Tm         = cal.Tm;

% plot basic information for initial fit
level = 2;
run('../MCMC/PlotFit.m')


%% *****  calibrate pseudo-components and melting point parameters  *******

cal_MORB_6c;  % read calibration file


% uncomment following lines to run MCMC again with previous best fit as initial guess
% T0_init = T0_best; A_init = A_best; B_init = B_best; r_init = r_best; cmp_mem_init = cmp_mem_best;
% m0     = [T0_init.';A_init.';B_init.';r_init.';dT_init.';cmp_mem_init(:).*indmem(:);];

% !!!  set MCMC parameters then Run Section to execute MCMC routine  !!!
Niter           = 1e3;              % number of samples to take
anneal.initstep = 0.1e-2;           % adjust step size to get reasonable acceptance ratio 20-30%
anneal.levels   = 1;                % select number of annealing levels
anneal.burnin   = max(1,Niter/10);  % set length of initial burn-in sequence
anneal.refine   = max(1,Niter/10);  % set length of final refinement sequence

% !!!  set data uncertainties to weight likelihood function  !!!
MLT_scl   = max(0.001,(MLT_oxd(:)-min(MLT_oxd(:)))./(max(MLT_oxd(:))-min(MLT_oxd(:))));
SOL_scl   = max(0.001,(SOL_mem(:)-min(SOL_mem(:)))./(max(SOL_mem(:))-min(SOL_mem(:))));
PHS_scl   = max(0.001,(PHS_frc(:)-min(PHS_frc(:)))./(max(PHS_frc(:))-min(PHS_frc(:))));
sigma_MLT =  0.1  * MLT_scl.^0.3;      % uncertainty of melt oxide composition
sigma_SOL =  0.2  * SOL_scl.^0.2;      % uncertainty of solid end-member composition
sigma_PHS =  0.3  * PHS_scl.^0.1;      % uncertainty of phase fractions
sigma_TSL =  0.5  * ones(size([Tsol(:);Tliq(:)])); % uncertainty of solidus/liquidus Temp
sigma = [sigma_MLT;sigma_SOL;sigma_PHS;sigma_TSL]; % combine all as in data vector

% load calibration data constraints into data vector
data   = [MLT_oxd(:);SOL_mem(:);PHS_frc(:);Tsol(:);Tliq(:)];

% construct lower and upper parameter bounds
dm     = [max(25,0.05*T0_init.'); ...
          max( 1,0.75* A_init.'); ...
          max( 1,0.50* B_init.'); ...
          max( 5,0.75* r_init.'); ...
        1*max(50,0.10*dT_init.'); ...
          max(10,1.00*cmp_mem_init(:)).*indmem(:)];
m0_lw  = m0 - dm;
m0_up  = m0 + dm;
mbnds  = [m0_lw(:),m0_up(:)]; % model parameter bounds
mbnds(1*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(1.0,                   mbnds(1*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(2*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(1.0,                   mbnds(2*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(3.0,                   mbnds(3*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(4*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:) = max(1000,                  mbnds(4*(cal.ncmp-1)+(1:cal.ncmp-1)       ,:));
mbnds(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:) = max(indmem(:)/100,min(99.9,mbnds(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),:)));
mbnds(m0==100 ,:) = 100;
mbnds(m0==T0_init(1),:) = T0_init(1);
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
LikeFunc  = @(dhat,model) ProbFuncs('LikeFuncSimplex',dhat,data,sigma,0.05,5,max(Psl)*3,model,cal);

bestfit = m0;  % initialise bestfit from initial conditions

%*****  RUN MCMC PARAMETER FITTING ROUTINE  *******************************
tic;
[models,prob,accept,bestfit] = mcmc(dhatFunc,PriorFunc,LikeFunc,ConstrFunc,5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem),m0,mbnds,anneal,Niter);
RunTime(1) = toc;
%**************************************************************************

% uncomment following line to plot likelihood histograms for fitted parameters (slow!)
plotmcmc(models, prob, [], mbnds, anneal, mNames); 

% extract best fit parameters
T0_best       = bestfit(               (1:cal.ncmp-1)).';
A_best        = bestfit(1*(cal.ncmp-1)+(1:cal.ncmp-1)).';
B_best        = bestfit(2*(cal.ncmp-1)+(1:cal.ncmp-1)).';
r_best        = bestfit(3*(cal.ncmp-1)+(1:cal.ncmp-1)).';
dT_best       = bestfit(4*(cal.ncmp-1)+(1:cal.ncmp-1)).';
cmp_mem_best  = reshape(bestfit(5*(cal.ncmp-1)+(1:cal.ncmp*cal.nmem)),cal.ncmp,cal.nmem);
cmp_oxd_best  = cmp_mem_best*cal.mem_oxd/100;

% evaluate forward model for best fit parameters
[dhat,MLT_oxdfit,SOL_oxdfit,SYS_oxdfit,SOL_memfit,PHS_oxdfit,PHS_frcfit,SOL_cmp,MLT_cmp,SYS_cmp,Tsolfit,Tliqfit,~] = dhatFunc(bestfit);
[Lbest,Vsimplex] = LikeFunc(dhat,bestfit);

% evaluate melting points as function of pressure
if isfield(cal,'Tsol'); cal = rmfield(cal,{'Tsol' 'Tliq'}); end
PP         = linspace(0.001,max(Psl)*3,50).';
var.m      = ones(size(PP))/2; var.x = var.m; var.f = 0*var.m;
cal.T0     = T0_best;
cal.A      = A_best;
cal.B      = B_best;
cal.r      = r_best;
cal.dTH2O  = dT_best;
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

%!!!  adjust desired level of detail to plot then run section  !!!
%     level = 1     only simple line plots
%     level = 2     add Harker diagrams and T-X pseudo-sections
%     level = 3     add mineral systems and display best fit parameters

level = 3;
PlotFit.m

%% save and display calibration
save('MORB_calibration_6cmp');
