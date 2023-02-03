clear; close all;

% load default parameters
run('../usr/par_default')

 % test decreasing grid step
NN = [50,100,200]; 

for nn = NN

    % set run parameters
    runID    =  'bnchm_mms';         % run identifier
    nop      =  1;                   % output frame plotted/saved every 'nop' time steps
    bnchm    =  1;                   % set flag for mms benchmark in fluidmech

    % set model domain parameters
    N        =  nn + 2;              % number of grid points in z-direction (incl. 2 ghosts)
    h        =  D/(N-2);             % grid spacing (equal in both dimensions, do not set) [m]

    % update inner indeces
    inz = 2:N-1;
    inx = 2:N-1;

    % create output directory
    if ~isfolder([opdir,'/',runID])
        mkdir([opdir,'/',runID]);
    end

    % run code
    run('../src/mms')

    figure(17); clf;
    colormap(ocean);
    subplot(2,3,1); imagesc(x_mms,zw_mms,-W(:,inx)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,2); imagesc(xu_mms,z_mms, U(inz,:)*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,3); imagesc(x_mms ,z_mms, P(inz,inx)/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. sol. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,4); imagesc(x_mms,zw_mms,-(W(:,inx)-W_mms(:,inx))*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $W$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,5); imagesc(xu_mms,z_mms, (U(inz,:)-U_mms(inz,:))*hr); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $U$ [m/hr]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    subplot(2,3,6); imagesc(x_mms ,z_mms, (P(inz,inx)-P_mms(inz,inx))/1e3); axis ij equal tight; box on; colorbar('TicklabelInterpreter','latex'); title('num. err. $P$ [kPa]','Interpreter','latex'); set(gca,'TicklabelInterpreter','latex')
    drawnow;

    % get solution error
    EW = norm(W(:  ,inx)-W_mms(:  ,inx),'fro')./norm(W_mms(:  ,inx),'fro');
    EU = norm(U(inz,:  )-U_mms(inz,:  ),'fro')./norm(U_mms(inz,:  ),'fro');
    EP = norm(P(inz,inx)-P_mms(inz,inx),'fro')./norm(P_mms(inz,inx),'fro');

    % plot error convergence
    fh18 = figure(18);
    p1 = loglog(h,EW,'rs','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    p2 = loglog(h,EU,'go','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    p3 = loglog(h,EP,'bv','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('grid step [m]','Interpreter','latex')
    ylabel('rel. numerical error [1]','Interpreter','latex')
    set(gca,'TicklabelInterpreter','latex')
    title('Numerical convergence in space','Interpreter','latex','FontSize',20)

    if nn == NN(1)
        p4 = loglog(D./NN,mean([EW,EU,EP]).*(NN(1)./NN).^2,'k-','LineWidth',2);  % plot linear trend for comparison
    end
    if nn == NN(end)
        legend([p1,p2,p3,p4],{'error W','error U','error P','quadratic'},'Interpreter','latex','box','on','location','southeast')
    end

    % plot error convergence
    fh19 = figure(19);
    DOFS = (NN+2).*(NN+2) + 2.*(NN+1).*(NN+2);
    dofs = (nn+2).*(nn+2) + 2.*(nn+1).*(nn+2);
    p5 = loglog(dofs,FMtime,'r+','MarkerSize',8,'LineWidth',2); axis xy tight; hold on; box on;
    set(gca,'TicklabelInterpreter','latex','FontSize',12)
    xlabel('\# dofs [1]','Interpreter','latex','FontSize',16)
    ylabel('time to solution [s]','Interpreter','latex','FontSize',16)
    title('Scaling of direct solver','Interpreter','latex','FontSize',20)

    if nn == NN(1)
        p6 = loglog(DOFS,0.95*FMtime*(DOFS./DOFS(1)).^1,'k-','LineWidth',2);  % plot linear trend for comparison
    end
    if nn == NN(end)
        legend([p5,p6],{'time to solution','linear'},'Interpreter','latex','box','on','location','southeast')
    end

end

name = [opdir,'/',runID,'/',runID,'_bnchm'];
print(fh18,name,'-dpng','-r300','-vector');

name = [opdir,'/',runID,'/',runID,'_sclng'];
print(fh19,name,'-dpng','-r300','-vector');
