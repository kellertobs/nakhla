
if ~exist('level','var'); level = 0; end

if level>0

    % plot fitted melting points
    figure(102); clf;

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

    % plot fitted phase fractions
    figure(103); clf; 
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

    % plotted fitted system components
    figure(104); clf;

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

end

if level>1

    % plot fitted liquid, solid, mixture compositions
    figure(105); clf;

    spz = ceil(sqrt(noxd-1));
    spx = ceil((noxd-1)/spz);

    kk = 2;
    for ix = 1:spx
        for iz = 1:spz
            if kk<=noxd
                subplot(spz,spx,kk-1);
                scatter(MLT_oxd  (:,1),MLT_oxd  (:,kk),25,Tmp,'o'); colormap('copper'); axis tight; hold on
                scatter(SOL_oxd  (:,1),SOL_oxd  (:,kk),25,Tmp,'s');
                scatter(SYS_oxd  (:,1),SYS_oxd  (:,kk),25,Tmp,'d');
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
    figure(106); clf; 

    spz = ceil(sqrt(noxd));
    spx = ceil((noxd)/spz);

    kk = 1;
    for ix = 1:spx
        for iz = 1:spz
            if kk<=noxd
                subplot(spz,spx,kk);
                scatter(MLT_oxd  (:,kk),Tmp,25,[0.7,0.7,0.7],'o'); axis tight; hold on
                scatter(SOL_oxd  (:,kk),Tmp,25,[0.7,0.7,0.7],'s');
                scatter(SYS_oxd  (:,kk),Tmp,25,[0.7,0.7,0.7],'d');
                scatter(MLT_oxdfit(:,kk),Tmp,25,[0.7,0.1,0.2],'o','filled');
                scatter(SOL_oxdfit(:,kk),Tmp,25,[0.2,0.1,0.7],'s','filled');
                scatter(SYS_oxdfit(:,kk),Tmp,25,[0.1,0.1,0.1],'d','filled');
                for iem = 1:cal.ncmp-1
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

end

if level>2

    figno = 107;
    % plot fitted mineral system compositions
    kmem = 1;
    for iph=2:nphs-1

        iox = find(hasoxd(iph,:)==1);
        nox = length(iox);

        figure(figno); clf; figno = figno+1;

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


    % !!! enter the below values into cal file !!!
    cmp_oxd = round(cmp_oxd_best,2)
    cmp_mem = round(cmp_mem_best,2)     % => cal.cmp_mem
    T0      = round(T0_best,0)          % => cal.T0
    A       = round(A_best,2)           % => cal.A
    B       = round(B_best,2)           % => cal.B
    r       = round(r_best,1)           % => cal.r
    dTH2O   = round(dT_best,0)          % => cal.dTH2O
    c0      = round([SYS_cmp(  1,1:end-1)./sum(SYS_cmp(  1,1:end-1),2),SYS_cmp(1,end)],3) % => cal.c0
    c1      = round([SYS_cmp(end,1:end-1)./sum(SYS_cmp(end,1:end-1),2),SYS_cmp(1,end)],3) % => cal.c1

end