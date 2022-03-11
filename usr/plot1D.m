runID    =  'fract_1D';            % run identifier
opdir    =  '../out/';             % output directory
frame    =  0:20:200;

name = [opdir,runID,'/',runID,'_par'];
load(name);

figure(1); clf;

for ff = frame
    name = [opdir,'/',runID,'/',runID,'_',num2str(ff)];
    load(name,'U','W','P','Pt','f','x','m','phi','chi','mu','H','C','V','T','c','v','cm','cx','vm','vf','IT','CT','SIm','SIx','SI','RIP','RID','it','ct','sim','six','si','rip','rid','dHdt','dCdt','dVdt','dITdt','dCTdt','dSImdt','dSIxdt','dfdt','dxdt','Gf','Gx','rho','eta','exx','ezz','exz','txx','tzz','txz','eII','tII','dt','time','step','hist');
    
    shade = max(0,min(0.95,(frame(end)-ff)./(frame(end)-frame(1))));
    
    figure(1); 
    subplot(1,5,1)
    plot(mean(T(2:end-1,2:end-1),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;
    subplot(1,5,2)
    plot(mean(c(2:end-1,2:end-1),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;    
    subplot(1,5,3)
    plot(mean(v(2:end-1,2:end-1),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;
    subplot(1,5,4)
    plot(mean(chi(2:end-1,2:end-1),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;    
    subplot(1,5,5)
    plot(mean(phi(2:end-1,2:end-1),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;
    
    figure(2);
    subplot(1,5,1)
    plot(mean(-W(:,2:end-1),2),Zfc.','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;
    subplot(1,5,2)
    plot(mean(-(chi(1:end-1,2:end-1)+chi(2:end,2:end-1))/2.*wx(:,2:end-1),2),Zfc.','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;    
    subplot(1,5,3)
    plot(mean(-(phi(1:end-1,2:end-1)+phi(2:end,2:end-1))/2.*wf(:,2:end-1),2),Zfc.','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;
    subplot(1,5,4)
    plot(mean(rho(2:end-1,2:end-1),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;    
    subplot(1,5,5)
    plot(mean(log10(eta(2:end-1,2:end-1)),2),Z(2:end-1).','LineWidth',1.5,'Color',[1,1,1].*shade); axis ij tight; hold on; box on;
end