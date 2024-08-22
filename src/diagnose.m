% print diagnostics
fprintf(1,'\n         total time to solution = %3.3f sec\n\n',toc(TTtime));
fprintf(1,'         thermo-chemical solve  = %1.3e sec\n'  ,TCtime/(iter-1));
fprintf(1,'         phase equilibr. solve  = %1.3e sec\n'  ,EQtime/(iter-1));
fprintf(1,'         coefficients update    = %1.3e sec\n'  ,UDtime/(iter-1));
fprintf(1,'         fluid-mechanics solve  = %1.3e sec\n\n',FMtime/(iter-1));

fprintf(1,'         min T   =  %4.1f;    mean T   = %4.1f;    max T   = %4.1f;   [degC]\n' ,min(T(:)-273.15),mean(T(:)-273.15),max(T(:)-273.15));
fprintf(1,'         min SiO2=  %1.4f;    mean SiO2= %1.4f;    max SiO2= %1.4f;   [wt]\n'   ,min(c_oxd(:,:,cal.Si)./sum(c_oxd(:,:,1:end-1),3),[],'all'),mean(c_oxd(:,:,cal.Si)./sum(c_oxd(:,:,1:end-1),3),'all'),max(c_oxd(:,:,cal.Si)./sum(c_oxd(:,:,1:end-1),3),[],'all'));
fprintf(1,'         min H2O =  %1.4f;    mean H2O = %1.4f;    max H2O = %1.4f;   [wt]\n\n' ,min(c_oxd(:,:,cal.H ),[],'all'),mean(c_oxd(:,:,cal.H ),'all'),max(c_oxd(:,:,cal.H ),[],'all'));

fprintf(1,'         min x   =  %1.4f;    mean x   = %1.4f;    max x   = %1.4f;   [wt]\n'   ,min(x(:)  ),mean(x(:)  ),max(x(:)  ));
fprintf(1,'         min f   =  %1.4f;    mean f   = %1.4f;    max f   = %1.4f;   [wt]\n'   ,min(f(:)  ),mean(f(:)  ),max(f(:)  ));
fprintf(1,'         min m   =  %1.4f;    mean m   = %1.4f;    max m   = %1.4f;   [wt]\n\n' ,min(m(:)  ),mean(m(:)  ),max(m(:)  ));

fprintf(1,'         min rho =  %4.1f;    mean rho = %4.1f;    max rho = %4.1f;   [kg/m3]\n'  ,min(rho(:)),mean(rho(:))   ,max(rho(:)));
fprintf(1,'         min eta =  %1.2e;  mean eta = %1.2e;  max eta = %1.2e; [Pas]\n\n',min(eta(:)),geomean(eta(:)),max(eta(:)));

fprintf(1,'         min U   = %1.4e;    mean U   = %1.4e;    max U   = %1.4e;   [m/s]\n'  ,min(U(:)  ),mean(U(:)  ),max(U(:)  ));
fprintf(1,'         min W   = %1.4e;    mean W   = %1.4e;    max W   = %1.4e;   [m/s]\n'  ,min(-W(:) ),mean(-W(:) ),max(-W(:) ));
fprintf(1,'         min P   =  %1.4e;    mean P   =  %1.4e;    max P   = %1.4e;  [Pa]\n\n'  ,min(Pt(:)),mean(Pt(:)),max(Pt(:)));
