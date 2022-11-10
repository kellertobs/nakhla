% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','TiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O'};
cal.oxd    = [  42.7    0.0    0.0    0.0   57.3    0.0    0.0    0.0     % for
                30.0    0.5    0.0   69.5    0.0    0.0    0.0    0.0     % fay
                55.0    0.05   7.5    3.0   33.95   0.4    0.1    0.0     % opx
                49.0    0.1    2.0   27.6    5.1   16.0    0.2    0.0     % cpx
                44.4    0.0   35.8    0.0    0.0   19.2    0.6    0.0     % ant
                67.3    0.0   20.2    0.0    0.0    0.8   11.2    0.5     % alb
                64.8    0.0   18.3    0.0    0.0    0.0    0.0   16.9     % kfs
                99.99   0.01   0.0    0.0    0.0    0.0    0.0    0.0];   % qtz

cal.cmpStr = {'for','fay','opx','cpx','ant','alb','kfs','qtz'};
cal.cmp    = [  95.0    5.0    0.0    0.0    0.0    0.0    0.0    0.0     % cphs0 => dunite
                12.0    8.0   34.0   20.0   24.0    2.0    0.0    0.0     % perCx => pyroxenite
                 1.0    1.5   15.0   34.0   28.0   20.0    0.5    0.0     % perCm => basalt
                 0.0    0.0    0.0    8.0    9.0   15.0   31.0   37.0];   % cphs1 => rhyolite

cal.oxd = cal.oxd./sum(cal.oxd,2)*100;
cal.cmp = cal.cmp./sum(cal.cmp,2)*100;

cal.cmp_oxd = cal.cmp*cal.oxd./100;

cal.rhox0 = [3270,4390,3000,3250,2730,2620,2580,2650];
cal.rhom0 = [2710,3580,2580,2850,2530,2310,2290,2360];
cal.rhof0 = 500;

% specify phase diagram parameters
cal.cphs0    =  cal.cmp_oxd(1,1)/100; % phase diagram lower bound composition [wt SiO2]
cal.cphs1    =  cal.cmp_oxd(4,1)/100; % phase diagram upper bound composition [wt SiO2]
cal.Tphs0    =  845;                  % phase diagram lower bound temperature [degC]
cal.Tphs1    =  1780;                 % phase diagram upper bound temperature [degC]
cal.PhDg     =  [8.0,4.0,1.2,1.2];    % phase diagram curvature factor (> 1)
cal.perCm    =  cal.cmp_oxd(3,1)/100; % peritectic liquidus composition [wt SiO2]
cal.perCx    =  cal.cmp_oxd(2,1)/100; % peritectic solidus  composition [wt SiO2]
cal.perT     =  1145;                 % peritectic temperature [degC]
cal.clap     =  1e-7;                 % Clapeyron slope for P-dependence of melting T [degC/Pa]
cal.dTH2O    =  [1300,1100,900];      % solidus shift from water content [degC/wt^0.75]

