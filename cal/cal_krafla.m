% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O'};
cal.oxds   = [  42.7    0.0    0.0   57.3    0.0    0.0    0.0     % for
                30.0    0.0   70.0    0.0    0.0    0.0    0.0     % fay
                55.0    7.5    3.0   34.0    0.4    0.1    0.0     % opx
                49.0    2.0   27.7    5.1   16.0    0.2    0.0     % cpx
                44.4   35.8    0.0    0.0   19.2    0.6    0.0     % ant
                67.3   20.2    0.0    0.0    0.8   11.2    0.5     % alb
                64.8   18.3    0.0    0.0    0.0    0.0   16.9     % kfs
               100.0    0.0    0.0    0.0    0.0    0.0    0.0];   % qtz

cal.cmpStr = {'for','fay','opx','cpx','ant','alb','kfs','qtz'};
cal.cmps   = [ 100.0    0.0    0.0    0.0    0.0    0.0    0.0    0.0     % cphs0 => dunite
                20.0   11.0   34.0   35.0    0.0    0.0    0.0    0.0     % perCx => pyroxenite
                 2.0    2.0   18.0   30.0   28.0   20.0    0.0    0.0     % perCm => basalt
                 0.0    0.0    0.0    8.0    9.0   15.0   31.0   37.0];   % cphs1 => rhyolite

cal.oxds = cal.oxds./sum(cal.oxds,2)*100;
cal.cmps = cal.cmps./sum(cal.cmps,2)*100;

p = cal.cmps*cal.oxds./100;

% specify phase diagram parameters
cal.cphs0    =  p(1,1)/100;          % phase diagram lower bound composition [wt SiO2]
cal.cphs1    =  p(4,1)/100;          % phase diagram upper bound composition [wt SiO2]
cal.Tphs0    =  845;                 % phase diagram lower bound temperature [degC]
cal.Tphs1    =  1870;                % phase diagram upper bound temperature [degC]
cal.PhDg     =  [7.0,4.0,1.0,1.0];   % Phase diagram curvature factor (> 1)
cal.perCm    =  p(3,1)/100;          % peritectic liquidus composition [wt SiO2]
cal.perCx    =  p(2,1)/100;          % peritectic solidus  composition [wt SiO2]
cal.perT     =  1145;                % peritectic temperature [degC]
cal.clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
cal.dTH2O    =  [1300,1100,900];     % solidus shift from water content [degC/wt^0.75]

