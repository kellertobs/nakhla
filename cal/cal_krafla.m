% set melting model end-member oxide compositions
cal.oxdStr =  {'SiO$_2$','Al$_2$O$_3$','FeO','MgO','CaO','Na$_2$O','K$_2$O'};
cal.oxds   = [  42.7    0.0    0.0   57.3    0.0    0.0    0.0     % cphs0 => dunite
                51.4   14.3   13.4    6.2   11.8    2.5    0.4     % perCm => basalt
                77.0   10.0    3.1    0.4    5.1    1.1    3.3];   % cphs1 => rhyolite
cal.oxds = cal.oxds./sum(cal.oxds,2)*100;

% specify phase diagram parameters
cal.cphs0    =  0.42;                % phase diagram lower bound composition [wt SiO2]
cal.cphs1    =  0.77;                % phase diagram upper bound composition [wt SiO2]
cal.Tphs0    =  845;                 % phase diagram lower bound temperature [degC]
cal.Tphs1    =  1870;                % phase diagram upper bound temperature [degC]
cal.PhDg     =  [7.0,4.0,1.0,1.0];   % Phase diagram curvature factor (> 1)
cal.perCm    =  0.517;               % peritectic liquidus composition [wt SiO2]
cal.perCx    =  0.477;               % peritectic solidus  composition [wt SiO2]
cal.perT     =  1145;                % peritectic temperature [degC]
cal.clap     =  1e-7;                % Clapeyron slope for P-dependence of melting T [degC/Pa]
cal.dTH2O    =  [1300,1100,900];     % solidus shift from water content [degC/wt^0.75]

