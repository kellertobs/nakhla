function  [H2Osat]  =  fluidsat(T,P,SiO2,cal)


cal.H2Osat_maf = (4.7773e-7.*P.^0.6 + 1e-11.*P) .* exp(2565*(1./(T+273.15)-1./(1473.15))); % Katz et al., 2003; Moore et al., 1998
cal.H2Osat_fls = (3.5494e-3.*P.^0.5 + 9.623e-8.*P - 1.5223e-11.*P.^1.5)./(T+273.15) + 1.2436e-14.*P.^1.5; % Liu et al., 2015

c      = max(0,min(1,(SiO2-45)./(75-45)));
H2Osat = (1-c).*cal.H2Osat_maf + c.*cal.H2Osat_fls;


end