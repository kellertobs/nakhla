function  [H2Osat]  =  fluidsat(var)

MW  = [60.0855, 79.88, 101.96, 71.85, 40.3, 56.08, 61.98, 94.2, 18.02];

P    = min(var.P,0.5)*1e9/1e5;
T    = var.T+273.15;
X    = var.X./MW ./ sum(var.X./MW,2);
X    = X./sum(X(:,1:end-1),2);
ln_fH2O = log(P);

a = 2565;
b = [-1.997,-0.9275,2.736];
c = 1.171;
d = -14.21;

XH2O   = exp((a./T + sum(b.*X(:,[3,4,7]),2).*(P./T) + c.*ln_fH2O + d)/2); % Moore et al., 1998

H2Osat = XH2O.*MW(end) ./ sum([X(:,1:end-1), XH2O].*MW,2);

end