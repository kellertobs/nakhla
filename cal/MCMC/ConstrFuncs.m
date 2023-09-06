function varargout = ConstrFuncs(varargin)
% functions to apply specified constraints to proposed parameter sets

[varargout{1:nargout}] = feval(varargin{:});
end


function [model] = SumConstr(model, nc, ne, scale)

wk = reshape(max(0,model(1:nc*ne)).',nc,ne,[]);
wk = wk./sum(wk,2)*scale;
model(1:nc*ne) = reshape(wk,nc*ne,[]).';

end

function [model] = NoConstr(model)
model = model;
end