function varargout = ConstrFuncs(varargin)
% functions to apply specified constraints to proposed parameter sets

[varargout{1:nargout}] = feval(varargin{:});
end


function [model] = SumConstr(model, nc, ne, scale)

model = reshape(max(0,model).',nc,ne,[]);
model = model./sum(model,2)*scale;
model = reshape(model,nc*ne,[]).';

end

function [model] = NoConstr(model)
model = model.'; % for internal consistency
end