function y=validate_solver_inputs(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid)
% Validate solver inputs and return the initial parameters as a column vector.
validateattributes(y,{'double','single'},{'real','vector','nonempty','finite'},mfilename,'y');
y=y(:);
validateattributes(bhat,{'numeric'},{'nonempty','finite'},mfilename,'bhat');
validateattributes(s,{'numeric'},{'numel',numel(bhat),'finite'},mfilename,'s');
validateattributes(s_aux,{'numeric'},{'numel',numel(bhat),'finite'},mfilename,'s_aux');
validateattributes(lambda,{'numeric'},{'real','scalar','finite','positive'},mfilename,'lambda');
validateattributes(beta,{'numeric'},{'real','scalar','finite','nonnegative'},mfilename,'beta');
validateattributes(gamma,{'numeric'},{'real','scalar','finite','positive'},mfilename,'gamma');
assert(isa(symbol,'function_handle') && isa(valid,'function_handle'), ...
    'MargRegVarPro:InvalidCallback','symbol and valid must be function handles.');
feasible=valid(y);
assert(islogical(feasible) && isscalar(feasible) && feasible, ...
    'MargRegVarPro:InvalidInitialPoint', ...
    'The initial blur parameters must satisfy valid(y), returning scalar logical true.');
end
