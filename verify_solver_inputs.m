function verify_solver_inputs
% Check solver input handling and export-figure resource management.
n=8; bhat=fft2(reshape(sin((1:n^2)/7),n,n)); s=ones(n);
symbol=@(y) gauss_psf_symbol_derivs(y,n); valid=@(y) all(y>0);
solvers={@solve_reconstruction,@solve_auxiliary_scale,@solve_gcv_blur};
for k=1:numel(solvers)
    args={bhat,s,.1,.05,.01,s,symbol,valid};
    if k==3, args{end+1}=100; end
    [yc,hc,ic,sc]=solvers{k}([2;3],args{:});
    [yr,hr,ir,sr]=solvers{k}([2 3],args{:});
    assert(isequal(yc,yr) && isequal(hc,hr) && ic==ir && sc==sr);
    must_error(@() solvers{k}([-1;2],args{:}),'MargRegVarPro:InvalidInitialPoint');
    must_error(@() solvers{k}([NaN;2],args{:}),'');
    must_error(@() solvers{k}([Inf;2],args{:}),'');
    must_error(@() solvers{k}(ones(2),args{:}),'');
    bad=args; bad{1}=bhat(1:end-1,:);
    must_error(@() solvers{k}([2;3],bad{:}),'');
    bad=args; bad{7}=@bad_derivatives;
    must_error(@() solvers{k}([2;3],bad{:}),'MargRegVarPro:InvalidInitialObjective');
    if k==3
        bad=args; bad{end}=0;
        must_error(@() solvers{k}([2;3],bad{:}),'');
    else
        must_error(@() solvers{k}([2;3],args{:},0),'');
    end
end

% A failing export must release its own hidden figure, preserving others.
sentinel=figure('Visible','off','Name','Existing user figure');
sentinel_cleanup=onCleanup(@() delete(sentinel)); %#ok<NASGU>
before=findall(groot,'Type','figure');
must_error(@failing_export,'MargRegVarPro:TestExportError');
assert(all(isgraphics(before)) && numel(before)==numel(findall(groot,'Type','figure')));
fprintf('PASS: solver input errors, row/column equivalence, and figure cleanup.\n');
end

function must_error(callback,identifier)
try
    callback();
catch err
    if ~isempty(identifier), assert(strcmp(err.identifier,identifier),err.message); end
    return
end
error('MargRegVarPro:MissingExpectedError','Expected invalid input to be rejected.');
end

function [c,dc,d2c]=bad_derivatives(~)
c=ones(64,1); dc=nan(64,2); d2c=zeros(64,2,2);
end

function failing_export
[~,cleanup]=output_figure('Visible','off'); %#ok<ASGLU>
error('MargRegVarPro:TestExportError','Simulated failed export.');
end
