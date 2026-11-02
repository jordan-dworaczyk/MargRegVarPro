function [alpha,y,lambda,score,trace,history,alternation_converged,search] = select_gcv_strength(y_start,bhat,s,eta)
% alternation_converged records only the relative y/lambda change tests.
% It does not certify an inner gradient test; history(:,7) stores inner stops.
% One bounded scalar search. Every trial uses the same initial guess for the blur parameters.
% Ground truth and reconstruction errors are absent from the selection.
trace=zeros(0,7);
options=optimset('TolX',2e-3,'MaxFunEvals',30,'MaxIter',30,'Display','off');
[alpha,~,exitflag,search]=fminbnd(@nested_score,0,1,options);
search.exitflag=exitflag;
[y,lambda,score,history,alternation_converged]=fit_gcv_sample(alpha,y_start,bhat,s,eta);

    function value=nested_score(trial_alpha)
        [trial_y,trial_lambda,value,trial_history,trial_alternation_converged]= ...
            fit_gcv_sample(trial_alpha,y_start,bhat,s,eta);
        trace(end+1,:)=[trial_alpha,value,trial_y.',trial_lambda, ...
            size(trial_history,1),trial_alternation_converged];
        fprintf('alpha=%.6g, GCV=%.6g, y_alpha*=[%.5g %.5g], alternations=%d, alternation converged=%d\n', ...
            trial_alpha,value,trial_y,size(trial_history,1),trial_alternation_converged);
    end
end
