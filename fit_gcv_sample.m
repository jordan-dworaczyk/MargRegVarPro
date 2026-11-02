function [y,lambda,score,history,alternation_converged] = fit_gcv_sample(alpha,y,bhat,s,eta)
% alternation_converged records only the relative y/lambda change tests.
% It does not certify an inner gradient test; history(:,7) stores inner stops.
% Alternate conditional GCV and blur estimation at a fixed prior exponent.
% The noise estimate is fixed, and gamma=1/eta for auxiliary delta=1, L=I.
validateattributes(y,{'double','single'},{'real','vector','numel',2,'finite','positive'},mfilename,'y');
y=y(:);
validateattributes(alpha,{'numeric'},{'real','scalar','finite','nonnegative'},mfilename,'alpha');
validateattributes(eta,{'numeric'},{'real','scalar','finite','positive'},mfilename,'eta');
n=size(bhat,1);
symbol=@(y) gauss_psf_symbol_derivs(y,n);
valid=@(y) all(isfinite(y)) && all(y>0);
beta=sqrt(alpha/eta); gamma=1/eta;
lambda_old=nan; alternation_converged=false;
history=zeros(200,7);
for k=1:200
    c=symbol(y);
    lambda=fminbnd(@(lambda) conditional_gcv(lambda,c,s,bhat), ...
        1e-8,1e3,optimset('TolX',1e-8,'Display','off'));
    [y_new,~,~,inner_stop]=solve_gcv_blur(y,bhat,s,lambda, ...
        beta,gamma,ones(n),symbol,valid,eta);
    relative_y_change=norm(y_new-y)/max(1,norm(y));
    relative_lambda_change=abs(lambda-lambda_old)/max(lambda,eps);
    history(k,:)=[k,y_new.',lambda,relative_y_change,relative_lambda_change,inner_stop];
    y=y_new;
    if relative_y_change<1e-4 && relative_lambda_change<1e-3
        alternation_converged=true;
        break
    end
    lambda_old=lambda;
end
history=history(1:k,:);
c=symbol(y);
[lambda,score]=fminbnd(@(lambda) conditional_gcv(lambda,c,s,bhat), ...
    1e-8,1e3,optimset('TolX',1e-8,'Display','off'));
end
