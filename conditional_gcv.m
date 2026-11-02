function score = conditional_gcv(lambda,c,s,bhat)
% GCV at fixed blur, using the exact trace of the Tikhonov filter.
% Parseval's factor 1/N cancels the normalization in the GCV denominator.
ell = lambda^2*abs(s(:)).^2;
one_minus_h = ell./(abs(c(:)).^2+ell);
score = sum(abs(one_minus_h.*bhat(:)).^2)/sum(one_minus_h)^2;
end
