function [f,g,B,c] = reconstruction_objective(y,bhat,s,lambda,beta,gamma,s_aux,symbol,valid)
% Residual-scaled reduced objective, gradient, and modified Newton matrix.
% Uniform base prior. symbol returns c, its gradient, and its Hessian.
N = numel(bhat);
p = numel(y);
if ~all(isfinite(y)) || ~valid(y)
    f = inf; g = nan(p,1); B = nan(p); c = [];
    return
end
[c,dc,d2c] = symbol(y);
bhat = bhat(:); s = s(:); s_aux = s_aux(:);
a = abs(c).^2;
ell = lambda^2 * abs(s).^2;
den = a + ell;
q = a + gamma * abs(s_aux).^2;
if any(den<=0) || any(q<=0) || any(~isfinite(den)) || any(~isfinite(q))
    f = inf; g = nan(p,1); B = nan(p);
    return
end

% VarPro term: Parseval normalization for MATLAB's unnormalized fft2.
r = lambda * abs(s) .* abs(bhat) ./ sqrt(den);
da_half = real(conj(c) .* dc);
J = -(lambda * abs(s) .* abs(bhat) ./ den.^1.5) .* da_half;
varpro = .5 * sum(r.^2) / N;
g_varpro = (J.' * r) / N;
B_varpro = (J.' * J) / N;

% Exact first and second derivatives of the auxiliary log-determinant.
dq = 2 * da_half;
g_logdet = sum(dq ./ q,1).';
H_logdet = zeros(p);
for i = 1:p
    for j = 1:p
        d2q = 2 * real(conj(dc(:,i)).*dc(:,j) + conj(c).*d2c(:,i,j));
        H_logdet(i,j) = sum(d2q./q - dq(:,i).*dq(:,j)./q.^2);
    end
end
f = varpro + beta^2/2 * sum(log(q));
g = g_varpro + beta^2/2 * g_logdet;
B = B_varpro + beta^2/2 * H_logdet;
end
