function [mu, dmu, d2mu] = defocus_psf_symbol_derivs(y, n)
%DEFOCUS_PSF_SYMBOL_DERIVS
% Fourier symbol and exact first/second derivatives for a smoothed
% defocus PSF with fixed edge smoothing.
%
% Parameter:
%   y = [R]
%   R = blur radius in pixels
%
% Model:
%   h(r;R) = 1 / (1 + exp((r - R)/eps_edge))
%   p      = h / sum(h)
%
% Output:
%   mu    : n^2 x 1
%   dmu   : n^2 x 1
%   d2mu  : n^2 x 1 x 1

    y = y(:);
    if numel(y) ~= 1
        error('defocus_psf_symbol_derivs: y must be [R].');
    end

    R = y(1);
    eps_edge = 5e-1;   % fixed smoothing width in pixels

    if ~isscalar(n) || n <= 0 || n ~= round(n)
        error('n must be a positive integer scalar.');
    end
    % Positivity and finiteness are checked before the objective calls this function.

    center = [floor(n/2) + 1, floor(n/2) + 1];
    [xg, yg] = ndgrid(-fix(n/2):ceil(n/2)-1, -fix(n/2):ceil(n/2)-1);
    r = sqrt(xg.^2 + yg.^2);

    % Base soft disk
    a = (r - R) / eps_edge;
    h = 1 ./ (1 + exp(a));

    % Logistic helper
    q = h .* (1 - h);

    % First derivative wrt R
    hR = (1 / eps_edge) * q;

    % Second derivative wrt R
    hRR = (1 / eps_edge^2) * q .* (1 - 2*h);

    % Normalize p = h / Z
    Z   = sum(h(:));
    ZR  = sum(hR(:));
    ZRR = sum(hRR(:));

    p = h / Z;

    % First derivative of normalized PSF
    pR = hR / Z - h * (ZR / Z^2);

    % Second derivative of normalized PSF
    pRR = hRR / Z ...
        - 2 * hR * (ZR / Z^2) ...
        - h  * (ZRR / Z^2) ...
        + 2 * h * (ZR^2 / Z^3);

    % Shift center to (1,1), then FFT
    p_shift   = circshift(p,   1 - center);
    pR_shift  = circshift(pR,  1 - center);
    pRR_shift = circshift(pRR, 1 - center);

    mu_mat  = fft2(p_shift);
    dR_mat  = fft2(pR_shift);
    dRR_mat = fft2(pRR_shift);

    mu = mu_mat(:);

    dmu = zeros(n*n, 1);
    dmu(:,1) = dR_mat(:);

    d2mu = zeros(n*n, 1, 1);
    d2mu(:,1,1) = dRR_mat(:);
end
