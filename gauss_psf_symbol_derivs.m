function [mu, dmu, d2mu] = gauss_psf_symbol_derivs(sigma, n)
%GAUSS_PSF_SYMBOL_DERIVS
% Analytic Fourier symbol and its first/second derivatives for a
% normalized anisotropic Gaussian PSF.
%
% INPUTS
%   sigma : [sigma1; sigma2] or [sigma1, sigma2], with sigma1,sigma2 > 0
%   n     : image size (assumed square n-by-n)
%
% OUTPUTS
%   mu    : n^2 x 1
%           Fourier symbol of the normalized PSF
%   dmu   : n^2 x 2
%           dmu(:,1) = d mu / d sigma1
%           dmu(:,2) = d mu / d sigma2
%   d2mu  : n^2 x 2 x 2
%           d2mu(:,s,r) = d^2 mu / (d sigma_s d sigma_r)
%
% PSF model:
%   g(x,y) = exp( -x^2/sigma1^2 - y^2/sigma2^2 )
%   p      = g / sum(g(:))
%
% The FFT symbol is computed from the centered PSF by shifting the PSF
% center to (1,1) before applying fft2.

    % ----------------------------
    % Input checks
    % ----------------------------
    sigma = sigma(:);
    if numel(sigma) ~= 2
        error('sigma must have exactly two entries: [sigma1; sigma2].');
    end
    if ~isscalar(n) || n <= 0 || n ~= round(n)
        error('n must be a positive integer scalar.');
    end

    s1 = sigma(1);
    s2 = sigma(2);
    % Positivity and finiteness are checked before the objective calls this function.

    % ----------------------------
    % Spatial grid and PSF center
    % ----------------------------
    center = [floor(n/2) + 1, floor(n/2) + 1];
    [xg, yg] = ndgrid(-fix(n/2):ceil(n/2)-1, -fix(n/2):ceil(n/2)-1);

    % ----------------------------
    % Unnormalized Gaussian PSF
    % ----------------------------
    g = exp(-(xg.^2) / (s1^2) - (yg.^2) / (s2^2));

    % ----------------------------
    % First derivatives of g
    % ----------------------------
    % dg/ds1 = g * (2 x^2 / s1^3)
    % dg/ds2 = g * (2 y^2 / s2^3)
    g1 = g .* (2 * xg.^2 / s1^3);
    g2 = g .* (2 * yg.^2 / s2^3);

    % ----------------------------
    % Second derivatives of g
    % ----------------------------
    % d2g/ds1^2 = g * (4 x^4 / s1^6 - 6 x^2 / s1^4)
    % d2g/ds2^2 = g * (4 y^4 / s2^6 - 6 y^2 / s2^4)
    % d2g/(ds1 ds2) = g * (4 x^2 y^2 / (s1^3 s2^3))
    g11 = g .* (4 * xg.^4 / s1^6 - 6 * xg.^2 / s1^4);
    g22 = g .* (4 * yg.^4 / s2^6 - 6 * yg.^2 / s2^4);
    g12 = g .* (4 * xg.^2 .* yg.^2 / (s1^3 * s2^3));

    % ----------------------------
    % Normalization and derivatives
    % ----------------------------
    Z   = sum(g(:));
    Z1  = sum(g1(:));
    Z2  = sum(g2(:));
    Z11 = sum(g11(:));
    Z22 = sum(g22(:));
    Z12 = sum(g12(:));

    % ----------------------------
    % Normalized PSF p = g / Z
    % ----------------------------
    p = g / Z;

    % First derivatives of normalized PSF
    p1 = g1 / Z - g .* (Z1 / Z^2);
    p2 = g2 / Z - g .* (Z2 / Z^2);

    % Second derivatives of normalized PSF
    p11 = g11 / Z ...
        - g1 .* (Z1 / Z^2) ...
        - g1 .* (Z1 / Z^2) ...
        - g  .* (Z11 / Z^2) ...
        + 2 * g .* (Z1 * Z1 / Z^3);

    p22 = g22 / Z ...
        - g2 .* (Z2 / Z^2) ...
        - g2 .* (Z2 / Z^2) ...
        - g  .* (Z22 / Z^2) ...
        + 2 * g .* (Z2 * Z2 / Z^3);

    p12 = g12 / Z ...
        - g1 .* (Z2 / Z^2) ...
        - g2 .* (Z1 / Z^2) ...
        - g  .* (Z12 / Z^2) ...
        + 2 * g .* (Z1 * Z2 / Z^3);

    % ----------------------------
    % Shift PSF center to (1,1), then FFT
    % ----------------------------
    p_shift   = circshift(p,   1 - center);
    p1_shift  = circshift(p1,  1 - center);
    p2_shift  = circshift(p2,  1 - center);
    p11_shift = circshift(p11, 1 - center);
    p22_shift = circshift(p22, 1 - center);
    p12_shift = circshift(p12, 1 - center);

    mu_mat  = fft2(p_shift);
    d1_mat  = fft2(p1_shift);
    d2_mat  = fft2(p2_shift);
    d11_mat = fft2(p11_shift);
    d22_mat = fft2(p22_shift);
    d12_mat = fft2(p12_shift);

    % ----------------------------
    % Pack outputs
    % ----------------------------
    mu = mu_mat(:);

    dmu = zeros(n*n, 2);
    dmu(:,1) = d1_mat(:);
    dmu(:,2) = d2_mat(:);

    d2mu = zeros(n*n, 2, 2);
    d2mu(:,1,1) = d11_mat(:);
    d2mu(:,2,2) = d22_mat(:);
    d2mu(:,1,2) = d12_mat(:);
    d2mu(:,2,1) = d12_mat(:);
end
