function [mu, dmu, d2mu] = motion_psf_symbol_derivs(y, n)
%MOTION_PSF_SYMBOL_DERIVS
% Analytic Fourier symbol and first/second derivatives for a smooth motion PSF.
%
% Convention:
%   y(1) = L      > 0
%   y(2) = theta
%
%   theta is measured from the horizontal image axis, matching
%   fspecial('motion', L, theta_deg).
%
% Smooth model:
%   u = cos(theta)*col - sin(theta)*row
%   v = sin(theta)*col + cos(theta)*row
%
%   q(v) = exp(-v^2/(2*sigma_perp^2))
%   w(u;L) = 0.5 * [tanh((u+L/2)/tau) - tanh((u-L/2)/tau)]
%   g = q .* w
%   p = g / sum(g(:))

    y = y(:);
    if numel(y) ~= 2
        error('motion_psf_symbol_derivs: y must be [L; theta].');
    end
    if ~isscalar(n) || n <= 0 || n ~= round(n)
        error('n must be a positive integer scalar.');
    end

    L = y(1);
    theta = y(2);

    if ~(isfinite(L) && isfinite(theta)) || L <= 0
        error('motion_psf_symbol_derivs: require L > 0 and finite theta.');
    end

    % Fixed transverse width and endpoint smoothing used in the paper.
    sigma_perp = 0.30;
    tau        = 0.15;

    center = [floor(n/2) + 1, floor(n/2) + 1];

    [row, col] = ndgrid(-fix(n/2):ceil(n/2)-1, ...
                        -fix(n/2):ceil(n/2)-1);

    c = cos(theta);
    s = sin(theta);

    % Image-coordinate angle convention:
    % theta = 0 horizontal, theta = pi/2 vertical.
    u =  c * col - s * row;
    v =  s * col + c * row;

    % theta derivatives
    u_t  = -s * col - c * row;
    v_t  =  c * col - s * row;

    % Since u_tt = -u and v_tt = -v
    u_tt = -u;
    v_tt = -v;

    % Smooth line window
    aplus  = (u + L/2) / tau;
    aminus = (u - L/2) / tau;

    Tplus  = tanh(aplus);
    Tminus = tanh(aminus);

    Splus  = sech2_local(aplus);
    Sminus = sech2_local(aminus);

    w = 0.5 * (Tplus - Tminus);

    % Derivatives of w wrt u and L
    wu   = (Splus - Sminus) / (2 * tau);
    wuu  = -(Splus .* Tplus - Sminus .* Tminus) / (tau^2);

    wL   = (Splus + Sminus) / (4 * tau);
    wLL  = -(Splus .* Tplus - Sminus .* Tminus) / (4 * tau^2);

    wuL  = -(Splus .* Tplus + Sminus .* Tminus) / (2 * tau^2);

    % Gaussian cross-profile
    sig2 = sigma_perp^2;

    q = exp(-(v.^2) / (2 * sig2));

    qv  = -(v / sig2) .* q;
    qvv = ((v.^2) / sig2^2 - 1 / sig2) .* q;

    % Unnormalized PSF
    g = q .* w;

    % First derivatives
    gL = q .* wL;

    gT = qv .* v_t .* w + q .* wu .* u_t;

    % Second derivatives
    gLL = q .* wLL;

    gLT = qv .* v_t .* wL + q .* wuL .* u_t;

    gTT = qvv .* (v_t.^2) .* w ...
        + qv  .* v_tt .* w ...
        + 2 * qv .* v_t .* wu .* u_t ...
        + q   .* wuu .* (u_t.^2) ...
        + q   .* wu  .* u_tt;

    % Normalization
    Z   = sum(g(:));
    ZL  = sum(gL(:));
    ZT  = sum(gT(:));
    ZLL = sum(gLL(:));
    ZLT = sum(gLT(:));
    ZTT = sum(gTT(:));

    if ~isfinite(Z) || Z <= 0
        error('Invalid PSF normalization in motion_psf_symbol_derivs.');
    end

    % Normalized PSF
    p = g / Z;

    % First derivatives of normalized PSF
    pL = gL / Z - g .* (ZL / Z^2);
    pT = gT / Z - g .* (ZT / Z^2);

    % Second derivatives of normalized PSF
    pLL = gLL / Z ...
        - 2 * gL .* (ZL / Z^2) ...
        - g  .* (ZLL / Z^2) ...
        + 2 * g .* (ZL^2 / Z^3);

    pTT = gTT / Z ...
        - 2 * gT .* (ZT / Z^2) ...
        - g  .* (ZTT / Z^2) ...
        + 2 * g .* (ZT^2 / Z^3);

    pLT = gLT / Z ...
        - gL .* (ZT / Z^2) ...
        - gT .* (ZL / Z^2) ...
        - g  .* (ZLT / Z^2) ...
        + 2 * g .* (ZL * ZT / Z^3);

    % Shift center to (1,1), then FFT
    p_shift   = circshift(p,   1 - center);
    pL_shift  = circshift(pL,  1 - center);
    pT_shift  = circshift(pT,  1 - center);
    pLL_shift = circshift(pLL, 1 - center);
    pTT_shift = circshift(pTT, 1 - center);
    pLT_shift = circshift(pLT, 1 - center);

    mu_mat  = fft2(p_shift);
    dL_mat  = fft2(pL_shift);
    dT_mat  = fft2(pT_shift);
    dLL_mat = fft2(pLL_shift);
    dTT_mat = fft2(pTT_shift);
    dLT_mat = fft2(pLT_shift);

    mu = mu_mat(:);

    dmu = zeros(n*n, 2);
    dmu(:,1) = dL_mat(:);
    dmu(:,2) = dT_mat(:);

    d2mu = zeros(n*n, 2, 2);
    d2mu(:,1,1) = dLL_mat(:);
    d2mu(:,2,2) = dTT_mat(:);
    d2mu(:,1,2) = dLT_mat(:);
    d2mu(:,2,1) = dLT_mat(:);
end


function s = sech2_local(x)
    % Squared hyperbolic secant used in the analytic derivatives.
    s = 1 ./ (cosh(x).^2);
end
