function run_reconstruction_defocus
% Defocus reconstruction study: three initial guesses and three determinant choices.
% Requires Image Processing Toolbox. All parameters are fixed in this study.
previous_rng=rng;
rng_cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(0,'twister');
experiment_folder = fileparts(mfilename('fullpath'));
output_folder = fullfile(experiment_folder,'output','reconstruction_defocus');
if ~exist(output_folder,'dir'), mkdir(output_folder); end

%% Image, blur parameters, and initial guesses
x_true = im2double(imread(fullfile(experiment_folder,'cameraman_reconstruction.png')));
n = size(x_true,1);
N = n^2;
y_true = 2.5;
y_start = [1.5 4 6];
symbol = @(y) defocus_psf_symbol_derivs(y,n);
valid = @(y) all(y>0);
[row,col] = ndgrid(-n/2:n/2-1);
eps_edge = .5;
g = 1./(1+exp((sqrt(row.^2+col.^2)-y_true)/eps_edge));

psf_true = g/sum(g(:));
c_true = fft2(circshift(psf_true,[-n/2,-n/2]));

%% Observations and image-prior parameters
b_true = real(ifft2(c_true.*fft2(x_true)));
noise_level = .01;
noise = randn(n); noise = noise/norm(noise,'fro');
b = b_true + noise_level*norm(b_true,'fro')*noise;
bhat = fft2(b);
sigma_noise = noise_level*norm(b_true,'fro')/sqrt(N);
eta = 1/sigma_noise^2;
delta = 100;
lambda = sqrt(delta/eta);

% Periodic second differences in both directions, in MATLAB FFT order.
omega = 2*pi/n*(0:n-1);
d2 = 2*cos(omega)-2;
s = sqrt(d2(:).^2+d2.^2);

%% The three objective choices, with a uniform base prior
% Order in saved arrays: No log-det, Proposed, Gaussian marginalization.
% Exact marginalization requires alpha=1 and the matched auxiliary precision.
% Proposed alpha minimizes mean true-image RRE over the three initial guesses on the oracle grid.
alpha = [0, 0.75, 1];
beta = sqrt(alpha/eta);
% Marginalization retains alpha=1; the proposed weight is shared across initial guesses.
gamma = [1/eta, 1/eta, delta/eta];
s_aux = cat(3,ones(n),ones(n),s);
method_names = {'No log-det','Proposed','Gaussian marginalization'};

% f = (1/(2N))*sum(ell*|bhat|^2/(|c|^2+ell))
%     + (beta^2/2)*sum(log(|c|^2+gamma*|s_aux|^2)).
% The solver evaluates this expression and its derivatives explicitly.
y_estimate = zeros(numel(y_true),3,3);
x_estimate = zeros(n,n,3,3);
psf_estimate = zeros(n,n,3,3);
history = cell(3,3);
iterations = zeros(3); stop = zeros(3);
parameter_error = zeros(3); RRE = zeros(3); PSF_error = zeros(3);
relative_residual = zeros(3); PSNR = zeros(3); SSIM = zeros(3);

%% Solve from each initial guess, then reconstruct x from the estimated blur
for start = 1:3
    for method = 1:3
        [y,history{method,start},iterations(method,start),stop(method,start)] = ...
            solve_reconstruction(y_start(:,start),bhat,s,lambda,beta(method), ...
                gamma(method),s_aux(:,:,method),symbol,valid);
        y_estimate(:,method,start) = y;
        c = reshape(symbol(y),n,n);
        xhat = conj(c).*bhat./(abs(c).^2+lambda^2*abs(s).^2);
        x = real(ifft2(xhat));
        psf = fftshift(real(ifft2(c)));
        fitted_data = real(ifft2(c.*fft2(x)));
        x_estimate(:,:,method,start) = x;
        psf_estimate(:,:,method,start) = psf;
        parameter_error(method,start) = norm(y-y_true);
        RRE(method,start) = norm(x-x_true,'fro')/norm(x_true,'fro');
        PSF_error(method,start) = norm(psf-psf_true,'fro')/norm(psf_true,'fro');
        relative_residual(method,start) = norm(fitted_data-b,'fro')/norm(b,'fro');
        PSNR(method,start) = 10*log10(1/mean((x(:)-x_true(:)).^2));
        SSIM(method,start) = ssim(x,x_true);
        fprintf('%s, initial guess %d: y*=%s, RRE=%.6g, PSF error=%.6g, stop=%d\n', ...
            method_names{method},start,mat2str(y.',6),RRE(method,start), ...
            PSF_error(method,start),stop(method,start));
    end
end
save(fullfile(output_folder,'results.mat'),'y_true','y_start','y_estimate', ...
    'x_true','b_true','b','bhat','psf_true','x_estimate','psf_estimate', ...
    'history','iterations','stop','parameter_error','RRE','PSF_error', ...
    'relative_residual','PSNR','SSIM','eta','delta','lambda','beta','alpha','gamma','s');

%% Display every reconstructed image; columns are initial guesses
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Position',[100 100 1000 850]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
for method = 1:3
    for start = 1:3
        nexttile; imagesc(x_estimate(:,:,method,start),[0 1]);
        axis image off; colormap gray;
        title(sprintf('%s, initial guess %d\nRRE %.3f, PSF error %.3f', ...
            method_names{method},start,RRE(method,start),PSF_error(method,start)));
    end
end
exportgraphics(fig,fullfile(output_folder,'all_reconstructions.png'),'Resolution',200);

%% Show all estimated PSFs on the same central crop
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Position',[100 100 1000 850]);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
crop = n/2+1+(-12:12);
for method = 1:3
    for start = 1:3
        nexttile; imagesc(psf_estimate(crop,crop,method,start));
        axis image off; colormap gray;
        title(sprintf('%s, initial guess %d',method_names{method},start));
    end
end
exportgraphics(fig,fullfile(output_folder,'all_psfs.png'),'Resolution',200);
end
