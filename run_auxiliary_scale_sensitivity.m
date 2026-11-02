function run_auxiliary_scale_sensitivity
% Defocus auxiliary-scale sensitivity with a fixed determinant strength.
previous_rng=rng;
rng_cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(0,'twister');
experiment_folder = fileparts(mfilename('fullpath'));
output_folder = fullfile(experiment_folder,'output','auxiliary_scale');
if ~exist(output_folder,'dir'), mkdir(output_folder); end

%% Image, smooth defocus kernel, and observed data
x_true = im2double(imread(fullfile(experiment_folder,'cameraman_reconstruction.png')));
n = size(x_true,1); N = n^2;
y_true = 3; y_start = 3.5;
[row,col] = ndgrid(-n/2:n/2-1);
radius = sqrt(row.^2+col.^2);
h = 1./(1+exp((radius-y_true)/.5));
psf_true = h/sum(h(:));
c_true = fft2(ifftshift(psf_true));
b_true = real(ifft2(c_true.*fft2(x_true)));
noise = randn(n); noise = noise/norm(noise,'fro');
b = b_true + .01*norm(b_true,'fro')*noise;
bhat = fft2(b);
sigma_noise = .01*norm(b_true,'fro')/sqrt(N);
eta = 1/sigma_noise^2;
delta = 100; lambda = sqrt(delta/eta);
omega = 2*pi/n*(0:n-1);
d2 = 2*cos(omega)-2;
s = sqrt(d2(:).^2+d2.^2);

%% Auxiliary prior: fixed strength, identity operator, varying scale
beta = 1/sqrt(N); alpha = eta*beta^2;
tilde_delta = 10.^(-12:12);
gamma = tilde_delta/eta;
symbol = @(y) defocus_psf_symbol_derivs(y,n);
valid = @(y) all(isfinite(y)) && all(y>0);
y_estimate = zeros(size(tilde_delta));
RRE = y_estimate; PSF_error = y_estimate;
iterations = y_estimate; stop = y_estimate; history = cell(size(tilde_delta));

%% Solve and evaluate the reconstruction at each auxiliary scale
for j = 1:numel(tilde_delta)
    [y,history{j},iterations(j),stop(j)] = solve_auxiliary_scale( ...
        y_start,bhat,s,lambda,beta,gamma(j),ones(n),symbol,valid,1e-10);
    y_estimate(j) = y;
    c = reshape(symbol(y),n,n);
    x = real(ifft2(conj(c).*bhat./(abs(c).^2+lambda^2*s.^2)));
    psf = fftshift(real(ifft2(c)));
    RRE(j) = norm(x-x_true,'fro')/norm(x_true,'fro');
    PSF_error(j) = norm(psf-psf_true,'fro')/norm(psf_true,'fro');
    fprintf('tilde delta=%g, y=%.9g, RRE=%.6g, PSF error=%.6g, stop=%d\n', ...
        tilde_delta(j),y,RRE(j),PSF_error(j),stop(j));
end
parameter_error = abs(y_estimate-y_true);
save(fullfile(output_folder,'results.mat'),'x_true','b_true','b','bhat', ...
    'eta','delta','lambda','s','beta','alpha','gamma','tilde_delta', ...
    'y_true','y_start','y_estimate','parameter_error','RRE','PSF_error', ...
    'history','iterations','stop');
T = table(tilde_delta.',y_estimate.',parameter_error.',RRE.',PSF_error.', ...
    iterations.',stop.','VariableNames', ...
    {'tilde_delta','y_estimate','parameter_error','RRE','PSF_error','iterations','stop'});
writetable(T,fullfile(output_folder,'results.csv'));

disp(T);
end
