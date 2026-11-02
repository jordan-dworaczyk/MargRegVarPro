function run_gcv_single_sample
% A single illustrative automatic-selection experiment, not a population study.
folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(folder,'output','gcv_single_sample');
if ~exist(output_folder,'dir'), mkdir(output_folder); end

%% Fixed sample: Cameraman, moderate Gaussian blur, 5% noise
x_true=im2double(imread(fullfile(folder,'cameraman_reconstruction.png')));
n=size(x_true,1); N=n^2;
y_true=[3;2]; y_start=[4;4];
previous_rng=rng;
rng_cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
seed=20271124; rng(seed,'twister'); % Fixed realization used in the paper.
c_true=reshape(gauss_psf_symbol_derivs(y_true,n),n,n);
psf_true=fftshift(real(ifft2(c_true)));
b_true=real(ifft2(c_true.*fft2(x_true)));
noise_level=.05;
noise=randn(n); noise=noise/norm(noise,'fro');
b=b_true+noise_level*norm(b_true,'fro')*noise;
bhat=fft2(b);
eta_reference=N/(noise_level^2*norm(b_true,'fro')^2);

%% Difference MAD noise estimate; first-difference image regularization
d=b(2:end,2:end)-b(1:end-1,2:end)-b(2:end,1:end-1)+b(1:end-1,1:end-1);
sigma_estimate=median(abs(d(:)-median(d(:))))/(2*.674489750196082);
eta=sigma_estimate^-2;
[u,v]=ndgrid(0:n-1);
s=sqrt(4*sin(pi*u/n).^2+4*sin(pi*v/n).^2);

%% Select alpha from data, with beta^2=alpha/eta and gamma=1/eta
timer=tic;
[alpha_selected,y_selected,lambda_selected,score_selected,trace, ...
    selected_history,selected_alternation_converged,search]=select_gcv_strength(y_start,bhat,s,eta);
[y_control,lambda_control,score_control,control_history,control_alternation_converged]= ...
    fit_gcv_sample(0,y_start,bhat,s,eta);
seconds=toc(timer);
alpha=[0 alpha_selected]; beta=sqrt(alpha/eta); gamma=1/eta;
y_estimate=[y_control y_selected]; lambda=[lambda_control lambda_selected];
gcv_score=[score_control score_selected];
alternation_converged=[control_alternation_converged selected_alternation_converged];
% Both convergence fields record the outer relative-change test.
converged=alternation_converged;

%% Ground truth is used only to evaluate the two final fits
x_estimate=zeros(n,n,2); psf_estimate=x_estimate;
parameter_error=zeros(1,2); RRE=parameter_error; PSF_error=parameter_error;
PSNR=parameter_error; SSIM=parameter_error;
for j=1:2
    c=reshape(gauss_psf_symbol_derivs(y_estimate(:,j),n),n,n);
    x=real(ifft2(conj(c).*bhat./(abs(c).^2+lambda(j)^2*s.^2)));
    psf=fftshift(real(ifft2(c)));
    x_estimate(:,:,j)=x; psf_estimate(:,:,j)=psf;
    parameter_error(j)=norm(y_estimate(:,j)-y_true);
    RRE(j)=norm(x-x_true,'fro')/norm(x_true,'fro');
    PSF_error(j)=norm(psf-psf_true,'fro')/norm(psf_true,'fro');
    PSNR(j)=10*log10(1/mean((x(:)-x_true(:)).^2));
    SSIM(j)=ssim(x,x_true);
end
T=table(["No log-det";"Automatic selection"],alpha.',beta.', ...
    lambda.',y_estimate(1,:).',y_estimate(2,:).',gcv_score.', ...
    RRE.',PSF_error.',PSNR.',SSIM.',alternation_converged.','VariableNames', ...
    {'method','alpha','beta','lambda','y1','y2','GCV','RRE','PSF_error','PSNR','SSIM','alternation_converged'});
disp(T);
writetable(T,fullfile(output_folder,'results.csv'));
Trace= array2table(trace,'VariableNames', ...
    {'alpha','GCV','y1','y2','lambda','alternations','alternation_converged'});
writetable(Trace,fullfile(output_folder,'search_trace.csv'));
save(fullfile(output_folder,'results.mat'),'x_true','b_true','b','bhat', ...
    'y_true','y_start','seed','noise_level','eta_reference','eta','sigma_estimate', ...
    's','alpha','beta','gamma','lambda','y_estimate','gcv_score','converged','alternation_converged', ...
    'x_estimate','psf_true','psf_estimate','parameter_error','RRE','PSF_error', ...
    'PSNR','SSIM','trace','selected_history','control_history','search','seconds','T');
fprintf('Noise precision ratio %.6g; elapsed %.1fs; search exit %d\n', ...
    eta/eta_reference,seconds,search.exitflag);
% Run plot_gcv_single_sample to export the paper figure from these results.
end
