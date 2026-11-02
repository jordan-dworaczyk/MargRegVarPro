function run_alpha_sweep_defocus
% Defocus reference-image weight selection over the three initial guesses.
% Requires Image Processing Toolbox. All parameters are fixed in this study.
previous_rng=rng;
rng_cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(0,'twister');
experiment_folder = fileparts(mfilename('fullpath'));
output_folder = fullfile(experiment_folder,'output','alpha_sweep_defocus');
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

%% Candidate weights for the proposed auxiliary prior
% Choose one alpha for the entire blur family using the mean image RRE.
% All three initial guesses contribute, including runs that reach the iteration limit.
% Ground truth is used for this assessment, but is never passed to the solver.
alpha_grid = unique([.05 .1 .2 .35 .5 eta/N .75 1 1.5 2]);
gamma = 1/eta;
s_aux = ones(n);
number_of_weights = numel(alpha_grid);
RRE = zeros(number_of_weights,3);
PSF_error = RRE;
parameter_error = RRE;
iterations = RRE;
stop = RRE;
history = cell(number_of_weights,3);
y_estimate = zeros(numel(y_true),number_of_weights,3);

%% Solve from each initial guess for every candidate alpha
for start = 1:3
    for j = 1:number_of_weights
        beta = sqrt(alpha_grid(j)/eta);
        [y,history{j,start},iterations(j,start),stop(j,start)] = ...
            solve_reconstruction(y_start(:,start),bhat,s,lambda,beta, ...
                gamma,s_aux,symbol,valid);
        c = reshape(symbol(y),n,n);
        x = real(ifft2(conj(c).*bhat ./ (abs(c).^2+lambda^2*abs(s).^2)));
        psf = fftshift(real(ifft2(c)));
        y_estimate(:,j,start) = y;
        RRE(j,start) = norm(x-x_true,'fro')/norm(x_true,'fro');
        PSF_error(j,start) = norm(psf-psf_true,'fro')/norm(psf_true,'fro');
        parameter_error(j,start) = norm(y-y_true);
        fprintf('Initial guess %d, alpha %.8g: RRE %.6g, PSF error %.6g, stop %d\n', ...
            start,alpha_grid(j),RRE(j,start),PSF_error(j,start),stop(j,start));
    end
end

%% Select the weight with the smallest mean image RRE and save the results
mean_RRE = mean(RRE,2);
[~,index] = min(mean_RRE);
selected_alpha = alpha_grid(index);
fprintf('Selected alpha %.12g (mean image RRE across all three initial guesses).\n',selected_alpha);
save(fullfile(output_folder,'results.mat'),'alpha_grid','selected_alpha', ...
    'RRE','PSF_error','parameter_error','y_estimate','iterations','stop','history', ...
    'x_true','b','bhat','psf_true','y_true','y_start','eta','delta','lambda','gamma','s');
T = table(alpha_grid(:),RRE(:,1),RRE(:,2),RRE(:,3),mean_RRE, ...
    'VariableNames',{'Alpha','RRE_initial_guess_1','RRE_initial_guess_2','RRE_initial_guess_3','Mean_RRE'});
disp(T);
writetable(T,fullfile(output_folder,'selection.csv'));

%% Display the three image-error curves and the selected family weight
[fig,figure_cleanup]=output_figure('Visible','off','Color','w');
semilogx(alpha_grid,RRE,'-o');
hold on;
xline(selected_alpha,'--k','Selected alpha','HandleVisibility','off');
xlabel('$\alpha$','Interpreter','latex');
ylabel('Image RRE');
legend('Initial guess 1','Initial guess 2','Initial guess 3','Location','best');
grid on;
exportgraphics(fig,fullfile(output_folder,'alpha_selection.png'),'Resolution',200);
end
