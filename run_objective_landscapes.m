function run_objective_landscapes
% Gaussian-blur objective landscapes, with and without the determinant prior.
% Generate all six panels. Requires Image Processing Toolbox.
previous_rng=rng;
rng_cleanup=onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(0,'twister');

%% Experiment parameters
n = 128;
N = n^2;
sigma_true = [3, 5];
sigma_prior = [4, 4];
prior_precision = [0, 100, 1000]; % uniform, weaker Gaussian, stronger Gaussian
noise_level = 0.01;
delta = 100;                    % image-prior precision
sigma_grid = linspace(0.01, 7, 80);
experiment_folder = fileparts(mfilename('fullpath'));
output_folder = fullfile(experiment_folder, 'output', 'landscapes');
if ~exist(output_folder, 'dir'), mkdir(output_folder); end

%% True image, Gaussian PSF, and noisy observations
x_true = imresize(im2double(imread(fullfile(experiment_folder, ...
    'cameraman.png'))), [n, n], 'bilinear');
[u, v] = ndgrid(-n/2:n/2-1);
% Sampled Gaussian parameterization: exp(-u^2/sigma1^2-v^2/sigma2^2).
psf = exp(-u.^2/sigma_true(1)^2 - v.^2/sigma_true(2)^2);
psf(psf < eps) = 0;
psf = psf / sum(psf(:));
c_true = fft2(circshift(psf, [-n/2, -n/2]));
b_true = real(ifft2(c_true .* fft2(x_true)));

noise = randn(n);
noise = noise / norm(noise, 'fro');
b = b_true + noise_level * norm(b_true, 'fro') * noise;
bhat = fft2(b);
sigma_noise = noise_level * norm(b_true, 'fro') / sqrt(N);
eta = 1 / sigma_noise^2;
lambda = sqrt(delta / eta);

%% Image regularizer and matched auxiliary parameters
% FFT frequencies start at zero. |s|^2 is the sum of squared separate
% second-difference symbols, not the square of the 2D Laplacian symbol.
omega = 2*pi/n * (0:n-1);
d2 = 2*cos(omega) - 2;
s = sqrt(d2(:).^2 + d2.^2);
ell = lambda^2 * abs(s).^2;

% Matched auxiliary choice: tilde Q = Q, with alpha = 1.
alpha = 1;
beta = sqrt(alpha / eta);
gamma = delta / eta;

%% VarPro and determinant contributions at each blur parameter
ngrid = numel(sigma_grid);
varpro_values = zeros(ngrid);
logdet_values = zeros(ngrid);
for i = 1:ngrid
    sigma1 = sigma_grid(i);
    for j = 1:ngrid
        sigma2 = sigma_grid(j);
        psf = exp(-u.^2/sigma1^2 - v.^2/sigma2^2);
        psf(psf < eps) = 0;
        psf = psf / sum(psf(:));
        c = fft2(circshift(psf, [-n/2, -n/2]));
        a = abs(c).^2;
        q = a + gamma * abs(s).^2;

        % Residual units: unnormalized fft2 requires the Parseval factor 1/N.
        varpro = sum(ell(:) .* abs(bhat(:)).^2 ./ ...
                     (a(:) + ell(:))) / (2*N);
        logdet = beta^2/2 * sum(log(q(:)));

        % Store the objective contributions in precision units.
        % log det Q = N log eta + sum log q; N log eta is constant in y.
        varpro_values(j,i) = eta * varpro;
        logdet_values(j,i) = eta * logdet + alpha*N/2 * log(eta);
    end
end

%% Add the base-prior penalties explicitly
[sigma1, sigma2] = meshgrid(sigma_grid, sigma_grid);
distance_squared = (sigma1-sigma_prior(1)).^2 + (sigma2-sigma_prior(2)).^2;
% These penalties are -log p0 in precision units, omitting constants.
weaker_prior = prior_precision(2)/2 * distance_squared;
stronger_prior = prior_precision(3)/2 * distance_squared;

map_uniform = varpro_values;
map_weaker = varpro_values + weaker_prior;
map_stronger = varpro_values + stronger_prior;
logdet_uniform = varpro_values + logdet_values;
logdet_weaker = varpro_values + logdet_values + weaker_prior;
logdet_stronger = varpro_values + logdet_values + stronger_prior;

%% Grid minimizers and saved numerical results
% Rows: ordinary MAP, determinant MAP. Columns: the three base priors.
landscapes = zeros(ngrid, ngrid, 2, 3);
landscapes(:,:,1,1) = map_uniform;
landscapes(:,:,1,2) = map_weaker;
landscapes(:,:,1,3) = map_stronger;
landscapes(:,:,2,1) = logdet_uniform;
landscapes(:,:,2,2) = logdet_weaker;
landscapes(:,:,2,3) = logdet_stronger;
minimizers = zeros(2, 2, 3); % parameter, objective row, prior column
for column = 1:3
    for row = 1:2
        values = landscapes(:,:,row,column);
        [~, index] = min(values(:));
        [j, i] = ind2sub(size(values), index);
        minimizers(:,row,column) = [sigma_grid(i); sigma_grid(j)];
        fprintf('Prior precision %g, log-det %d: (%.4f, %.4f)\n', ...
            prior_precision(column), row-1, sigma_grid(i), sigma_grid(j));
    end
end
save(fullfile(output_folder, 'landscape_results.mat'), ...
    'landscapes', 'minimizers', 'sigma_grid', 'sigma_true', 'sigma_prior', ...
    'prior_precision', 'eta', 'delta', 'lambda', 'alpha', 'beta', 'gamma');

%% Plot the six landscapes in the manuscript residual units
% Convert the saved precision-unit arrays to residual units for display.
display_landscapes = landscapes/eta;
display_landscapes(:,:,2,:) = display_landscapes(:,:,2,:) - alpha*N/(2*eta)*log(eta);
% The affine transformation preserves grid minimizers.
titles = {'Uniform prior', 'Gaussian prior, $w_{hp}=10^2/\eta$', ...
          'Gaussian prior, $w_{hp}=10^3/\eta$'};
row_labels = {'No log-det', 'With log-det'};
[fig,figure_cleanup]=output_figure('Visible','off','Color','w', 'Units','inches', 'Position',[1 1 15.2 8]);
tiles = tiledlayout(fig, 2, 3, 'Padding','compact', 'TileSpacing','compact');
for row = 1:2
    row_values = display_landscapes(:,:,row,:);
    limits = [min(row_values(:)), max(row_values(:))];
    filled_levels = linspace(limits(1), limits(2), 24);
    line_levels = linspace(limits(1), limits(2), 10);
    for column = 1:3
        ax = nexttile(tiles, (row-1)*3 + column);
        values = display_landscapes(:,:,row,column);
        contourf(ax,sigma_grid,sigma_grid,values,filled_levels,'LineStyle','none');
        hold(ax,'on');
        contour(ax,sigma_grid,sigma_grid,values,line_levels, ...
            'LineColor',[.2 .2 .2], 'LineWidth',.7);
        h_min = plot(ax,minimizers(1,row,column),minimizers(2,row,column), ...
            'ko','MarkerSize',8,'LineWidth',1.4,'MarkerFaceColor','w');
        h_true = plot(ax,sigma_true(1),sigma_true(2),'kp', ...
            'MarkerSize',13,'LineWidth',1.2,'MarkerFaceColor',[1 .85 .15]);
        if column > 1
            % Cross polygon: white fill with a thin black outline.
            cross_x = [-1 -.72 0 .72 1 .28 1 .72 0 -.72 -1 -.28];
            cross_y = [-.72 -1 -.28 -1 -.72 0 .72 1 .28 1 .72 0];
            cross_scale = 0.15;
            h_prior = patch(ax,sigma_prior(1)+cross_scale*cross_x, ...
                sigma_prior(2)+cross_scale*cross_y,'w', ...
                'EdgeColor','k','LineWidth',0.7);
        end
        axis(ax,'square'); xlim(ax,[0 7]); ylim(ax,[0 7]);
        colormap(ax,parula(256)); clim(ax,limits);
        set(ax,'FontSize',18,'LineWidth',1,'TickDir','out','Layer','top', ...
            'Box','on','TickLabelInterpreter','latex');
        xlabel(ax,'$\sigma_1$','Interpreter','latex','FontSize',18);
        ylabel(ax,'$\sigma_2$','Interpreter','latex','FontSize',18);
        if row == 1
            title(ax,titles{column},'Interpreter','latex','FontWeight','normal','FontSize',18);
        end
        if column == 1
            text(ax,-.30,.50,row_labels{row},'Units','normalized','Rotation',90, ...
                'HorizontalAlignment','center','VerticalAlignment','middle', ...
                'Interpreter','latex','FontSize',18,'Clipping','off');
        end
        if column == 3
            cb = colorbar(ax);
            cb.TickLabelInterpreter = 'latex'; cb.FontSize = 18;
        end
    end
end
lgd = legend([h_min,h_true,h_prior], ...
    {'Minimizer','True parameters','Gaussian-prior mean'}, ...
    'Interpreter','latex','Orientation','horizontal','FontSize',18, ...
    'Location','southoutside','Box','off');
lgd.Layout.Tile = 'south';
% MATLAB replaces patch shapes with rectangles in its standard legend.
% Draw the key explicitly so it shows the same outlined cross as the panels.
lgd.Visible = 'off';
key = axes(fig,'Position',[0 0 1 .045],'Visible','off');
xlim(key,[0 15.2]); ylim(key,[0 .36]); hold(key,'on');
plot(key,3.8,.18,'ko','MarkerSize',8,'LineWidth',1.4,'MarkerFaceColor','w');
text(key,4.0,.18,'Minimizer','Interpreter','latex','FontSize',18);
plot(key,6.1,.18,'kp','MarkerSize',13,'LineWidth',1.2,'MarkerFaceColor',[1 .85 .15]);
text(key,6.3,.18,'True parameters','Interpreter','latex','FontSize',18);
patch(key,8.95+.07*cross_x,.18+.07*cross_y,'w','EdgeColor','k','LineWidth',.7);
text(key,9.15,.18,'Gaussian-prior mean','Interpreter','latex','FontSize',18);
figure_name = 'F_sigma2D_three_priors_raw_objectives_2x3_shared_row_colorbars_large_font';
exportgraphics(fig,fullfile(output_folder,[figure_name '.png']),'Resolution',400);
exportgraphics(fig,fullfile(output_folder,[figure_name '.pdf']),'ContentType','vector');
fprintf('Saved results and figures in %s\n',output_folder);
end
