function plot_auxiliary_scale_sensitivity
% Plot saved sweep results and evaluate the objective directly in Fourier space.

folder = fileparts(mfilename('fullpath'));
output_folder = fullfile(folder,'output','auxiliary_scale');
if ~exist(output_folder,'dir'), mkdir(output_folder); end
data=load(fullfile(folder,'output','auxiliary_scale','results.mat'));
parameter_error=data.parameter_error;
RRE=data.RRE;
PSF_error=data.PSF_error;
tilde_delta=data.tilde_delta;
gamma=tilde_delta/data.eta;
b=data.b;
lambda=data.lambda;
s=data.s;
eta=data.eta;
bhat=data.bhat;
alpha=data.alpha;
y_true=data.y_true;
y_estimate=data.y_estimate;
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.8 2.75]);
set(fig,'defaultAxesFontName','Times','defaultTextFontName','Times');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
metrics={parameter_error,RRE,PSF_error};
labels={'Parameter error','Reconstruction RRE','Relative PSF error'};
for j=1:3
    ax=nexttile;
    semilogx(gamma,metrics{j},'-o','MarkerSize',4,'LineWidth',1.2);
    if j~=2, ax.YScale='log'; end
    grid on; ax.GridAlpha=.12; ax.Toolbar.Visible='off';
    xlabel('$\gamma$','Interpreter','latex'); title(labels{j});
    xticks(10.^(-16:4:4)); xlim([gamma(1) gamma(end)]);
end
exportgraphics(fig,fullfile(output_folder,'paper_delta_sweep_smooth_only_main.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(output_folder,'paper_delta_sweep_smooth_only_main.png'),'Resolution',250);

%% Contributions in residual units: f, dropping the constant in log(eta)
n=size(b,1); N=n^2;
y_grid=unique([logspace(-12,log10(.75),200),linspace(.001,6.75,650)]);
selected=[1e-4 1 1e4];
varpro=zeros(size(y_grid)); logdet=zeros(3,numel(y_grid));
ell=lambda^2*s.^2;
for k=1:numel(y_grid)
    c=reshape(defocus_psf_symbol_derivs(y_grid(k),n),n,n);
    a=abs(c).^2;
    varpro(k)=1/(2*N)*sum(ell(:).*abs(bhat(:)).^2./(a(:)+ell(:)));
    for j=1:3
        logdet(j,k)=(alpha/eta)/2*sum(log(a(:)+selected(j)/eta));
    end
end
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.9 3]);
set(fig,'defaultAxesFontName','Times','defaultTextFontName','Times');
tl=tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
for j=1:3
    ax=nexttile; hold on;
    h=plot(y_grid,varpro,'--',y_grid,logdet(j,:),'-.', ...
        y_grid,varpro+logdet(j,:),'-','LineWidth',1.2);
    xline(y_true,':k','HandleVisibility','off');
    idx=find(tilde_delta==selected(j));
    c=reshape(defocus_psf_symbol_derivs(y_estimate(idx),n),n,n); a=abs(c).^2;
    f=1/(2*N)*sum(ell(:).*abs(bhat(:)).^2./(a(:)+ell(:))) ...
        +(alpha/eta)/2*sum(log(a(:)+selected(j)/eta));
    plot(y_estimate(idx),f,'s','MarkerFaceColor',[.85 .2 .15],'MarkerEdgeColor','k');
    ax.Toolbar.Visible='off'; grid on; ax.GridAlpha=.12;
    xlabel('$y$','Interpreter','latex');
    title(sprintf('$\\gamma=10^{%d}/\\eta$',round(log10(selected(j)))),'Interpreter','latex');
    xlim([0 6.75]);
    if j==1, ylabel('Objective contribution'); end
end
lg=legend(h,{'$\varphi(\mathbf{y})$', ...
    '$(\beta^2/2)\mathcal{R}_{\gamma}(\mathbf{y})$', ...
    '$\mathbf{f}_{\beta,\gamma}(\mathbf{y})$'}, ...
    'Interpreter','latex','Orientation','horizontal');
lg.Layout.Tile='south';
exportgraphics(fig,fullfile(output_folder,'paper_delta_sweep_smooth_only_objective_decomposition.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(output_folder,'paper_delta_sweep_smooth_only_objective_decomposition.png'),'Resolution',250);
end
