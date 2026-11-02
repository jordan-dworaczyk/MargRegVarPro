function plot_gcv_single_sample
% Display the saved single-sample fits; no optimization is repeated here.

folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(folder,'output','gcv_single_sample');
if ~exist(output_folder,'dir'), mkdir(output_folder); end
data=load(fullfile(folder,'output','gcv_single_sample','results.mat'));
x_true=data.x_true;
x_estimate=data.x_estimate;
psf_true=data.psf_true;
psf_estimate=data.psf_estimate;
trace=data.trace;
alpha=data.alpha;
gcv_score=data.gcv_score;
lambda=data.lambda;
y_estimate=data.y_estimate;
RRE=data.RRE;
PSF_error=data.PSF_error;
SSIM=data.SSIM;
n=size(x_true,1);
% Match the reconstruction composites' normal-weight LaTeX headings (8.5 pt).
% Three paired groups: equally sized image and PSF squares with the same
% mid-gray separator proportion (4/84) as Figures 3--5.
figure_width=7.5; figure_height=1.48;
panel_side=1.10; separator_width=panel_side*4/84;
pair_width=2*panel_side+separator_width;
group_gap=0.20; panel_bottom=0.08;
left_margin=(figure_width-3*pair_width-2*group_gap)/2;
pair_left=left_margin+(0:2)*(pair_width+group_gap);
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches', ...
    'Position',[1 1 figure_width figure_height]);
set(fig,'defaultAxesFontName','Times','defaultTextFontName','Times');
images={x_true,x_estimate(:,:,1),x_estimate(:,:,2)};
% Display-only crop matching Figures 3--5; full-image metrics are unchanged.
zoom_rows=51:138; zoom_cols=64:151;
kernels={psf_true,psf_estimate(:,:,1),psf_estimate(:,:,2)};
names={'Truth','No log-det','Automatic selection'};
% Show 25-by-25 PSF crops with independent automatic color limits.
crop=n/2+1+(-12:12);
for j=1:3
    ax=axes(fig,'Units','inches','Position', ...
        [pair_left(j) panel_bottom panel_side panel_side]);
    imagesc(ax,images{j}(zoom_rows,zoom_cols),[0 1]); axis(ax,'image'); axis(ax,'off');
    colormap(ax,gray); ax.Toolbar.Visible='off';

    annotation(fig,'rectangle', ...
        [(pair_left(j)+panel_side)/figure_width panel_bottom/figure_height ...
         separator_width/figure_width panel_side/figure_height], ...
        'FaceColor',[.5 .5 .5],'EdgeColor','none');
    ax=axes(fig,'Units','inches','Position', ...
        [pair_left(j)+panel_side+separator_width panel_bottom panel_side panel_side]);
    imagesc(ax,kernels{j}(crop,crop)); axis(ax,'image'); axis(ax,'off');
    colormap(ax,gray); ax.Toolbar.Visible='off';

    heading_ax=axes(fig,'Units','inches','Position', ...
        [pair_left(j) panel_bottom pair_width panel_side],'Color','none');
    axis(heading_ax,'off'); heading_ax.Toolbar.Visible='off';
    title(heading_ax,names{j},'FontName','Times','FontSize',8.5, ...
        'FontWeight','normal','Interpreter','latex');
end
% A white boundary retains the small outer margins when R2024b tightly crops
% the export. It is invisible on the white background and adds no data panel.
annotation(fig,'rectangle',[0 0 1 1],'Color','w','LineWidth',0.1);
% Retain embedded vector fonts (print/painters outlines LaTeX text in R2024b).
exportgraphics(fig,fullfile(output_folder,'gcv_single_sample_reconstruction.pdf'), ...
    'ContentType','vector','Resolution',450);
exportgraphics(fig,fullfile(output_folder,'gcv_single_sample_reconstruction.png'),'Resolution',250);

[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Position',[100 100 650 400]);
set(fig,'defaultAxesFontName','Times','defaultTextFontName','Times');
[sorted_alpha,order]=sort(trace(:,1));
plot(sorted_alpha,trace(order,2),'-o'); hold on;
plot(alpha(2),gcv_score(2),'s','MarkerFaceColor',[.85 .2 .15]);
xlabel('$\alpha$','Interpreter','latex'); ylabel('Final conditional GCV'); grid on;
title('Evaluations from the bounded search');
exportgraphics(fig,fullfile(output_folder,'gcv_search_trace.png'),'Resolution',200);

fid=fopen(fullfile(output_folder,'table_rows.tex'),'w');
names={'No log-det','Automatic selection'};
for j=1:2
    fprintf(fid,'%s & %.3f & $%.3g$ & $(%.3f,%.3f)$ & %.3f & %.3f & %.3f \\\\ \n', ...
        names{j},alpha(j),lambda(j),y_estimate(:,j),RRE(j),PSF_error(j),SSIM(j));
end
fclose(fid);
end
