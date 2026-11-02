function plot_blurred_data
% Plot the four saved blurred and noisy observations used in the paper.
folder=fileparts(mfilename('fullpath'));
outfolder=fullfile(folder,'output','reconstruction_figures');
if ~exist(outfolder,'dir'), mkdir(outfolder); end
studies={'reconstruction_gaussian','reconstruction_defocus', ...
    'reconstruction_motion','gcv_single_sample'};
labels={'Gaussian blur, 1\% noise','Defocus blur, 1\% noise', ...
    'Motion blur, 1\% noise','Gaussian blur, 5\% noise'};
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 1.85]);
tiles=tiledlayout(fig,1,4,'Padding','compact','TileSpacing','compact');
for k=1:4
    a=load(fullfile(folder,'output',studies{k},'results.mat'),'b');
    ax=nexttile(tiles); imagesc(ax,a.b,[0 1]);
    axis(ax,'image','off'); colormap(ax,gray); ax.Toolbar.Visible='off';
    title(ax,labels{k},'Interpreter','latex','FontWeight','normal','FontSize',8.5);
end
annotation(fig,'rectangle',[0 0 1 1],'Color','w','LineWidth',0.1);
exportgraphics(fig,fullfile(outfolder,'paper_fig_blurred_data.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outfolder,'paper_fig_blurred_data.png'),'Resolution',450);
close(fig);
end
