function plot_ground_truth_reference
% The paper's image and true PSFs, from the three reconstruction runs.
folder=fileparts(mfilename('fullpath'));
outfolder=fullfile(folder,'output','reconstruction_figures');
if ~exist(outfolder,'dir'), mkdir(outfolder); end
types={'gaussian','defocus','motion'}; titles={'Gaussian PSF','Defocus PSF','Motion PSF'};
a=load(fullfile(folder,'output','reconstruction_gaussian','results.mat'));
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.1 1.85]);
tl=tiledlayout(fig,1,4,'Padding','compact','TileSpacing','compact');
ax=nexttile(tl); imagesc(ax,a.x_true); axis(ax,'image','off'); colormap(ax,gray);
title(ax,'True image','Interpreter','latex','FontWeight','normal','FontSize',8.5);
for k=1:3
    a=load(fullfile(folder,'output',['reconstruction_' types{k}],'results.mat'));
    P=a.psf_true; w=41; first=floor((size(P,1)-w)/2)+1;
    ax=nexttile(tl); imagesc(ax,P(first:first+w-1,first:first+w-1));
    axis(ax,'image','off'); colormap(ax,gray);
    title(ax,titles{k},'Interpreter','latex','FontWeight','normal','FontSize',8.5);
end
exportgraphics(fig,fullfile(outfolder,'paper_fig_ground_truth_reference.pdf'),'ContentType','vector');
exportgraphics(fig,fullfile(outfolder,'paper_fig_ground_truth_reference.png'),'Resolution',450);
close(fig);
end
