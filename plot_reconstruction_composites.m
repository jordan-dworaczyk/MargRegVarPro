function plot_reconstruction_composites
% Paper figures from the saved reconstruction results, including every initial guess.

experiment_folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(experiment_folder,'output','reconstruction_figures');
if ~exist(output_folder,'dir'), mkdir(output_folder); end
types={'gaussian','defocus','motion'};
method_order=[1 3 2]; % No log-det, Gaussian marginalization, Proposed
method_labels={'No log-det','Gaussian marginalization','Proposed'};
all_rows=cell(0,13);

for type=1:3
    a=load(fullfile(experiment_folder,'output',['reconstruction_' types{type}],'results.mat'));
    n=size(a.x_true,1);
    zoom_cols=round(.25*n):round(.25*n)+round(.34*n);
    zoom_rows=round(.20*n):round(.20*n)+round(.34*n);
    halfwidth=[14 10 15];
    psf_crop=n/2+1+(-halfwidth(type):halfwidth(type));
    [fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.55 4.9]);
    tiles=tiledlayout(fig,3,3,'Padding','compact','TileSpacing','compact');
    for row=1:3
        method=method_order(row);
        for start=1:3
            image_crop=a.x_estimate(zoom_rows,zoom_cols,method,start);
            image_crop=max(0,min(1,image_crop)); % display only; metrics use raw x
            psf=a.psf_estimate(psf_crop,psf_crop,method,start);
            psf=max(0,psf); psf=psf/max(psf(:)); % independent peak normalization
            ir=round(linspace(1,size(image_crop,1),84));
            ic=round(linspace(1,size(image_crop,2),84));
            pr=round(linspace(1,size(psf,1),84));
            pc=round(linspace(1,size(psf,2),84));
            panel=[image_crop(ir,ic),.5*ones(84,4),psf(pr,pc)];
            ax=nexttile(tiles); imagesc(ax,panel,[0 1]); axis(ax,'image'); axis(ax,'off');
            ax.Toolbar.Visible='off';
            colormap(ax,gray);
            if row==1
                initial=sprintf('%.3g,',a.y_start(:,start)); initial=initial(1:end-1);
                if numel(a.y_true)>1, initial=['(' initial ')']; end
                title(ax,sprintf('Initial Guess %d: $\\mathbf y_0=%s$',start,initial), ...
                    'FontName','Times','FontSize',8.5,'FontWeight','normal','Interpreter','latex');
            end
            if start==1
                text(ax,-.18,.5,method_labels{row},'Units','normalized','Rotation',90, ...
                    'HorizontalAlignment','center','VerticalAlignment','middle', ...
                    'Interpreter','latex','FontSize',8.2,'Clipping','off');
            end
            all_rows(end+1,:)={types{type},method_labels{row},start,mat2str(a.y_estimate(:,method,start).',12), ...
                a.parameter_error(method,start),a.relative_residual(method,start), ...
                a.RRE(method,start),a.PSNR(method,start),a.SSIM(method,start), ...
                a.PSF_error(method,start),a.iterations(method,start),a.stop(method,start),a.alpha(method)};
        end
    end
    name=fullfile(output_folder,['paper_fig_composite_recon_psf_' types{type}]);
    % Keep the full figure dimensions when vector export trims margins.
    annotation(fig,'rectangle',[0 0 1 1],'Color','w','LineWidth',0.1);
    drawnow;
    % Vector export preserves text and avoids hidden-figure raster clipping in R2024b.
    exportgraphics(fig,[name '.pdf'],'ContentType','vector','Resolution',450);
    exportgraphics(fig,[name '.png'],'Resolution',450);
    close(fig);
end
metrics=cell2table(all_rows,'VariableNames',{'PSF','Method','Initial_guess','Estimate', ...
    'Parameter_error','Relative_residual','RRE','PSNR','SSIM','PSF_error','Iterations','Stop','Alpha'});
writetable(metrics,fullfile(output_folder,'all_initial_guess_metrics.csv'));
disp(metrics);
end
