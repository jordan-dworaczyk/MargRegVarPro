function plot_reconstruction_landscapes
% Objective landscapes and all optimization trajectories for the saved fits.
% Use the paper's grid sizes, range rules, and N*f display scaling.

folder=fileparts(mfilename('fullpath'));
output_folder=fullfile(folder,'output','reconstruction_figures');
if ~exist(output_folder,'dir'), mkdir(output_folder); end
types={'gaussian','motion','defocus'};
row_labels={'Gaussian','Motion','Defocus'};
method_order=[1 3 2];
method_labels={'No log-det','Gaussian marginalization','Proposed'};
colors=[0 .45 .74;.20 .63 .17;.80 .20 .10];
[fig,figure_cleanup]=output_figure('Visible','off','Color','w','Units','inches','Position',[1 1 7.05 6.35]);
tiles=tiledlayout(fig,3,3,'Padding','compact','TileSpacing','compact');
tiles.Units='normalized'; tiles.OuterPosition=[.065 .095 .91 .87];
% Color anchors for the objective-landscape palette.
anchors=[.10 .23 .36;.12 .49 .61;.62 .77 .62;.95 .83 .48;.98 .94 .68];
palette=interp1(linspace(0,1,5),anchors,linspace(0,1,256),'pchip');
colormap(fig,max(0,min(1,palette)));
for row=1:3
    a=load(fullfile(folder,'output',['reconstruction_' types{row}],'results.mat'));
    n=size(a.x_true,1); N=n^2;
    [r,c]=ndgrid(-n/2:n/2-1);
    ell=a.lambda^2*abs(a.s).^2;
    for column=1:3
        method=method_order(column);
        ax=nexttile(tiles,(row-1)*3+column);
        all_y=[a.y_true,a.y_start,reshape(a.y_estimate(:,method,:),numel(a.y_true),3)];
        if row==3
            xlimits=[max(.001,min(all_y)-.5),max(all_y)+.8];
            xgrid=linspace(xlimits(1),xlimits(2),900);
            XX=xgrid; YY=zeros(size(xgrid));
        else
            if row==2
                xlimits=[.001,max(all_y(1,:))+2.4]; ylimits=[-pi/6,pi/6]; ng=30;
            else
                xlimits=[max(.25,min(all_y(1,:))-1.5),max(all_y(1,:))+1.5];
                ylimits=[max(.25,min(all_y(2,:))-1.5),max(all_y(2,:))+1.5]; ng=32;
            end
            xgrid=linspace(xlimits(1),xlimits(2),ng);
            ygrid=linspace(ylimits(1),ylimits(2),ng);
            [XX,YY]=meshgrid(xgrid,ygrid);
        end
        if method==3, s_aux=a.s; else, s_aux=ones(n); end
        F=zeros(size(XX));
        for point=1:numel(XX)
            if row==1
                g=exp(-r.^2/XX(point)^2-c.^2/YY(point)^2);
            elseif row==2
                u=cos(YY(point))*c-sin(YY(point))*r;
                v=sin(YY(point))*c+cos(YY(point))*r;
                g=exp(-v.^2/(2*.30^2)).*.5.* ...
                    (tanh((u+XX(point)/2)/.15)-tanh((u-XX(point)/2)/.15));
            else
                g=1./(1+exp((sqrt(r.^2+c.^2)-XX(point))/.5));
            end
            psf=g/sum(g(:));
            symbol=fft2(circshift(psf,[-n/2,-n/2]));
            spectrum=abs(symbol).^2;
            q=spectrum+a.gamma(method)*abs(s_aux).^2;
            varpro=sum(ell(:).*abs(a.bhat(:)).^2./(spectrum(:)+ell(:)))/(2*N);
            logdet=a.beta(method)^2/2*sum(log(q(:)));
            F(point)=N*(varpro+logdet);
        end
        Z=log10(1+max(F-min(F(:)),0));
        if row==3
            plot(ax,xgrid,Z,'Color',[.18 .18 .18],'LineWidth',1.15); hold(ax,'on');
            ylim(ax,[0 max(.5,prctile(Z,96))]); yticks(ax,[]);
            truth_z=interp1(xgrid,Z,a.y_true,'linear','extrap');
            plot(ax,a.y_true,truth_z,'kp','MarkerFaceColor',[1 .9 0],'MarkerSize',7.5);
        else
            imagesc(ax,xgrid,ygrid,Z); set(ax,'YDir','normal'); hold(ax,'on');
            contour(ax,XX,YY,Z,8,'Color',[.15 .15 .15],'LineWidth',.28);
            clim(ax,[0 prctile(Z(:),94)]); ylim(ax,ylimits);
            plot(ax,a.y_true(1),a.y_true(2),'kp','MarkerFaceColor',[1 .9 0],'MarkerSize',7.5);
        end
        for start=1:3
            trajectory=[a.history{method,start},a.y_estimate(:,method,start)];
            initial=a.y_start(:,start); endpoint=a.y_estimate(:,method,start);
            if row==2
                trajectory(2,:)=atan2(sin(trajectory(2,:)),cos(trajectory(2,:)));
                initial(2)=atan2(sin(initial(2)),cos(initial(2)));
                endpoint(2)=atan2(sin(endpoint(2)),cos(endpoint(2)));
            end
            if row==3
                trajectory(2,:)=interp1(xgrid,Z,trajectory(1,:),'linear','extrap');
                initial(2)=interp1(xgrid,Z,initial(1),'linear','extrap');
                endpoint(2)=interp1(xgrid,Z,endpoint(1),'linear','extrap');
                outside=trajectory(1,:)<xlimits(1) | trajectory(1,:)>xlimits(2);
            else
                outside=trajectory(1,:)<xlimits(1) | trajectory(1,:)>xlimits(2) | ...
                    trajectory(2,:)<ylimits(1) | trajectory(2,:)>ylimits(2);
            end
            trajectory(:,outside)=nan;
            plot(ax,trajectory(1,:),trajectory(2,:),'.-','Color',colors(start,:), ...
                'LineWidth',.72,'MarkerSize',3.5);
            plot(ax,initial(1),initial(2),'o','Color',colors(start,:), ...
                'MarkerFaceColor','w','MarkerSize',4.2);
            if a.stop(method,start)~=0
                plot(ax,endpoint(1),endpoint(2),'ks','MarkerFaceColor',colors(start,:),'MarkerSize',5.2);
            else
                plot(ax,endpoint(1),endpoint(2),'ks','MarkerFaceColor','w','MarkerSize',7);
                plot(ax,endpoint(1),endpoint(2),'kx','MarkerSize',7.2);
            end
        end
        xlim(ax,xlimits);
        set(ax,'FontName','Times','FontSize',7.4,'LineWidth',.65,'TickDir','out', ...
            'Box','off','TickLabelInterpreter','latex');
        if row==1, title(ax,method_labels{column},'Interpreter','latex','FontWeight','normal'); end
        if column==1
            if row==1, ylabel(ax,{row_labels{row},'$\sigma_2$'},'Interpreter','latex');
            elseif row==2, ylabel(ax,{row_labels{row},'$\theta$'},'Interpreter','latex');
            else, ylabel(ax,{row_labels{row},'$\log_{10}(1+\Delta F)$'},'Interpreter','latex'); end
        end
        if row==3, xlabel(ax,'$R$','Interpreter','latex'); end
        if row==2
            xticks(ax,[5 10 15]); yticks(ax,[-pi/6 0 pi/6]);
            yticklabels(ax,{'$-\pi/6$','0','$\pi/6$'});
        end
        label=sprintf('(%c)',char('a'+(row-1)*3+column-1));
        text(ax,.02,.98,label,'Units','normalized','VerticalAlignment','top', ...
            'FontName','Times','FontSize',7.4,'BackgroundColor','w');
    end
end
key=axes(fig,'Position',[.07 .015 .86 .025],'Visible','off'); hold(key,'on');
xlim(key,[0 1]); ylim(key,[0 1]);
plot(key,.03,.5,'kp','MarkerFaceColor',[1 .9 0],'MarkerSize',7);
text(key,.055,.5,'Truth','FontSize',8);
plot(key,.21,.5,'ko','MarkerFaceColor','w','MarkerSize',4.2);
text(key,.235,.5,'Initial Guess','FontSize',8);
plot(key,.46,.5,'ks','MarkerFaceColor',colors(1,:),'MarkerSize',5.2);
text(key,.485,.5,'Stopped run','FontSize',8);
plot(key,.73,.5,'ks','MarkerFaceColor','w','MarkerSize',7);
plot(key,.73,.5,'kx','MarkerSize',7.2);
text(key,.76,.5,'Iteration limit','FontSize',8);
name=fullfile(output_folder,'paper_fig_landscape_grid_siam');
annotation(fig,'rectangle',[0 0 1 1],'Color','w','LineWidth',0.1);
drawnow;
% Vector export avoids hidden-figure raster clipping in R2024b.
exportgraphics(fig,[name '.pdf'],'ContentType','vector','Resolution',450);
exportgraphics(fig,[name '.png'],'Resolution',450);
end
