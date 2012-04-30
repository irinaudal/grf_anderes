function   plot_mcmcresult(draws, id, nf, z_grid_Pz, Pz, zS, x_obs,z_obs,s_obs, f_true, x_gridM,z_gridM,res_dir,add_obs_values,plot_type,rate)


if(nargin < 16)
    rate = 1;
end
if(nargin < 15)
    plot_type = 2;
end
if(nargin < 14)
    add_obs_values = 0;
end
if(add_obs_values==0)
    add_obs_values_type = 0;
    %add_obs_values_type = 1;  %plot z_spectroscopic
    %add_obs_values_type = 2;  %plot z_mcmc
end




N = numel(zS);
niter = size(draws,1);
[m2 m1] = size(x_gridM);
draws_z = draws(:,1:N);
draws_f = draws(:,(N+1):end);

Pz = Pz(id,:);
dz = z_grid_Pz(2)-z_grid_Pz(1);
zS = zS(id);

plotN = 0;

if(plot_type == 1)
    fig = figure();
    disp('Printing Trace plot')
    if(length(id) > 4)
        id = 1:4;
    end
    for j=1:length(id)

            plotN = plotN+1;
            subplot(length(id),2,plotN);

            plot(1:niter,draws_z(:,j));
            xlabel('Iterations');  title(['Trace of ','y',num2str(id(j))]);

            plotN = plotN+1;
            subplot(length(id),2,plotN);

            [hi zi] = hist(draws_z(:,j));
            norm_hist = hi./sum(hi)/(zi(2)-zi(1));
            bar(zi, norm_hist, 1,'FaceColor',[0.7 0.7 0.9]);

    %        hist(draws_y(:,j));
    %        h = findobj(gca,'Type','patch');
    %        set(h,'FaceColor',[0.7 0.7 0.9]);

            xlabel(['y',num2str(id(j))]);  title(['Posterior of ','y',num2str(id(j))]);
            hold on;

            [fi xi] = ksdensity(draws_z(:,j));
            plot(xi,fi,'b','LineWidth',2);

            plot(z_grid_Pz, Pz(j,:)/dz,'r','LineWidth',2);

            plot(zS(j),0,'o','MarkerFaceColor','g')

            id1 = Pz(j,:)>10^-4;
            z_range1 = z_grid_Pz( id1 )';
            id2 = fi > 10^-4;
            z_range2 = xi( id2 );
            axis([min([z_range1 z_range2]) max([z_range1 z_range2])  0 max([max(Pz(j,:))/dz, max(fi), max(norm_hist)])]);

            legend('Posterior Hist','Posterior Smooth','Prior',['Truth y=',num2str(zS(j))] ,'Location','NorthEast');
            hold off;

    %         subplot(2,3,3);
    %         plot(z_grid_Pz, Pz_z(1,:));
    %         xlabel('z_1'); ylabel('density'); title('Prior  p(z)');
    %         hold on;
    %         plot(zS(1),0,'o','MarkerFaceColor','g')
    %         hold off;
    %         subplot(2,3,4);
    %         plot(1:niter,first_few_draws_z(:,2));
    %         xlabel('iter'); ylabel('z_2');
    %         subplot(2,3,5);
    %         hist(first_few_draws_z(:,2))
    %         xlabel('z_2'); ylabel('frequency');
    %         hold on;
    %         plot(zS(2),0,'o','MarkerFaceColor','g')
    %         hold off;
    %         subplot(2,3,6);
    %         plot(z_grid_Pz, Pz(2,:));
    %         xlabel('z_2'); ylabel('density');
    %         hold on;
    %         plot(zS(2),0,'o','MarkerFaceColor','g')
    %         hold off;
    %print(fig,'-depsc','-tiff','-r300',['result',num2str(niter),'_',num2str(N),'_',num2str(m1),'/fig_mcmc_y.eps']);

    end
    close(fig);
else

    %find number of subplots to make
    v = [1 3 5 7 8 11 15 19];
    if (nf > v(end))
        id = v(end)+1; 
        nrow = 4;
        ncol = 5;
    elseif (nf == 1)
        id = nf+1;
        nrow = 1;
        ncol = 2;
    else
        for k = 2:length(v)
            if(v(k-1)<=nf && nf<v(k))
                id = v(k-1)+1;
                if(id <= 8)
                    nrow = 2;
                    ncol = id/2;
                elseif(id <= 12)
                    nrow = 3;
                    ncol = id/3;
                else
                    nrow = 4;
                    ncol = id/4;
                end
            end
        end
    end
    %disp([id nrow ncol])

    %%%%%% Plot GRF overlayed with true observed locations
    if(add_obs_values == 1)
        add_obs_values_type=1;
    end 
    
    %plot truth
    fig = figure();
    subplot(nrow,ncol,1);    
    pcolor(x_gridM,z_gridM,f_true);   %surf() plot viewed from above
    shading flat;      
    xlabel('x'); ylabel('z'); zlabel('f(x,z)');  
    title('True GRF, Equidist Grid, no noise');
    h2 = caxis;
    colorbar('East');

    if(add_obs_values > 0)
        hold on;
        scatter(x_obs,z_obs,12,s_obs,'filled'); %plot: a=vector of location, b=vector of redshift, size, color choice, 
        scatter(x_obs,z_obs,15,'o','MarkerEdgeColor',[.5 .5 .5]);        
        title('True GRF, Equidist Grid & z_s, no noise.'); 
        hold off;
    end
    
    for j=2:id
        %plot as many iterations of simulated grf as possible
            %starting from end, plot every k-th, where k=rate
            iter = niter-id*rate+(j)*rate;
        disp(iter);
        subplot(nrow,ncol,j);           
        pcolor(x_gridM,z_gridM,reshape(draws_f(iter,:),m2,m1));   %surf() plot viewed from above
        shading flat;      
        xlabel('x'); ylabel('z'); zlabel('f(x,z)');  
        title(['MCMC draw of f.  Iter = ',num2str(iter)]); 
        caxis(h2);

        if(add_obs_values_type == 1)
            hold on;
            scatter(x_obs,z_obs,12,s_obs,'filled'); %plot: a=vector of location, b=vector of redshift, size, color choice, 
            scatter(x_obs,z_obs,15,'o','MarkerEdgeColor',[.5 .5 .5]);
            title(['MCMC draw of f & true z_s.  Iter = ',num2str(iter)]); 
            hold off;
        end

        %starting from end, plot every 10th    
    end
    colorbar('East');

    fullscreen = get(0,'ScreenSize');
    set(fig, 'Position', [0 0 fullscreen(3) fullscreen(4) ] );
    
    if(add_obs_values < 2)
        export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/fig_result.pdf'])
    end
%     if(add_obs_values == 2)
%         export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/fig_resultPhotZ.pdf'])
%     end
    %print(fig,'-depsc','-tiff','-r300',['result',num2str(niter),'_',num2str(N),'_',num2str(m1),'/fig_mcmc_f.eps']);

    
    if(add_obs_values > 0) % Do additional similar plot
        %%%%%% Plot GRF overlayed with MCMC sampled locations
        if(add_obs_values == 1)
            add_obs_values_type=2;
        end    
        %plot truth
        fig = figure();
        subplot(nrow,ncol,1);    
        pcolor(x_gridM,z_gridM,f_true);   %surf() plot viewed from above
        shading flat;      
        xlabel('x'); ylabel('z'); zlabel('f(x,z)');  
        title('True GRF, Equidist Grid, no noise');
        colorbar('East');
        h2 = caxis;

        if(add_obs_values > 0)
            hold on;
            scatter(x_obs,z_obs,12,s_obs,'filled'); %plot: a=vector of location, b=vector of redshift, size, color choice, 
            scatter(x_obs,z_obs,15,'o','MarkerEdgeColor',[.5 .5 .5]);        
            title('True GRF, Equidist Grid & z_s, no noise.'); 
            hold off;
        end

        for j=2:id
            %plot as many iterations of simulated grf as possible
                %starting from end, plot every k-th, where k=rate
                iter = niter-id*rate+(j)*rate;
            disp(iter);
            subplot(nrow,ncol,j);           
            pcolor(x_gridM,z_gridM,reshape(draws_f(iter,:),m2,m1));   %surf() plot viewed from above
            shading flat;      
            xlabel('x'); ylabel('z'); zlabel('f(x,z)');  
            title(['MCMC draw of f.  Iter = ',num2str(iter)]); 
            caxis(h2);

            if(add_obs_values_type == 2)
                hold on;
                scatter(x_obs,draws_z(iter,:),12,s_obs,'filled'); %plot: a=vector of location, b=vector of redshift, size, color choice, 
                scatter(x_obs,draws_z(iter,:),15,'o','MarkerEdgeColor',[.5 .5 .5]);
                title(['MCMC draws of f & z.  Iter = ',num2str(iter)]); 
                hold off;
            end

            %starting from end, plot every 10th    
        end
        colorbar('East');

        fullscreen = get(0,'ScreenSize');
        set(fig, 'Position', [0 0 fullscreen(3) fullscreen(4) ] );


        addpath('export_fig');
    %     if(add_obs_values < 2)
    %         export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/fig_result.pdf'])
    %     end
        if(add_obs_values_type == 2)
            export_fig(fig,'-painters', '-r300', '-q101', [res_dir,'/fig_resultPhotZ.pdf'])
        end
        %print(fig,'-depsc','-tiff','-r300',['result',num2str(niter),'_',num2str(N),'_',num2str(m1),'/fig_mcmc_f.eps']);
    end
end