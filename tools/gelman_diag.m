
% 
%     % %Gelman & Rubin Diagnostic
%     addpath('mcmcstat/');
%     res_dir = 'result_simple_10000_100_21';
% 
%     Nchains=3;
%     
%     MAINchain = cell(Nchains,1);
%     dchain    = cell(Nchains);
%     schain    = cell(Nchains);
%     
%     for jj=1:Nchains
%         %run Gibbs sampler
%         disp('Reading data result for MCMC test...');
%         [ resmc ] = load([res_dir,'/mcmcdraws.txt']);
%         
%         
%         % compute statistics
%         disp('Checking MCMC draws...');
%         z_true = z_obs'; 
%         f_true = f_gridM(:)';
% %    [delta delta_true d c P P_true xsi stat_foci foci_true stat_ex ex_true stats stats_true] = diagnose_grf(resmc, ...
%         [stats stats_true delta delta_true] = diagnose_grf(resmc, ...
%             z_true, f_true, x_gridM(1,:),z_gridM(:,1), res_dir);
%         
%         %store all chains
%         MAINchain{jj} = resmc;
%         dchain{jj}    = delta;
%         schain{jj}    = stats;
% 
%     end
%     
%     %Make Diagnostic plots (script)
%     gelman_diag  
%  

% 
% temp = load([res_dir,'/mcmcdraws.txt']);
% chain{1} = temp(1:iter_plot,:);
% temp = load([res_dir,'/mcmcdraws.txt']);
% chain{2} = temp(1:iter_plot,:);
% temp = load([res_dir,'/mcmcdraws.txt'])/2;
% chain{3} = temp(1:iter_plot,:);
% 
% temp=[];

 
%evaluate the diagnostics ahead of time
mychain = cell(3,1);

lengthDraws = size(rchain{1},1);
iterations=[2:(lengthDraws-1)];
lengthR = length(iterations);
saveR1 = zeros(lengthR,size(MAINchain{1},2) +1);  % (# iterations-2, # diagnostics)
saveR2 = zeros(lengthR,size(schain{1},2) +1);  
saveR3 = zeros(lengthR,size(pchain{1},2) +1);  

% saveR1=saveR1(1:end,:);
% saveR2=saveR2(1:end,:);
% saveR3=saveR3(1:end,:);

for loop=1:(lengthR)
         
    ii = iterations(loop);

    % grf and 
    mychain{1} = MAINchain{1}(1:ii,:);
    mychain{2} = MAINchain{2}(1:ii,:);
    mychain{3} = MAINchain{3}(1:ii,:);
    %Gelman & Rubin Diagnostic
    R1 = psrf_gelrub(mychain); %need all 3 chains
%    R1 = psrf2(mychain);
%    R1 = psrf(mychain{1}, mychain{2}, mychain{3}); %need all 3 chains

    mychain{1} = schain{1}(1:ii,:);    %Stats
    mychain{2} = schain{2}(1:ii,:);
    mychain{3} = schain{3}(1:ii,:);
    R2 = psrf_gelrub(mychain); %need all 3 chains
    
    mychain{1} = pchain{1}(1:ii,:);    %Power Spectrum
    mychain{2} = pchain{2}(1:ii,:);
    mychain{3} = pchain{3}(1:ii,:);
    R3 = psrf_gelrub(mychain); %need all 3 chains

    
    saveR1(loop,:) = [ii R1];
    saveR2(loop,:) = [ii R2];
    saveR3(loop,:) = [ii R3];
end   


%%%%%%% PSRF plot %%%%%%%%%      %need all 3 chains   
    DD=20;  %starting point for PLOTS
    fig = figure();
       subplot(2,2,1);
    R_toplot = saveR1(DD:end,1+((N+1):(N+M)));
    plot(DD:lengthR, R_toplot); 
    axis([-1 lengthR 0.95 1.12]); 
    hold on;  hline = refline([0 1]);set(hline,'Color','k');  hold off; 
    xlabel('Iteration, t');
    ylabel('R');
    title('Potential Scale Reduction Factor of all Draws of f');
       subplot(2,2,3);
    if (fixed_z==false)
        R_toplot = saveR1(DD:end,1+(1:N)); %1st column are IDs
        plot(DD:lengthR, R_toplot); 
        axis([-1 lengthR 0.95 1.12]); 
        hold on;  hline = refline([0 1]);set(hline,'Color','k');  hold off; 
        xlabel('Iteration, t');
        ylabel('R');
        title('PSRF of all Draws of z');
    else
        plot(0);
        text(0.1,0,['Fixed redshift at: ',fixed_z],'BackgroundColor',[.98 .98 .98], 'FontSize',14, 'FontName', 'Arial');
        axis off;
    end
        subplot(2,2,2);
    R_toplot = saveR2(DD:end,((10-3+1):10)); %1st column are IDs
    plot(DD:lengthR, R_toplot); 
    axis([-1 lengthR 0.95 1.12]); 
    hold on;  hline = refline([0 1]);set(hline,'Color','k');  hold off; 
    xlabel('Iteration, t');
    ylabel('R');
    title('PSRF of Histograms');    
    
       subplot(2,2,4);
    R_toplot = saveR2(DD:end,2:(10-3)); %1st column are IDs
    plot(DD:lengthR, R_toplot); 
    axis([-1 lengthR 0.95 1.12]); 
    hold on;  hline = refline([0 1]);set(hline,'Color','k');  hold off; 
    xlabel('Iteration, t');
    ylabel('R');
    title('PSRF of statistics');
    
%         fullscreen = get(0,'ScreenSize');
%         set(fig, 'Position', [0 0 fullscreen(3) fullscreen(4) ] );
    print(fig,'-djpeg','-r300',[res_dir,'/psrf_ALL.jpg'])
%     export_fig(fig,'-painters', '-r300', '-q101', [res_dir,'/psrf.pdf'])
      
    
    fig=figure();
        R_toplot = saveR3(DD:end,2:end); %1st column are IDs
        hline = refline([0 1]);set(hline,'Color','k');
        axis([-1 lengthR 0.95 1.12]); 
        hold on;  plot(DD:lengthR, R_toplot);  hold off; 
    xlabel('Iteration, t');
    ylabel('R');
    title('PSRF of Power Spectrum at each wavelength');
    print(fig,'-djpeg','-r300',[res_dir,'/psrf_PowSpec.jpg'])
 
    
save([res_dir,'/psrf_draws.txt'],'saveR1','-ASCII');
save([res_dir,'/psrf_stats.txt'],'saveR2','-ASCII');
save([res_dir,'/psrf_statsP.txt'],'saveR3','-ASCII');
 
