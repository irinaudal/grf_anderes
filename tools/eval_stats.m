%function  [delta delta_true d c P P_true xsi stat_foci foci_true stat_ex ex_true stats stats_true] = diagnose_grf(draws, ...
function  [stats stats_true statsP stats_trueP] = eval_stats(draws,IDchain, mu_f, ...
z_true, f_true, xgrid,zgrid, res_dir)

    % Focal-point average, depends on focal size 'q' and number of regions 'K' 
    K = 4;
    q = 4;
    % Focal-point histogram, depends on bin definition: minimum boundary of bin 'lowbin' and bin increment 'dbin' 
    dbin   = 5;
    lowbin = mu_f-1.5*dbin;
    maxbin = mu_f+1.5*dbin;
    % Extremes, depends on threshold 'h' - will be sensitive to model
    h = [3 4 5]+mu_f;
    
    niter = size(draws,1);
    N     = length(z_true);
    M     = length(f_true);
    m1    = length(xgrid);
    m2    = length(zgrid);
    dx    = xgrid(2)-xgrid(1);
    dz    = zgrid(2)-zgrid(1);
    
    % separate variables: rows=MCMC iterations, cols=Variables
    y = draws(:,1:N);
    f = draws(:,(N+1):(N+M));
        
    
%     % density contrast, point-to-point statistic
%     mean_f = mean(f,2);
%     delta = ( f-mean_f(:,ones(M,1)) )./mean_f(:,ones(M,1));
%     delta_true = (f_true-mean(f_true))/mean(f_true);
%     
%     % distance between signals ..for variance convergence check, averaged over pixels, length=niter
%     d = sqrt(1/M * sum( (delta - delta_true(ones(niter,1),:) ).^2 , 2) );

%%%%%%%%%%%%%%%
        %%% Getting fixed wavelength from frequencies from pixel sampling rates
        Fs_x = 1;          % pixels per centimeter
        Fs_z = 1;
        % We can then compute the pixel sizes in each direction as the reciprocal of the pixel sampling rates:
        dxx = dx/Fs_x;     % centimeters per pixel
        dzz = dz/Fs_z;
        % We can then compute a linear distance scale for the X- and Y-axes of the cropped image:
        x = dxx*(0:m1-1)';  % centimeters,   where m1,m2: % pixels
        z = dzz*(0:m2-1)'; 
        % Next, we compute the frequency increments for both X and Y:
        dFx = Fs_x/m1;     % cycles per centimeter
        dFz = Fs_z/m2;
        % Finally, we can create a spatial frequency domain array for both X and Y:
        Fx = (-Fs_x/2:dFx:Fs_x/2-dFx)';     % cycles per centimeter
        Fz = (-Fs_z/2:dFz:Fs_z/2-dFz)';
        
        [mFx, mFz] = meshgrid(Fx,Fz);
        wavelength = sqrt(mFx.^2 + mFz.^2);
        lambda = unique(wavelength);
        nl = length(lambda);

    % Power spectrum of matter field realization
    P = zeros(niter, M);
    for j=1:niter
        f_gridM = reshape(f(j,:),m1,m2);
        F = fft2(  f_gridM  );
        DC = F(1,1);
        shiftF = fftshift(F);
        Amplitude = abs(shiftF); % Amplitude Spectrum
        PowerSpec = Amplitude.^2;  %Power Spectrum
         
        for k=1:nl
            id = (wavelength==lambda(k));
            P(j,k) = sum(PowerSpec(id));
        end
        
    end
    P = P(:,1:nl); %shorten P to proper elements
        P_true = [lambda, zeros(nl,1)]; 
        % True Power Spectrum
        f_gridM = reshape(f_true,m1,m2);
        F = fft2(  f_gridM  );
        DC = F(1,1);
        shiftF = fftshift(F);
        Amplitude_true = abs(shiftF); % Amplitude Spectrum
        PowerSpec_true = Amplitude_true.^2;  %Power Spectrum
         
        for k=1:nl
            id = (wavelength==lambda(k));
            P_true(k,2) = sum(PowerSpec_true(id));
        end
        
        
    
% fftransform = abs(fft2(f_gridM) ).^2;
% image(abs(fft2(f_gridM) ).^2)
% image(log ( abs(fft2(f_gridM) ).^2))
% imagesc(log ( abs(fft2(f_gridM) ).^2))
% imagesc(log ( abs(fftshift(fft2(f_gridM)) ).^2))
% vec=fft2(f_gridM);
% hist(real(vec(:)))
% vec=abs(fft2(f_gridM));
% hist(real(vec(:)))
% hist(log(vec(:)))
    
% imagesc(log ( abs(fftshift(fft2(f_gridM)) ).^2))
% vec=fft2(f_gridM);
% hist(real(vec(:)))
% vec=abs(fft2(f_gridM));
% hist(real(vec(:)))
% hist(log(vec(:)))    

%     P_true = abs(fft(f_true));	  %% absolute value of the fft 
%     P_true = P_true(1:np).^2/np;  %% take the power of positve freq. half
%     
%     % Deviation of power spectrum
%     denom = sum(P_true(ones(niter,1),:));
%      xsi = (P - P_true(ones(niter,1),:))./denom(ones(niter,1),:);
%     %xsi = sum((P - P_true(ones(niter,1),:)).^2,2);
%%%%%%%%%%%%%%%%%%%    

%     % Correlation factor
%     c = zeros(niter,1);
%     for t = 1:niter
%         delta_mean = mean(delta(1:t,:),1);
%         c(t) = sum(delta_true.*delta_mean) / ( sqrt(sum(delta_true.^2)) * sqrt(sum(delta_mean.^2)) );
%     end
%     %correlation coefficient between delta and true delta
%     c_delta= zeros(niter,1); 
%     for kk=1:niter
%         c_delta(kk) = corr(delta(kk,:)',delta_true');
%     end
%     
%     fig=figure();
%     subplot(2,1,1);
%     plot(c)
%     title('Correlation factor between delta and true delta');
%     xlabel('iteration t');
%     ylabel('Cor');
%     
%     subplot(2,1,2);
%     plot(c_delta)
%     title('Correlation coefficient between delta and true delta');
%     xlabel('iteration t');
%     ylabel('r');
%     print(fig,'-djpeg','-r300',[res_dir,'/fig_cor_delta',num2str(IDchain),'.jpg'])
%     close(fig);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Focal-point average, depends on focal size 'q' and number of regions 'K' 
    % randompy select K regions
    % +
    % Histogram of field content: count of field within each color group,
    % per each focal region. f_region /in color_[a,b]

seed = 24312; 
rand('state',seed); % set arbitrary seed for uniform draws
randn('state',seed); % set arbitrary seed for normal draws
    foci = [randsample(m1-q+1,K), randsample(m2-q+1,K)];
    
    % Identify color boundaries (for histogram bins)
    color_bin = lowbin:dbin:maxbin;
    nbin      = length(color_bin);
    [a b] = meshgrid(1:K,1:nbin);
    coords = [a(:) b(:)];
    
    f_com = [f; f_true];
    % Retrieve statistic for each foci separately, for all iterations at once
    stat_foci = zeros(niter+1,K);
    stat_hist = zeros(niter+1,nbin);
    
    for k = 1:K
        % get all pixel pocations for current loci
        id_row = max(1,foci(k,2)-floor(q/2)+1):min(m2,foci(k,2)+floor(q/2));
        id_col = max(1,foci(k,1)-floor(q/2)+1):min(m1,foci(k,1)+floor(q/2));
        id_full = zeros(length(id_col),length(id_row));
        for r = 1:length(id_row)
            id_full(:,r) = r*m1 + id_col-1;
        end
        id_full = id_full(:);
        id_full = id_full(id_full<(size(f_com,2)+1));
        
        id_length = length(id_full);

        % compute average
        stat_foci(:,k) = mean(f_com(:,id_full),2);  
        
        if(k==1)
          % compute histograms
          [H] = histc(f_com(:,id_full),color_bin,2);
          stat_hist(:,:) = H/id_length;
        end
        
%         figure()
%         bar(color_bin,H','BarWidth',1)
    end
    % compute truth
    foci_true = stat_foci(end,:);
    stat_foci(end,:) = [];   %remove the last true row from stats
    hist_true = stat_hist(end,:);
    stat_hist(end,:) = [];   %remove the last true row from stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Extremes integrals, depends on threshold 'h' - will be sensitive to model
    f_com = [f; f_true];
    stat_ex = zeros(niter+1,length(h));
    for j = 1:length(h)
        I = abs(f_com-mu_f(ones(size(f_com))))>h(j);
    %     % how many extremes for each iteration t
    %     nex = sum(I,2);
    %     % index of popularity within each iteration
    %     id = [1; cumsum(nex)+1];
        for t = 1:niter
            stat_ex(t,j) = sum(abs(f_com(t,I(t,:))))*dx*dz; %% Total area of |f|  - no negative integrals
        end
    end
    % compute truth
    ex_true = stat_ex(end,:);
    stat_ex(end,:) = [];   %remove the last true row from stats

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
     % Combine statistics together
     stats = [stat_ex stat_foci stat_hist];
     statsP = P;
     stats_true = [ex_true foci_true hist_true];
     stats_trueP = P_true(:)';
    

    % plot results
    
    %%% Power spectrum - should converge to zero
    fig=figure();
    subplot(2,2,1);
    imagesc(Fx, Fz, log ( PowerSpec ))
    shading flat;
    colorbar('East');
    xlabel('freq_x'); ylabel('freq_z');          
    title('log(Power Spectrum) of f, iter t=1000');

    subplot(2,2,2);
    plotID = niter*(1:10)/10;
    plot(lambda, log( P(plotID,:) ) );
    axis([lambda(1) lambda(nl) min(min(log( P(plotID,1:nl) ))) max(max(log( P(plotID,1:nl) ))) ]);
    xlabel('Wavelength'); ylabel('log(sum(Power Spectrum))');
    title('100th sample of log(sum(Power Spectrum)) of f');

    subplot(2,2,3);
    waveID = 2;
    plot(1:niter, log( P(:,waveID) ) );
    xlabel('Iteration 1'); ylabel('log(sum(Power Spectrum))');
    title(['log(Power Spectrum) of f at Wavelength = ',num2str(lambda(waveID))]);
    subplot(2,2,4);
    waveID = min(50,nl);
    plot(1:niter, log( P(:,waveID) ) );
    xlabel('Iteration, t'); ylabel('log(sum(Power Spectrum))');
    title(['log(Power Spectrum) of f at Wavelength = ',num2str(lambda(waveID))]);

   
    print(fig,'-djpeg','-r300',[res_dir,'/fig_power_spectrum_',num2str(IDchain),'.jpg'])
    close(fig)
    
    %%% Convergence diagnostics - should converge to 1
%     fig=figure();
%     subplot(2,2,1);
%     plot(d)
%     title('Convergence diagnostics for Matter Density Field');
%     xlabel('t, iteration number');
%     ylabel('Eucl.dist b/w samples');
% 
%     subplot(2,2,2);
%     if (niter < 6)
%         xsiID = [1 2 3 4 5];
%     else
%         xsiID = [1:20:500];
%     end
% %     plot(xsi(xsiID,:)')
%     plot(xsi')
%     xlabel('frequency');
%     ylabel('rel.dist b/w Power spectrums');
%     
%     subplot(2,3,4);
%       iter_plot = 1;
%     plot(delta_true, delta(iter_plot,:),'.','MarkerSize',2)
%     hold on;
%     hline = refline(1,0);
%     set(hline,'Color',[0.5 0.5 0.5])
%     xlabel('delta_{true}');
%     ylabel('delta_{MCMC}');
%     title([num2str(iter_plot),' iteration of MCMC']);
%     subplot(2,3,5);
%     if (niter < 6)
%         iter_plot = 2;
%     else
%         iter_plot = 10;
%     end      
%     plot(delta_true, delta(iter_plot,:),'.','MarkerSize',2)
%     hold on;
%     hline = refline(1,0);
%     set(hline,'Color',[0.5 0.5 0.5])
%     xlabel('delta_{true}');
%     ylabel('delta_{MCMC}');
%     title([num2str(iter_plot),' iteration of MCMC']);
%     subplot(2,3,6);
%     if (niter < 6)
%         iter_plot = 5;
%     else
%         iter_plot = 100;
%     end      
%     plot(delta_true, delta(iter_plot,:),'.','MarkerSize',2)
%     hold on;
%     hline = refline(1,0);
%     set(hline,'Color',[0.5 0.5 0.5])
%     xlabel('delta_{true}');
%     ylabel('delta_{MCMC}');
%     title([num2str(iter_plot),' iteration of MCMC']);
%         
%     %addpath('export_fig');
%     %export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/fig_diag.pdf'])
%     %print(fig,'-depsc','-tiff','-r300',[res_dir,'/fig_diag.jpg'])
%     print(fig,'-djpeg','-r300',[res_dir,'/fig_diag_',num2str(IDchain),'.jpg'])
    
end