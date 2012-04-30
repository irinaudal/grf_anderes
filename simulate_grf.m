function [x_gridM, x_obs, z_gridM, z_obs, S11, f_eqM, s_gal_wNoise] = ...
    simulate_grf(N,m1,m2, xmin,xmax,zmin,zmax, sd_e, mu_f, ...
                 zS,niter, nu,sig2,rho, cov_type,verbose,doplots,a,res_dir) 

%####################################################################################################
%# simulate_grf   Function simulates data for redshift-random field project
%                 Generate stationary centered Gaussian Random Field on R^2
%                 625 observation locations in (0,1)^2
%
% Input: N  = Number of observed galaxies
%        m1 = Number of x-locations for GRF on a grid
%        m2 = Number of z-locations for GRF on a grid
%        xmin = minimum value of x-grid
%        xmax = maximum value of x-grid
%        zmin = minimum value of z-grid (for now, same as posterior z-grid)
%        zmax = maximum value of z-grid (for now, same as posterior z-grid)
%        sd_e = Standard deviation of observed measurement error.
%        mu_f = Mean of GRF (constant, for now)
%        zS   = z-location(spectroscopic redshift) for observed instance of GRF
%        niter = total number of iterations (for plotting only) 
%        nu   = order of modified Bessel function of the second kind
%        sig2 = sill parameter
%        rho  = range parameter
%        cov_type = 'matern' or 'dblexp' or 'inverse'
%        verbose = (T/F) display progress of program
%        doplots = (T/F) create plots for results
%
% Output: Simulated parameters: 
%         x_gridM, x_obs, z_gridM, z_obs, S11, f_eqM, s_gal_wNoise
%####################################################################################################

%% Initialize
M = m1*m2;

method = 1;  %generate random field (f,z) together at once. Need to compute main covariance structure Sig.
%method = 2;   % generate random field (f), then compute (z) by using this (f)
              % and measurement error.

%Prepare equiditant locations for GRF f
x_grid = linspace(xmin,xmax,m1)';  %equidistant timepoints  %seq(1,4,length=m1)
z_grid = linspace(zmin,zmax,m2)';  %equidistant timepoints  
[x_gridM z_gridM] = meshgrid(x_grid,z_grid); 

%Prepare observation locations for observed field s
x_obs = unifrnd(xmin,xmax,N,1);%sort(unifrnd(xmin,xmax,N,1));
x_obs = sort(x_obs);
z_obs = zS;       %spectroscopic (true) redshift, NOT sorted

%%
if (method == 1)
%% Method 1: generate random field (f,z) together at once.

    if(verbose)
        disp('Calculating the covariance matrices of Gaussian Random Field data, f&z');
    end
    %Covariance matrix Sigma of random field
    S11 = cov_f(x_gridM(:),z_gridM(:),nu,sig2,rho,cov_type,verbose,a);
    if(verbose == 3)
        disp('Check that Sigma_f (or S11) matrix is pos.def. in simulate_grf'); 
        E=sort(eig(S11));
        S=sort(svd(S11));
        disp('Eigenvalues and Singular values must be positive and equal in Left-Right of 10 elements:');
        disp([E(1:10) S(1:10)]);                 %  <-  these must be equal!
        disp([E((end-9):end) S((end-9):end)]);   %  <-  these must be equal!
    end
    
    S12 = integrate1_cov_f(x_gridM(:),x_obs,z_gridM(:),z_obs,M,N,nu,sig2,rho,cov_type,verbose,a);
    S22 = integrate2_cov_f(x_obs,z_obs,N,nu,sig2,rho,cov_type,verbose,a);

 %    S22 = S22 + eye(N)*sd_e^2;     %Add noise variance 
 %note: noise will be added as a separate component
 
    Sig = [S11 S12; S12' S22];

    if(verbose)
        disp('Check that Sigma matrix is pos.def. in simulate_grf'); 
        E=sort(eig(Sig));
        S=sort(svd(Sig));
        disp('Eigenvalues and Singular values must be positive and equal in Left-Right of 10 elements:');
        disp([E(1:10) S(1:10)]);                 %  <-  these must be equal!
        disp([E((end-9):end) S((end-9):end)]);   %  <-  these must be equal!
    end
    
    %Cholesky decomposition  Sigma=L*L'
    L = chol(Sig,'lower');

    if(verbose)
        disp('Generating Gaussian Random Field data');
    end
    
    %Generate white noise for random field
    W = randn((M+N),1);     

    %Generate white noise for additive noise
    epsilon = sd_e*randn(N,1);    
  
    %Gaussian Random Field realization
    if(strcmp(cov_type,'simple') || strcmp(cov_type,'simpex'))
        mu = [mu_f(ones(M,1)); mu_f(ones(N,1))]; %stacked vector of means for f and s
    else
        mu = [mu_f(ones(M,1)); z_obs*mu_f]; %stacked vector of means for f and s
    end
    fs = L*W + mu;

    %Separate GRF into components of f (stored as matrix) and s (as vector)
    f_eqM = reshape(fs(1:M),m2,m1);    
    s_gal_wNoise = fs((M+1):(M+N)) + epsilon;

    %Integrate f_eq over z : required for plots only within this function
    if(strcmp(cov_type,'simple') || strcmp(cov_type,'simpex') || strcmp(cov_type,'simpdi'))
        s_tilde_eqM = f_eqM;
    else
        s_tilde_eqM = integrate_f(z_gridM,f_eqM);
    end
    
elseif (method == 2)
%% Method 2: generate random field (f), then compute (s) by using this (f)
%  and measurement error.
 
    if(verbose)
        disp('Calculating the covariance matrix of Gaussian Random Field data, f');
    end
    %Covariance matrix Sigma of GRF on grid
    S11 = cov_f(x_gridM(:),z_gridM(:),nu,sig2,rho,cov_type,verbose);

    %Cholesky decomposition  Sigma=L'*L
    L = chol(S11,'lower');

    if(verbose)
        disp('Generating Gaussian Random Field data');
    end

    %Generate white noise for random field
    Wf = normrnd(0,1,M,1);
    Wz = normrnd(0,sd_e,N,1);

    %Gaussian Random Field realization
    f = L*Wf + mu_f(ones(M,1));
    f_eqM = reshape(f, m2,m1);
      
    %inegrate f_gal over z_gal
    s_gal_wNoise = integrate_f(z_obs,f_eqM,x_obs,x_gridM,z_gridM) + Wz;

    %Integrate f_eq over z : required for plots only within this function
    s_tilde_eqM = integrate_f(z_gridM,f_eqM);

end
%%
%    %construct s_tilde at galaxy locations, as well
%    f_gal = interp2(x_gridM,z_gridM,f_eqM,x_obs,z_obs);
%    s_tilde_gal = integrate_f(z_obs,f_gal);

%% Print results
if(doplots)
        if(verbose)
            disp('Working on plots...');
        end
    %plot results
    fig = figure();
       subplot(2,2,1);    
    pcolor(x_gridM,z_gridM,s_tilde_eqM);   %surf() plot viewed from above
    shading flat;      
    xlabel('x'); ylabel('z'); zlabel('s^{tilde}(x,z)');  
    title('\int GRF over z, Equidist Grid, no noise');
    colorbar('East')
    h = caxis;
       subplot(2,2,2);           
    stem(x_obs,zS);
    hold on;
    scatter(x_obs,z_obs,20,s_gal_wNoise,'filled'); %plot: a=vector of location, b=vector of redshift, size, color choice, 
    hold off;
    xlabel('x'); ylabel('z');
    title(['\int GRF over z, True.Redshift vs U(0,4), with e ~ N(0,',num2str(sd_e),'^2)']);
    caxis(h);
       subplot(2,2,3);    
    pcolor(x_gridM,z_gridM,f_eqM);   %surf() plot viewed from above
    shading flat;      
    xlabel('x'); ylabel('z'); zlabel('f(x_0,z_0)');  
    title('True GRF, Equidist Grid, no noise');
    h2 = caxis;
    colorbar('East')
   
    export_fig(fig,'-painters', '-r300', '-q101', [res_dir,'/data.pdf'])
%    print(fig,'-depsc','-tiff','-r300',[res_dir,'/data.eps'])
    close(fig);
 end     
%%
   