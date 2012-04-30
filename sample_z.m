function  z_t = sample_z(f_tM, s_obs, x_obs, x_gridM, z_gridM, sd_e, Pz, z_grid_Pz, x0,z0,cov_type, verbose,a)

%%%%%%%%%%%%%%%%%%%%%
% sample_y.m   function takes MCMC draws of redshift, z.
%               Methods: numerical integration, cdf sampling
% Input: f_tM  = matrix of GRF at current iteration t, observed on x-z-grid

%        s_obs = observed data, integrated field over z location
%        x_obs = observed x location
%        x_gridM = matrix grid of x locations for GRF
%        z_gridM = matrix grid of z locations for GRF
%        sd_e  = standard deviation of observed measurement error.
%        Pz        = priors of redshift z
%        z_grid_Pz = grid of redshift z for which Pz is defined
%        x0        = mesh-grid of x_obs  vs z for which Pz is defined
%        z0        = mesh-grid of z for which Pz is defined vs x_obs
%        verbose = (T/F) display progress of program
%
% Output: draws  = posterior samples of parameters: redshift (z)
%%%%%%%%%%%%%%%%%%%%%
    
    if(verbose)
        disp('Sampling z_t');
    end
    
    %initialize
    N = length(x_obs);   % Number of observed galaxies
    np = length(z_grid_Pz);
    if( size(z_grid_Pz,1)==1 )
        z_grid_Pz = z_grid_Pz'; %coerse to column vector
    end
    
%     if(verbose)
%         disp('Kriging f0|f on a denser grid acording to redshift prior');
%     end
%     do_cov = 0; do_banding = 0;
%     f_tM_Pz = kriging(x0(:), x_gridM(:), z0(:), z_gridM(:), mu_ff, f_tM(:), 0, do_cov, Sig,invU22, do_banding);
%     f_tM_Pz = reshape(f_tM_Pz,n,N);

    %need f_tM_Pz to be a matrix: length(z_grid_Pz) X length(x_obs)

    if(verbose)
        disp('Interpolating f to another grid');
    end
    if(verbose)
        tic;
    end
    % Interpolate GRF to another grid
    f0_tM = interp2(x_gridM, z_gridM, f_tM, x0, z0);  %use linear interpolation
    if(verbose)
        toc;
    end
    
    % Integrate f along z, get cumulate integral
    if(strcmp(cov_type,'simple') || strcmp(cov_type,'simpex') || strcmp(cov_type,'simpdi'))
        mu_zM = f0_tM;
    else
        mu_zM = integrate_f(z_grid_Pz, f0_tM);
    end

    % Evaluate proportional part of posterior p(z|.)~p(z)*p(s|f,z)
    if(sd_e~=0)
        pdf = Pz' .* normpdf(s_obs(:,ones(np,1))', mu_zM, sd_e);
    else
        pdf = Pz'; %no error noise: use prior directly
    end
    % Compute normalizing consant, sum across z
    norm = sum(pdf, 1);
    
    % Compute normalized pdf of posterior p(z|.)
    pdf = pdf ./ norm(ones(np,1),:);

    % Adjust impossible draws (NaNs) normalized by 0  (good idea??)
    pdf(1:np,norm==0) = 0;
    
    % Compute cdf of posterior p(z|.)
    cdf = cumsum(pdf, 1);
    
    % Sample y using cdf method
    u = unifrnd(0,1,N,1);
    u = u(:,ones(np,1))';
    id = sum(cdf <= u, 1) + 1;   %allow for zeros, i.e. 0.001
    z = [z_grid_Pz; max(z_grid_Pz)];
    z_tmp = z(id);
%    acceptID = norm > 10^-10;   %do not update z_t if normalizing constant is =0.
%    z_t(acceptID) = z_tmp(acceptID);
    z_t = z_tmp;
    
end