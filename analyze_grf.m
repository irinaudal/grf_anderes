function  draws = analyze_grf(niter, s_obs, x_obs, x_gridM, z_gridM, Pz, z_grid_Pz,...
                              sd_e, mu_f, S11, nu,sig2,rho, cov_type, verbose,doplots,res_dir,...
                              fixed_z,zS,zPM,zMode,zMed,fixed_f,f_gridM, a)

%%%%%%%%%%%%%%%%%%%%%
% analyze_grf.m   function controls the MCMC draws from joint posterior for Gaussian
%             random field, f.
%             Methods: Gibbs sampler and Metropolis-Hastings algorithm
% Input: 
%        niter = total number of iterations per chain
%        s_obs = observed data, integrated field over z location
%        x_obs = observed x location
%        x_gridM = matrix grid of x locations for GRF
%        z_gridM = matrix grid of z locations for GRF
%        f_gridM = matrix of GRF f
%        Pz         = priors of redshift z
%        z_grid_Pz  = grid of redshift z for which Pz is defined
%        sd_e = Standard deviation of observed measurement error.
%        mu_f = Mean of GRF (constant, for now)
%        S11  = Covariance of GRF f on grid
%        nu   = order of modified Bessel function of the second kind
%        sig2 = sill parameter
%        rho  = range parameter
%        cov_type = 'matern' or 'dblexp' or 'inverse'
%        verbose = (T/F) display progress of program
%
% Output: draws  = posterior samples of parameters (z,f)
%%%%%%%%%%%%%%%%%%%%%

if (verbose)
    disp('Initializing...');
end
% Initialize
N = numel(x_obs);        % Number of observed galaxies
[m2 m1] = size(x_gridM); % Number of z-locations and x-locations for GRF on a grid
M = m1*m2;               % Number of grid evaluations of whole GRF
draws = zeros(niter,N+M);

if(verbose)
    disp('Computing meshgrid of between x_obs and z_posterior');    
end
[x0 z0] = meshgrid(x_obs, z_grid_Pz);

if(verbose)
    disp('Generate starting values');    
end
% Generate starting values ***
if(fixed_f)
    f_tM = f_gridM;
else
    %f_t = mvnrnd(mu_f(ones(M,1)),S11);  %simple starting values for grf, correct model.
    %f_tM = reshape(f_t, m2,m1);
    %f_tM = f_gridM + reshape(normrnd(0,0.1,M,1),m2,m1);

    f_tM = ones(m2,m1)*normrnd(mu_f, 4); % f is a constant/flat field    
end

if(strcmp(fixed_z,'zS'))
    z_t = zS;
    keep_fixed_z=true;
elseif(strcmp(fixed_z,'zPM'))
    z_t = zPM;
    keep_fixed_z=true;
elseif(strcmp(fixed_z,'zMode'))
    z_t = zMode;
    keep_fixed_z=true;
elseif(strcmp(fixed_z,'zMed'))
    z_t = zMed;
    keep_fixed_z=true;
else
    keep_fixed_z=false;
    z_t = normrnd(2,0.5,N,1); %generate a normal random variable    
    %generate from prior distribution of Pz...NOT DONE
end

for iter = 1:niter
    if(any(iter==[2,500])) %modify plotting of GRF f
        doplots=true;
        verbose=false;
    else
        doplots=false;
        verbose=false;
    end
    
    %sample f|.    
    if(~fixed_f)
        f_tM = sample_f(z_t,f_tM, s_obs, x_obs, x_gridM, z_gridM, mu_f, sd_e, S11, nu,sig2,rho,cov_type, verbose,doplots,iter,res_dir,a);
    end
    %sample z|.    
    if(~keep_fixed_z)
        z_t = sample_z(f_tM, s_obs, x_obs, x_gridM, z_gridM, sd_e, Pz, z_grid_Pz, x0,z0,cov_type, verbose,a);
    end
    
    %update draws
    draws(iter,1:N) = z_t;
    draws(iter,(N+1):(N+M)) = f_tM(:);

    if(verbose)
         disp(['Done with iteration #: ', num2str(iter)]);
    end
    if(~mod(iter,50))
        disp('Iter '); disp(iter);
    end
end

