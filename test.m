%%% test function for sampling from posterior distribution of (y,f|z)
%   in order to produce the sample of the Gaussian random field, f.
%%%%%%%%%%%%%%%%%%%%%


% Simulation specs
seed = 243231512; 
rand('state',seed); % set arbitrary seed for uniform draws
randn('state',seed); % set arbitrary seed for normal draws
niter = 200;
burnin = 0;
verbose = false; %true; %true %false;
doplots = false; %true  %see analyze_grf.m  how to adjust patrial printing of information
saveimage = true;
email = false;
gelman_rubin_diag = true;
Nchains = 3;

% Set parameters
N = 100;          %Number of observed galaxies (divisible by 4)
m1 = 21;          %Number of locations in 1-D for GRF on a regular grid (for now, same as posterior z-grid)
m2 = 21;
M = m1*m2;       %Number of GRF on a regular grid
xmin = 0;        %minimum value of x-grid
xmax = 4.005;    %maximum value of x-grid
zmin = 0.005;    %minimum value of z-grid (for now, same as posterior z-grid)
zmax = 4.005;    %maximum value of z-grid (for now, same as posterior z-grid)

sd_e = 0.1;%0.02;     %Standard deviation of observed measurement error.
mu_f = 4;        %mean of GRF (constant, for now)
Posterior_filename = 'mockPosterior1.txt';
BPZ_filename = 'mockGen1.bpz';
type = 'bpz';  % 'normal';   %what prior distribution to use?
nu = 2;        %Autocovariance parameters
rho = sqrt(2)/2;  % 5;  % sqrt(2)/2/10;
sig2 = 6;  % 1.25/100;% 1/100;
%TOL = 10^(-5);
%v_f = 0.01;   %sd for proposal of MH step
%v_z = 0.1;
a = 1/1000;
sd_z = 0.1; %redshift SD for normal prior distribution


cov_type = 'simple'; %'uncorr';  %'inverse'; 'dblexp'; 'matern'; 'simple'; 'simpex'; 'simpdi';

fixed_z = false;
fixed_f = false;

%set(0,'RecursionLimit',M*2)

% Add new directory to store results
res_dir = ['result_',cov_type,'_',num2str(niter),'_',num2str(N),'_',num2str(m1)];
status = mkdir(res_dir);
addpath(res_dir);
save([res_dir,'/params.txt'],'seed','cov_type','niter','burnin','N','m1','m2','sd_e','mu_f','nu','sig2','rho','-ASCII');

%%%%%%%%%%%%%%%%%%%%%
% simulate_new_data
    %read posterior file
    [id_p, Pz,z_grid_Pz,zS,zPM,zMode,zMed,~,zVar] = posterior_redshift_access(N,Posterior_filename,BPZ_filename,res_dir,verbose,sd_z,type);
    %z_grid_Pz = fix(z_grid_Pz*10000)/10000;

    %generate data
    disp('Running a simple test...');
    [x_gridM, x_obs, z_gridM, z_obs, S11, f_gridM, s_obs] ...
        = simulate_grf(N,m1,m2, xmin,xmax,zmin,zmax, sd_e, mu_f, ...
                       zS,niter, nu,sig2,rho, cov_type,...
                       verbose,doplots,a,res_dir);
    disp('All done. :)');

    %Note: z_obs, f_gal and s_tilde_gal are NOT really observed.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_gibbs
if (gelman_rubin_diag == false)
    %run Gibbs sampler
    disp('Running a simple MCMC test...');
    [ resmc ] = analyze_grf(niter, s_obs, x_obs, x_gridM, z_gridM, Pz, z_grid_Pz,... 
                            sd_e, mu_f, S11, nu,sig2,rho, cov_type, verbose,doplots,res_dir,fixed_z,zS,fixed_f,f_gridM,a);
    disp('All done with MCMC. :)');

    %remove burnin
    if(burnin>0)
        resmc(1:burnin,:) = [];
    end        

elseif (gelman_rubin_diag == true)

%Gelman & Rubin Diagnostic
addpath('mcmcstat/');
    
    MAINchain = cell(Nchains);
    dchain    = cell(Nchains);
    schain    = cell(Nchains);
    
    for jj=1:Nchains
        %run Gibbs sampler
        disp('Running a simple MCMC test...');
        [ resmc ] = analyze_grf(niter, s_obs, x_obs, x_gridM, z_gridM, Pz, z_grid_Pz,... 
                            sd_e, mu_f, S11, nu,sig2,rho, cov_type, verbose,doplots,res_dir,fixed_z,zS,fixed_f,f_gridM,a);
        disp(['All done with MCMC of chain: ',num2str(jj),'. :)']);

        %remove burnin
        if(burnin>0)
            resmc(1:burnin,:) = [];
        end
                
        if(saveimage)
            %store the mcmc results to file
            save([res_dir,'/mcmcdraws',num2str(jj),'.txt'], 'resmc', '-ASCII');
            storage = [id_p,s_obs,x_obs,zS,zPM,zMode,zMed,zVar];
            save([res_dir,'/data.txt'], 'storage', '-ASCII');
            save([res_dir,'/true_f.txt'],'f_gridM','-ASCII');
            clear('storage');
        end
        
        % compute statistics
        disp('Checking MCMC draws...');
        z_true = z_obs'; 
        f_true = f_gridM(:)';
%    [delta delta_true d c P P_true xsi stat_foci foci_true stat_ex ex_true stats stats_true] = diagnose_grf(resmc, ...
        [stats stats_true delta delta_true] = diagnose_grf(resmc, ...
            z_true, f_true, x_gridM(1,:),z_gridM(:,1), res_dir);
        
        %store all chains
        MAINchain{jj} = resmc;
        dchain{jj}    = delta;
        schain{jj}    = stats;

    end
    
    %Make Diagnostic plots (script)
    gelman_diag
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%

if(saveimage)
    %store the mcmc results to file
     save([res_dir,'/mcmcdraws.txt'], 'resmc', '-ASCII');
     storage = [id_p,s_obs,x_obs,zS,zPM,zMode,zMed,zVar];
     save([res_dir,'/data.txt'], 'storage', '-ASCII');
     save([res_dir,'/truth.txt'],'f_gridM','-ASCII');
    clear('storage');
end






%     disp('Checking MCMC draws...');
%     z_true = z_obs'; 
%     f_true = f_gridM(:)';
% %    [delta delta_true d c P P_true xsi stat_foci foci_true stat_ex ex_true stats stats_true] = diagnose_grf(resmc, ...
%     [stats stats_true delta delta_true] = diagnose_grf(resmc, ...
%             z_true, f_true, x_gridM(1,:),z_gridM(:,1), res_dir);
        
if(doplots)
    few = 1:4;
    nf = 5;%15;
    %first_few_draws = resmc(:,1:few);
    if (strcmp(cov_type,'simpex') || strcmp(cov_type,'simple') || strcmp(cov_type,'simpdi'))
        plot_mcmcresult(resmc, few, nf, z_grid_Pz, Pz, zS, x_obs,z_obs,s_obs, f_gridM, x_gridM, z_gridM,res_dir,1,[] );
    else
        plot_mcmcresult(resmc, few, nf, z_grid_Pz, Pz, zS, x_obs,z_obs,s_obs, f_gridM, x_gridM, z_gridM,res_dir,0,[] );
    end
    %few = 1:4 + N;
    %plot_mcmcresult(resmc, z_grid_Pz, Pz, zS, few);
    %plot_mcmcresult(resmc, few, nf, z_grid_Pz, Pz, zS, f_gridM, x_gridM, z_gridM,1 );

end

if(email)
%e-mail me when finished running script
!pwd | mail -s finished_with_gibbs_storage isudaltsova@ucdavis.edu
end

% quit matlab when job is done
%quit;