%%% test function for sampling from posterior distribution of (y,f|z)
%   in order to produce the sample of the Gaussian random field, f.
%%%%%%%%%%%%%%%%%%%%%

% Simulation specs
seed = 2432814; 
rand('state',seed); % set arbitrary seed for uniform draws
randn('state',seed); % set arbitrary seed for normal draws
niter = 100000;
burnin = 95aZ000;
verbose = false; %true; %true %false;
doplots = false; %true  %see analyze_grf.m  how to adjust patrial printing of information
saveimage = true;
email = false;
gelman_rubin_diag = true;%false;
Nchains = 3;

simulate = true;        %simulate MCMC chains for gelman_rubin_diag
readfromfile = false;   %read MCMC chains results from file
examine  = true;        %examine MCMC chains results already in matlab memory


% Set parameters
N = 100;          %Number of observed galaxies (divisible by 4)
m1 = 21;          %Number of locations in 1-D for GRF on a regular grid (for now, same as posterior z-grid)
m2 = 21;
M = m1*m2;       %Number of GRF on a regular grid
xmin = 0;        %minimum value of x-grid
xmax = 4.005;    %maximum value of x-grid
zmin = 0.005;    %minimum value of z-grid (for now, same as posterior z-grid)
zmax = 4.005;    %maximum value of z-grid (for now, same as posterior z-grid)

sd_e = 0.02;%0.02;     %Standard deviation of observed measurement error.
mu_f = 1;        %mean of GRF (constant, for now)
addpath('data');
addpath('tools');
addpath('tools/mcmcstat/');
addpath('tools/mcmcdiag/');
addpath('tools/export_fig/');
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

fixed_z = false;   %'zS';
fixed_f = false;

%set(0,'RecursionLimit',M*2)

% Add new directory to store results
addpath('results');
res_dir = ['results/result_',cov_type,'_',num2str(niter),'_',num2str(N),'_',num2str(m1)];
status = mkdir(res_dir);
addpath(res_dir);
save([res_dir,'/params.txt'],'seed','cov_type','niter','burnin','N','m1','m2','sd_e','mu_f','nu','sig2','rho','-ASCII');

%%%%%%%%%%%%%%%%%%%%%
% simulate_new_data
    %read posterior file
    [id_p, Pz,z_grid_Pz,zS,zPM,zMode,zMed,~,zVar] = posterior_redshift_access(N,Posterior_filename,BPZ_filename,res_dir,verbose,sd_z,type);

    %generate data
    disp('Running a simple test...');
    [x_gridM, x_obs, z_gridM, z_obs, S11, f_gridM, s_obs] ...
        = simulate_grf(N,m1,m2, xmin,xmax,zmin,zmax, sd_e, mu_f, ...
                       zS,niter, nu,sig2,rho, cov_type,...
                       verbose,doplots,a,res_dir);
    disp('All done. :)');

    %Note: z_obs, f_gal and s_tilde_gal are NOT really observed.
    % Create storage variable ob current observed data
    storage = [id_p,s_obs,x_obs,zS,zPM,zMode,zMed,zVar];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (gelman_rubin_diag == false)
    %run Gibbs sampler (once)
    Nchains = 1;
end 

% run_gibbs
if (simulate) 
    [MAINchain] = mult_chains_sim(Nchains, niter, burnin, saveimage, storage, s_obs, x_obs, x_gridM, z_gridM, Pz, z_grid_Pz,... 
                            sd_e, mu_f, S11, nu,sig2,rho, cov_type, verbose,doplots,res_dir,...
                            fixed_z,zS,zPM,zMode,zMed,fixed_f,f_gridM,a);
end  %END simulating MAINchain of MCMC draws

if (readfromfile)  
    [MAINchain] = mult_chains_read(Nchains,res_dir);      
end  %END reading MAINchain of MCMC draws from previously stores file
  
if (gelman_rubin_diag == true) %Gelman & Rubin Diagnostic   

if (examine)
    %Gelman & Rubin Diagnostic   
    [rchain schain pchain saveR1 saveR2 saveR3] = mult_chains_examine(MAINchain, mu_f,z_obs,f_gridM, x_gridM,z_gridM, fixed_z,fixed_f, res_dir);      
end  %END examine MCMC draws: resmc and their statistics
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

if(doplots)
    disp('Checking MCMC draws...');
    if(fixed_z==false)
        MCMCtraceplot(resmc,res_dir); %save trace of all MCMC
    else
        MCMCtraceplot(resmc(:,(N+1):end),res_dir); %save trace of f MCMC
    end
end
        
% if(doplots)
%     few = 1:4;
%     nf = 5;%15;
%     %first_few_draws = resmc(:,1:few);
%     if (strcmp(cov_type,'simpex') || strcmp(cov_type,'simple') || strcmp(cov_type,'simpdi'))
%         plot_mcmcresult(resmc, few, nf, z_grid_Pz, Pz, zS, x_obs,z_obs,s_obs, f_gridM, x_gridM, z_gridM,res_dir,1,[] );
%     else
%         plot_mcmcresult(resmc, few, nf, z_grid_Pz, Pz, zS, x_obs,z_obs,s_obs, f_gridM, x_gridM, z_gridM,res_dir,0,[] );
%     end
%     %few = 1:4 + N;
%     %plot_mcmcresult(resmc, z_grid_Pz, Pz, zS, few);
%     %plot_mcmcresult(resmc, few, nf, z_grid_Pz, Pz, zS, f_gridM, x_gridM, z_gridM,1 );
% end
% few = 1:4;
% nf  = 5;
% plot_mcmcresult(MAINchain{1}, few, nf, z_grid_Pz, Pz, zS, x_obs,z_obs,s_obs, f_gridM, x_gridM, z_gridM,res_dir,1,1);

if(email)
%e-mail me when finished running script
!pwd | mail -s finished_with_gibbs_storage isudaltsova@ucdavis.edu
end

% quit matlab when job is done
%quit;