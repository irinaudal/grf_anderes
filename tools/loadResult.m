    % %Gelman & Rubin Diagnostic
    addpath('mcmcstat/');
    res_dir = 'test_mcmc';

    
%save([res_dir,'/params.txt'],'seed','niter','burnin','N','m1','m2','sd_e','mu_f','nu','sig2','rho','-ASCII');

%Re-make data
        [ params ] = load([res_dir,'/params.txt']);



% Simulation specs
seed = params(1); 
rand('state',seed); % set arbitrary seed for uniform draws
randn('state',seed); % set arbitrary seed for normal draws
niter = params(2);
burnin = params(3);
verbose = true; %true %false;
doplots = true;  %see analyze_grf.m  how to adjust patrial printing of information
saveimage = true;
email = false;
gelman_rubin_diag = true;
Nchains = 3;

% Set parameters
N = params(4);          %Number of observed galaxies (divisible by 4)
m1 = params(5);          %Number of locations in 1-D for GRF on a regular grid (for now, same as posterior z-grid)
m2 = params(6);
M = m1*m2;       %Number of GRF on a regular grid
xmin = 0;        %minimum value of x-grid
xmax = 4.005;    %maximum value of x-grid
zmin = 0.005;    %minimum value of z-grid (for now, same as posterior z-grid)
zmax = 4.005;    %maximum value of z-grid (for now, same as posterior z-grid)

sd_e = params(7);%0.02;     %Standard deviation of observed measurement error.
mu_f = params(8);        %mean of GRF (constant, for now)
Posterior_filename = 'mockPosterior1.txt';
BPZ_filename = 'mockGen1.bpz';
type = 'bpz';  % 'normal';   %what prior distribution to use?
nu = params(9);        %Autocovariance parameters
rho = params(11);  % 5;  % sqrt(2)/2/10;
sig2 = params(10);  % 1.25/100;% 1/100;
%TOL = 10^(-5);
%v_f = 0.01;   %sd for proposal of MH step
%v_z = 0.1;
a = 1/1000;
sd_z = 0.1; %redshift SD for normal prior distribution


cov_type = 'simple'; %'uncorr';  %'inverse'; 'dblexp'; 'matern'; 'simple'; 'simpex'; 'simpdi';

fixed_z = false;
fixed_f = false;

%set(0,'RecursionLimit',M*2)


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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
 % READ MCMC output   
    MAINchain = cell(Nchains,1);
    dchain    = cell(Nchains);
    schain    = cell(Nchains);
    
    for jj=1:Nchains
        %run Gibbs sampler
        disp('Reading data result for MCMC test...');
        [ resmc ] = load([res_dir,'/mcmcdraws.txt']);
        
        
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