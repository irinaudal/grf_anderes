%%% This file performs testing of gelman/rubin diagnostic when random
%%% multivariate samples of GRF are drawn

% Simulation specs
seed = 2431326151; 
rand('state',seed); % set arbitrary seed for uniform draws
randn('state',seed); % set arbitrary seed for normal draws
niter = 1000;
burnin = 0;
Nchains = 3;
verbose = false;

% Set parameters
res_dir = 'test_mcmc_r';
status = mkdir(res_dir);
addpath(res_dir);
%Gelman & Rubin Diagnostic
addpath('mcmcstat/');
m1 = 5;          %Number of locations in 1-D for GRF on a regular grid (for now, same as posterior z-grid)
m2 = 5;
M = m1*m2;       %Number of GRF on a regular grid
xmin = 0;        %minimum value of x-grid
xmax = 3.995;    %maximum value of x-grid
zmin = 0.005;    %minimum value of z-grid (for now, same as posterior z-grid)
zmax = 3.995;    %maximum value of z-grid (for now, same as posterior z-grid)

sd_e = 0.1;%0.02;     %Standard deviation of observed measurement error.
mu_f = 4;        %mean of GRF (constant, for now)
nu = 2;        %Autocovariance parameters
rho = sqrt(2)/2;  % 5;  % sqrt(2)/2/10;
sig2 = 6;  % 1.25/100;% 1/100;
cov_type = 'simple'; %'uncorr';  %'inverse'; 'dblexp'; 'matern'; 'simple'; 'simpex'; 'simpdi';

pixelPlot = 20;


%Prepare equiditant locations for GRF f
x_grid = linspace(xmin,xmax,m1)';  %equidistant timepoints  %seq(1,4,length=m1)
z_grid = linspace(zmin,zmax,m2)';  %equidistant timepoints  
[x_gridM z_gridM] = meshgrid(x_grid,z_grid); 

%Covariance matrix Sigma of random field
S11 = cov_f(x_gridM(:),z_gridM(:),nu,sig2,rho,cov_type,verbose,9999);

%Cholesky decomposition  Sigma=L*L'
L = chol(S11,'lower');
  
    
MAINchain  = cell(Nchains,1);
SHORTchain = cell(Nchains,1);
mychain    = cell(Nchains,1);
schain    = cell(Nchains,1);
draws = zeros(niter,M);

for ichain = 1:Nchains
    for iter = 1:niter

        %Generate white noise for random field
        W = randn(M,1);     

        %Generate white noise for additive noise
        epsilon = sd_e*randn(M,1);    

        %Gaussian Random Field realization
        mu = mu_f(ones(M,1)); %stacked vector of means for f
        f = L*W + mu;

        %update draws
        draws(iter,:) = f;

        if(~mod(iter,1000))
            disp(['Completed Iter ',num2str(iter)]);
        end
    end
    
    MAINchain{ichain} = draws;
    SHORTchain{ichain} = draws(:,pixelPlot);
    %store the mcmc results to file
    save([res_dir,'/mcmcdraws',num2str(ichain),'.txt'], 'draws', '-ASCII');
end


iterations=[3:(niter-1)];
lengthR = length(iterations);
saveR = zeros(lengthR,M +1);  % M = # diagnostics
saveR2 = zeros(lengthR,1 +1);  % 1 = # diagnostics

for loop=1:lengthR
         
    ii = iterations(loop);

    mychain{1} = MAINchain{1}(1:ii,:);
    mychain{2} = MAINchain{2}(1:ii,:);
    mychain{3} = MAINchain{3}(1:ii,:);
    %Gelman & Rubin Diagnostic
    R  = psrf(mychain); %need all 3 chains

    schain{1} = SHORTchain{1}(1:ii,:);
    schain{2} = SHORTchain{2}(1:ii,:);
    schain{3} = SHORTchain{3}(1:ii,:);
    R2 = psrf(schain); %need all 3 chains
    
    saveR(loop,:) = [ii R];
    saveR2(loop,:) = [ii R2];
end   


    %Separate GRF into components of f (stored as matrix) and s (as vector)
    f_eqM = reshape(f,m2,m1);  

    %plot results
    fig = figure();
    pcolor(x_gridM,z_gridM,f_eqM);   %surf() plot viewed from above
    shading flat;      
    xlabel('x'); ylabel('z'); zlabel('f(x_0)');  
    title('True GRF, Equidist Grid, no noise');
    colorbar('East')

    fig = figure();
       subplot(2,1,1);
    plot(saveR(:,2:end))  %1st column are IDs
    hold on;   
    xlabel('Iteration, t');
    ylabel('R');
    title('SPRF of all Pixels');
       subplot(2,1,2);
    plot(saveR2(:,2))%1st column are IDs
    hold on;   
    xlabel('Iteration, t');
    ylabel('R');
    title(['SPRF of Pixel ',num2str(pixelPlot)]);


    print(fig,'-djpeg','-r600',[res_dir,'/testing',num2str(4),'.jpg'])
    addpath('export_fig');
    export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/testing',num2str(33),'.pdf'])

    
    
    fig = figure();
       subplot(2,1,1);
    plot(saveR(1:500,2:end))  %1st column are IDs
    hold on;   
    xlabel('Iteration, t');
    ylabel('R');
    title('SPRF of all Pixels');
       subplot(2,1,2);
    plot(saveR2(1:500,2))%1st column are IDs
    hold on;   
    xlabel('Iteration, t');
    ylabel('R');
    title(['SPRF of Pixel ',num2str(pixelPlot)]);


    print(fig,'-djpeg','-r600',[res_dir,'/testingShort',num2str(4),'.jpg'])
    addpath('export_fig');
    export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/testingShort',num2str(33),'.pdf'])
    