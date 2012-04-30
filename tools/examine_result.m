% perform analysis on results

addpath('result1100_101_51/');

Posterior_filename = 'mockPosterior1.txt';
BPZ_filename = 'mockGen1.bpz';

pars = load('params.txt');
data = load('data.txt');
resmc = load('mcmcdraws.txt');
f_gridM = load('truth.txt');

%storage = id_p,z_obs,x_obs,zS,zPM,zMode,zMed,zVar
id_p  = data(:,1);
z_obs = data(:,2);
x_obs = data(:,3);
zS    = data(:,4);
[~, Pz_y] = hdrload(Posterior_filename);
Pz_y = Pz_y(id_p,:);
Pz_y(:,1) = [];            %remove ID column of Posterior        


%initialize other variables
niter = pars(1) - pars(2);
burnin = pars(2);
N    = pars(3); 
m1   = pars(4);
m2   = pars(5); 
sd_e = pars(6);
mu_f = pars(7);
nu   = pars(8);
sig2 = pars(9); 
rho  = pars(10);
xmin = 0;        %minimum value of x-grid
xmax = 4;        %maximum value of x-grid
ymin = 0.005;    %minimum value of y-grid (for now, same as posterior y-grid)
ymax = 4.005;    %maximum value of y-grid (for now, same as posterior y-grid)
y_grid_Pz = linspace(ymin,ymax,401)';    %Redshift grid for Posterior
x_grid = linspace(xmin,xmax,m1)';  %equidistant timepoints  %seq(1,4,length=m1)
y_grid = linspace(ymin,ymax,m2)';  %equidistant timepoints  
[x_gridM y_gridM] = meshgrid(x_grid,y_grid); 


    %Integrate f_eq over y
    z_tilde_eqM = integrate_f(y_gridM,f_gridM);

disp('Checking data...');
    %plot results
    fig = figure();
       subplot(2,2,1);    
    pcolor(x_gridM,y_gridM,z_tilde_eqM);   %surf() plot viewed from above
    shading flat;      
    xlabel('x'); ylabel('y'); zlabel('z^{tilde}(x,y)');  
    title('\int GRF over y, Equidist Grid, no noise');
    colorbar('East')
    h = caxis;
       subplot(2,2,2);           
    stem(x_obs,zS);
    hold on;
    scatter(x_obs,zS,20,z_obs,'filled'); %plot: a=vector of location, b=vector of redshift, size, color choice, 
    hold off;
    xlabel('x'); ylabel('y');
    title(['\int GRF over y, True.Redshift vs U(0,4), with e ~ N(0,',num2str(sd_e),'^2)']);
    caxis(h);
       subplot(2,2,3);    
    pcolor(x_gridM,y_gridM,f_gridM);   %surf() plot viewed from above
    shading flat;      
    xlabel('x'); ylabel('y'); zlabel('f(x,y)');  
    title('True GRF, Equidist Grid, no noise');
    h2 = caxis;
    colorbar('East')
   
    
disp('Checking MCMC draws...');
few = 1:4;
nf = 20;
%first_few_draws = resmc(:,1:few);
plot_mcmcresult(resmc, few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,1);
few = 11:14;
%plot_mcmcresult(resmc, few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,1);
few = 91:94;
%plot_mcmcresult(resmc, few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,1);
plot_mcmcresult(resmc, few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,2,1 );
%plot_mcmcresult(resmc, few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,2,4 );
plot_mcmcresult(resmc(1:19,:), few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,2,1 );
%plot_mcmcresult(resmc(1:38,:), few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,2,1 );
%plot_mcmcresult(resmc(1:991,:), few, nf, y_grid_Pz, Pz_y, zS, f_gridM, x_gridM, y_gridM,2,55 );

