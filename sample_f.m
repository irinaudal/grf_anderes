function  f_tM_new = sample_f(z_t,f_tM, s_obs, x_obs, x_gridM, z_gridM, mu_f, sd_e, ...
                        S11, nu,sig2,rho,cov_type, verbose,doplots,iter,res_dir,a)

%%%%%%%%%%%%%%%%%%%%%
% sample_f.m   function takes MCMC draws of GRF, f.
%              Methods: cond.expectation of multivariate normal distr
%
% Input: z_t   = vector of z locations at current iteration
%        f_tM  = matrix of GRF at current iteration t, observed on x-z-grid
%        s_obs = observed data, integrated field over z location
%        x_obs = observed x location
%        x_gridM = matrix grid of x locations for GRF
%        z_gridM = matrix grid of z locations for GRF
%        mu_f  = mean of GRF (constant, for now)
%        sd_e  = standard deviation of observed measurement error.
%        S11   = covariance matrix of GRF f on grid
%        nu    = parameter of Matern cov, smoothness, order of Bessel function
%        sig2  = parameter of Matern cov, sill
%        rho   = parameter of Matern cov, range
%        cov_type = 'matern' or 'dblexp' or 'inverse'
%        verbose = (T/F) display progress of program
%
% Output: f_tM_new  = posterior samples of parameters: field (f) over reg.grid
%%%%%%%%%%%%%%%%%%%%%
    
method = 1;  % sample directly from conditional normal
%method = 2;  % use conditioning simulation using kriging: requires
              %knodledge of full covariance Sig [MNxMN].

    if(verbose)
        disp('Sampling f_t');
    end

    % Initialize
    N = numel(x_obs);   % Number of observed galaxies
    [m2 m1] = size(x_gridM); % Number of y-locations and x-locations for GRF on a grid
    M = m1*m2;            % Number of grid evaluations of whole GRF

if(method == 1)
    % Set-up covariance matrices
    S12 = integrate1_cov_f(x_gridM(:),x_obs,z_gridM(:),z_t,M,N,nu,sig2,rho,cov_type,verbose,a);
    S22 = integrate2_cov_f(x_obs,z_t,N,nu,sig2,rho,cov_type,verbose,a); 
    S22 = S22 + eye(N)*sd_e^2;     %Add noise variance 
    U = chol(S22);    %upper triangular Cholesky, L=U'.
    invL = inv(U)';   %MATLAB is faster and more stable inverting U.
    invS22 = invL'*invL;

    % Get conditional parameters
    if(strcmp(cov_type,'simple'))
        mu_bar = mu_f(ones(M,1)) + S12*invS22*(s_obs - mu_f(ones(N,1)));
    else
        mu_bar = mu_f(ones(M,1)) + S12*invS22*(s_obs - mu_f*z_t);
    end
    Sig_bar = S11 - S12*invS22*S12';
 
%    tic;
    L = chol(Sig_bar, 'lower'); %lower triangular Cholesky, L=U'.
    
    % Generate GRF from conditional normal distr
    W = normrnd(0,1,M,1);     
    f_t_new = L*W + mu_bar;
    
%    toc;
    
   if(verbose)
        disp('Check that Sigma_bar matrix is pos.def. in sample_f'); 
        E=sort(eig(Sig_bar));
        S=sort(svd(Sig_bar));
        disp('Eigenvalues and Singular values must be positive and equal in Left-Right of 10 elements:');
        disp([E(1:10) S(1:10)]);                 %  <-  these must be equal!
        disp([E((end-9):end) S((end-9):end)]);   %  <-  these must be equal!
   end
 %   newS = L'*L;
 %   tic;
 %   f_t_new = mvnrnd(mu_bar,newS);
 %   toc;
    
    if(verbose)
        disp('reshaping GRF to matric form');
    end
    f_tM_new = reshape(f_t_new, m2,m1);
    
elseif(method == 2)
    % Set-up covariance matrices
    S12 = integrate1_cov_f(x_gridM(:),x_obs,z_gridM(:),z_t,M,N,nu,sig2,rho,cov_type,verbose,a);
    S22 = integrate2_cov_f(x_obs,z_t,N,nu,sig2,rho,cov_type,verbose,a); 
    S22 = S22 + eye(N)*sd_e^2;     %Add noise variance 
    Sig = [S11 S12; S12' S22];

    if(verbose)
        disp('Check that Sigma matrix is pos.def.');
        e=eig(Sig);
        sigval=sort(diag(e));
        disp('minimum eigen value = ');
        disp(sigval(1:20));
        [ ~,s,v] = svd(Sig);
        sigval=sort(diag(s));
        disp('minimum simgular value = ');
        disp(sigval(1:20));
        U = chol(S22); %upper triangular Cholesky, L=U'.
        invL = inv(U)'; %MATLAB is faster and more stable inverting U.
        [ ~,s,v] = svd([S11 - S12*(invL'*invL)*S12']);
        sigval=sort(diag(s));
        disp('minimum simgular value = ');
        disp(sigval(1:20));
    end

    %Cholesky decomposition  Sigma=L*L'=U'U.
    L = chol(Sig,'lower');

    %subtract the mean
    s_obs_minus_mean = s_obs - z_t;
    
    % Generate observations from the conditional distribution f|z
    f_t_new = cond_dist_draw(x_gridM,x_obs,z_gridM,z_t,s_obs,mu_f,sd_e,nu,sig2,rho,Sig,cov_type,verbose,L,a);
    f_tM_new = reshape(f_t_new, m2,m1);
%     f_t_new = cond_dist_draw(x_gridM,x_obs,z_gridM,z_t,s_obs_minus_mean,0,sd_e,nu,sig2,rho,Sig,L);
%     f_tM_new = reshape(f_t_new + mu_f, m2,m1);
end

% check results
    if(doplots)
        if(verbose)
            disp('Working on plots...');
        end
        %disp('f_new: '); disp(f_tM_new(:));
        
        %plot results
        fig=figure();
        subplot(1,2,1);           
        pcolor(x_gridM,z_gridM,f_tM);   %surf() plot viewed from above
        shading flat;      
        xlabel('x'); ylabel('z');
        title(['f at iteration t=',num2str(iter)]);
        h = caxis;
        colorbar('East');

        subplot(1,2,2); 
        pcolor(x_gridM,z_gridM,f_tM_new);   %surf() plot viewed from above
        shading flat;      
        xlabel('x'); ylabel('z');
        title(['f at iteration t+1=',num2str(iter+1)]);
        caxis(h);
        
        export_fig(fig,'-painters', '-r300', '-q101', [res_dir,'/fig',num2str(iter),'.pdf'])
        close(fig);
    end 

end