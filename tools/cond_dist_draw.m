function  f_cond = cond_dist_draw(x_0,x,y_0,y,z_obs, mu_f,sd_e,nu,sig2,rho,Sig,verbose,L,a)

%%%
% cond_dist_draw  draws samples of GRF from conditional multivariate normal
%                 distribution one field given other field
% Input: x_0 = x location: sky/latitude, to-be-predicted
%        x   = x location: sky/latitude, known
%        y_0 = y location: redshift, to-be-predicted
%        y   = y location: redshift, known
%        mu  = whole mean vector
%        f = vector of GRF evaluated at equidistant timepoints 
% Output: vector of integrated GRF
%%%
  
%initialize
M = numel(x_0);
n = numel(x);
N = n+M;
mu = mu_f(ones(M,1)); 
mu_z = y*mu_f;

%%%% Generate observations from the conditional distribution Y|X
%1. find spatial prediction of GRF based on data (conditional expectation)
do_cov = 0;
f_hat = kriging(x_0(:),x,y_0(:),y,mu,mu_z,z_obs,sd_e,nu,sig2,rho,do_cov,Sig,0.000001,a);
S11 = Sig(1:M,1:M);
S12 = Sig(1:M,(M+1):N);
S22 = Sig((M+1):N,(M+1):N);
%l = chol(S22,'lower');
[v,d] = eig(S22);
invS22 = v*(d^-1)*v';
condCov = S11 - S12*invS22*S12';   %conditional covariance matrix

   %Cholesky decomposition  Sigma=U*U'
    sqrtCov = chol(condCov,'lower');


%consider computing the square-root via SVD.
%     [ u,s,v] = svd(condCov);
%     sqrtCov = u*(s.^0.5)*v';

W2 = normrnd(0,1,M,1);  %Generate white noise for random field
f_cond = sqrtCov*W2 + f_hat;



if(~verbose)
disp('Sig_fz . Sig_zz^-1 . z_obs       _____     Sig_fz . Sig_zz^-1 . y_t   ___  f_draw    ___   f_hat');
disp([S12*invS22*z_obs ,  S12*invS22*y,  f_cond, f_hat]);

disp('S12*invS22');
disp(S12*invS22);
disp('z_obs, z_obs-mu_z');
disp([z_obs, z_obs-mu_z]);

disp('Var(f|z)');
disp(S11 + S12*invS22*S12');
end

%disp('');
% %2. generate unconditional spatial process
% W2 = normrnd(0,1,N,1);  %Generate white noise for random field
% 
% Cov_new = [Cov]....
% 
% L = chol(Cov_new,'lower');
% F_tilde = L*W2+[mu; mu_z];      %GRF realization  ***check mean
% z_tilde_gal = F_tilde((M+1):N); %prepare GRF submatrix
% %3. estimate the spatial process based on unconditional spatial process
% do_cov = 0;
% f_tilde_hat = kriging(x_0(:),x,y_0(:),y,mu,mu_z,z_tilde_gal,sd_e,nu,sig2,rho,do_cov,Sig); 
% %4. find difference between siulated process and its prediction
% e = F_tilde(1:M) - f_tilde_hat;
% %5. take a draw from conditional distribution
% f_cond = f_hat + e;