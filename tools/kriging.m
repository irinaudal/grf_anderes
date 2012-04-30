function  [f_hat K] = kriging(x0,x,z0,z,mu_f,mu_s,s,sd_e,nu,sig2,rho,cov_type,verbose,do_cov_est,Sig,do_banding,TOL,a)

%%%
% kriging  computes the conditional mean (value of Gaussian Random
%          Field for specified location) E(f(x0,z0) | s(x,z))
% Input: x0 = x location: sky/latitude, to-be-predicted
%        x  = x location: sky/latitude, known
%        z0 = z location: redshift, to-be-predicted
%        z  = z location: redshift, known
%        mu_b = whole mean vector of both field-on-grid and s-obs
%        s    = vector of known values of random field
%        sd_e = the s.d. of added noise, to be added to the
%               diagonal of Sig 
%        do_cov = indicator to re-estimate covariance matrix, if the
%                 locations have changed from original calculation
%                 = 1, then do not give Sig matrix
%                 = 0, then give Sig matrix
%        Sig = covariance matrix for all locations
%        invU22 = inverse of uper triangular matrix from cholesky(S22)
%        do_banding = indicator to do banding, in case Sig matrix is too large
%        TOL = cut-off tolerance for seting entries of Sig to 0 around band
% Output: f(x0,z0) = predicted values of vector of random field
%         K = (optional) covariance matrix of the field  
%%%

if (nargin < 16)
    do_banding = 0;
end

%Merge locations, to-be-predicted and known, into vector
%t1com = [x0; x];
%t2com = [y0; y];
M = numel(x0);  %number of locations  to-be-predicted
N = numel(x); %number of locations  of conditional field

%Covariance matrix Sigma of random field
S11 = Sig(1:M,1:M);
if(do_cov_est==0)
    S12 = Sig(1:M,(M+1):(N+M));
    S22 = Sig((M+1):(N+M),(M+1):(N+M));
    
elseif(do_cov_est==1 || nargin < 8)
    S12 = integrate1_cov_f(x0,x,z0,z,M,N,nu,sig2,rho,cov_type,verbose,a);
    S22 = integrate2_cov_f(x,z,N,nu,sig2,rho,cov_type,verbose,a); 
    S22 = S22 + eye(N)*sd_e^2;     %Add noise variance 
end 
if(do_banding) %set some entries of Sigma to zero, away from diagonal band.
   band_width = sum(sqrt( (x-x(ones(N,1))).^2 + (z-z(ones(N,1))).^2 ) <= TOL);
   S22 = S22-triu(S22,band_width+1)-tril(S22,-(band_width+1));
end
    U = chol(S22); %upper triangular Cholesky, L=U'.
    invL = inv(U)'; %MATLAB is faster and more stable inverting U.
%    S=sparse(A)
%    C=inv(S)
%    B=full(C)

f_hat = S12 * (invL'*invL) * (s-mu_s) + mu_f;    %cond.exp. E(f(x0,z0) | s(x,z))
K = [S11 S12; S12' S22];
