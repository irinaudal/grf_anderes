function  K = cov_f(t1,t2,nu,sig2,rho,cov_type,verbose,a)

%%%
% cov_f    computes the variance-covariance matrix of
%          Gaussian Random Field for specified location
%          Covariance model: nu = smoothness,
%          sig2 = sill,
%          rho = range
%
% Input: t1 = column vector of every x location paired with y
%        t2 = column vector of every y location paired with x
%        n = number of locations 
%        nu = order of modified Bessel function of the second kind
%        sig2 = sill parameter
%        rho = range parameter
%        cov_type = 'matern' or 'dblexp' or 'inverse'
%
% Output: Covariance matrix K, (n x n) matrix
%%%

n = length(t1(:));
if (n ~= length(t2(:)))
    error('Two location vectors must be of the same length!');
end

% Autocovariance function parameters
if(nargin < 4)
    nu = 2;
    sig2 = 0.1/8;
    rho = sqrt(2)/2;
    %B = 0.1; A = 4;
end

if(verbose)
    tic;
end

% Compute distances b/w each location
if(verbose)
    disp('Computing distances between locations...');
end
DT = pdist([t1(:) t2(:)]);
if(verbose)
    disp('Evaluating cov_f...');
end


if(strcmp(cov_type,'matern') || strcmp(cov_type,'simple'))
 
     %%% Compute the Matern covariance
     % Autocovariance function simple parameters 
     A = 2*sqrt(nu)/rho;
     B = (sig2*2*sqrt(nu)^nu/(gamma(nu)*rho^nu));
     
     % Evaluate the covariance for distances DT
     if(verbose)
         disp('Evaluating matern class with besselk()...');
     end
     K = B*DT.^nu.*besselk(nu,A*DT);
     K = squareform(K);
     K = K + eye(n)*sig2;    %the diagonal element is limit{dt->0}[K] = sig2

elseif(strcmp(cov_type,'dblexp'))
    
    %%% Compute Double Exponential covariance
    K = exp(-0.5/a*DT.^2);
    K = squareform(K);
    K = K + eye(n)*1;  %the diagonal element is = 1

elseif(strcmp(cov_type,'inverse'))

    %%% Compute Inverse Square covariance
    K = DT.^(-2);
    K = squareform(K);
    K = K + eye(n)*1;  %the diagonal element is = 1
    
elseif(strcmp(cov_type,'uncorr'))
 
    %Compute the Double Exponential Covariance
    K = eye(n);
         
elseif(strcmp(cov_type,'simpex')) % 'expone')) %'simple'))
     % Autocovariance function simple parameters 
     A = 2*sqrt(nu)/rho;
 
    %Compute the Exponential Covariance
    K = sig2*exp(-A*DT);
    K = squareform(K);
    K = K + eye(n)*sig2;    %the diagonal element is limit{dt->0}[K] = sig2

elseif(strcmp(cov_type,'simpdi'))    
    %Compute the Block-diagonal Covariance
    K = zeros(size(DT));
    K(DT<0.5) = rho*DT(DT<0.5);

    K = squareform(K);
    K = K + eye(n)*sig2;    %the diagonal element is limit{dt->0}[K] = sig2
end

if(verbose)
    toc;
end
% %     %Check if Sigma is positive (semi) definite
% %     E=eig(K);
% %     S=sort(svd(K));
% %     disp([E(1:10) S(1:10)]);   %  <-  these must be equal!
