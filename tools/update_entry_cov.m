function  S = update_entry_cov(S,J,t1,t2,nu,A,B,integral)

%%%
% update_entry_cov  updates J-th row and column of the variance-covariance 
%          matrix of Gaussian Random Field for specified location
%          Covariance model: Matern, nu=2
%          alpha = 2*sqrt(nu)/rho = 4,
%          beta = sigma^2*alpha^2/[2^(nu-1) Gamma(nu)] = 0.1
% Input: S = covariance matrix
%        j = row/column index to be updated 
%        t1 = vector of x location: sky/latitude
%        t2 = vector of y location: redshift
%        integral = {0 = no integration}, {1 = 1D integration}, 
%                   {2 = 2D integration}
% Output: Covariance matrix S
%%%

%Autocovariance function parameters
if(nargin < 5)
nu = 2;
B = 0.1;
A = 4;
end

n = length(t1);

if(integral == 0)
    %Compute distances b/w each location for jth entry
    dt = [t1(J*ones(n,1)), t2(J*ones(n,1))] - [t1, t2]; 
    DT = sqrt(dt(:,1).^2 + dt(:,2).^2);    %vector norms of each dt

    %Autocovariance function
    K = B*DT.^nu.*besselk(nu,A*DT);
    K(DT==0) = nu*B/(A)^2;    %the diagonal element is limit{dt->0}[K] = 2
    S(J,:)=K;
    S(:,J)=K;

elseif (integral == 1)
    
    [node, w] = gaussquad(20, 0, yk(J));  %Gaussian quadrature with 7 nodes;
    int_K = sum(w(:,ones(n,1)) .* K_covfun(nu,A,B, xj, xk(J), yj, node)); 
    S(:,J)=int_K;
        
    for k=1:m
        [node, w] = gaussquad(20, 0, yk(k));  %Gaussian quadrature with 7 nodes;
        int_K = sum(w .* K_covfun(nu,A,B, xj(J), xk(k), yj(J), node)); 
        S(J,k)=int_K;
    end
    
    
elseif (integral == 2)
        
end