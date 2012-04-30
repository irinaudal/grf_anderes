function cov = cov_evalfun(nu,sig2,rho,cov_type,xj,xk,zj,zk,dont_mesh)

%%%% 
% K_covfun  evaluates the matern class autocovariance K(|t-s|) for given 
%           set of coordinates t=(xj,zj) and r=(xk,zk).
%           This is a Tool function for integration of covariance matrix
%           over y-coordinate
%
% Input: nu = order of Modified Bessel function of the second kind
%        a = inner constant multiple of the variables in the Matern model
%        b = outer constant in the Materm model
%        xj = vector (mx1) of fixed coordinates of first dimension of t
%        xk = value (1x1) of fixed coordinates of first dimension of r
%        zj = vector (mx1) of variable coordinates of second dimention of t
%        z(zk) = vector (nx1) of variable coordinates of second dimention
%                of r to be integrated over
%        dont_mesh = 1 <- Do NOT change the dimension of xj,xk,zj,zk
%
% Output:  vector of K(|t-s|) for each combination of coordinates of s and t

if(nargin < 9)
    dont_mesh = 0;     % by default: do the meshing of coordinates
end

if(dont_mesh==1) %dimension of x and y matrices is already pre-fixed, do not change.
    xjm = xj;
    zjm = zj;
    zkn = zk;

    Dx = xjm-xk;
    Dz = zjm-zkn;

elseif(dont_mesh==-1) % special mesh, mesh all coordinates
    [xkn,xjm] = meshgrid(xk(:),xj(:));
    [zkn,zjm] = meshgrid(zk(:),zj(:));

    Dx = xjm-xkn;
    Dz = zjm-zkn;
    
else %create mesh grid between xj,xk and zj,zk, special meshing: zk may be a single value
    if(length(xj) > 1)
        [xjm] = meshgrid(xj,zk);
    else
        xjm = xj;
    end
    [zjm zkn] = meshgrid(zj,zk);

    Dx = xjm-xk;
    Dz = zjm-zkn;
end
DT = sqrt(Dx.^2+Dz.^2);

% if(dont_mesh==1) %original coordinates, used to be: nargin>9 ...
% 
%     if(strcmp(cov_type,'matern'))
%         cov = (sig2*2*sqrt(nu)^nu/(gamma(nu)*rho^nu)).*(sqrt(Dx.^2+zj.^2)).^nu.*besselk(nu,(2*sqrt(nu)/rho).*sqrt(Dx.^2+zj.^2));
% 
%     elseif(strcmp(cov_type,'simple'))    
%         cov = exp(-sig2*sqrt(Dx.^2+Dz.^2));
%     end
%     
% elseif(dont_mesh==0) %meshed coordinates

    if(strcmp(cov_type,'matern') || strcmp(cov_type,'simple'))
        %%% Compute the Matern covariance
        % Autocovariance function simple parameters 
        A = 2*sqrt(nu)/rho;
        B = (sig2*2*sqrt(nu)^nu/(gamma(nu)*rho^nu));
     
        % Evaluate the covariance for distances DT
        cov = B.*DT.^nu.*besselk(nu,A.*DT);
        cov(isnan(cov)) = sig2;
        
    elseif(strcmp(cov_type, 'simpex'))  %'expone'))   %'simple')) 
        % Autocovariance function simple parameters 
        A = 2*sqrt(nu)/rho;
 
        %Compute the Exponential Covariance
        cov = sig2*exp(-A*DT);
        
        cov(DT==0) = sig2;
    elseif(strcmp(cov_type,'simpdi'))    
        %Compute the Block-diagonal Covariance
        cov = zeros(size(DT));
        cov(DT<0.5) = rho*DT(DT<0.5);
        
        cov(DT==0) = sig2;
    end
    %cov=reshape(cov,:,1)   
% end
