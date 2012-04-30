function  int_K = integrate1_cov_f(xj,xk,zj,zk,m,n,nu,sig2,rho,cov_type,verbose,nNodes)

%%%
% integrate1_cov_f  computes the numerical integral of variance-covariance 
%          matrix of Gaussian Random Field for specified location
%          integrated(0->zk) w.r.t. location zk.
%          Covariance model: Matern, nu = smoothness,
%          sig2 = sill,
%          rho = range
% Method:  version 1: Adaptive Simpson quardature (MATLAB built-in function) 
%          evaluated at every point
%            NOT vectorized
%          version 2: Gaussian quadrature
%            IS vectorized
%
% Input: xj location: sky/latitude of GRF f
%        xk location: sky/latitude of GRF s
%        zj location: redshift of GRF f
%        zk location: redshift of GRF s
%        m = number of fixed locations of xj
%        n = number of integrated locations of zk
%        nu = order of modified Bessel function of the second kind
%        sig2 = sill parameter
%        rho = range parameter
%        cov_type = 'matern' or 'dblexp' or 'inverse'
%
% Output: 1-dim integral of Covariance matrix K
%%%
    if(verbose)
        disp('Integrating 1D covariance');
    end

    if(verbose)
        tic;
    end

if(strcmp(cov_type,'matern'))

    %%% Compute the Matern covariance

    % version = 2;  %use Gaussian quadrature - it is much faster, b/c/ vectorized.
    if(nargin < 5)
        m = length(xj);
        if (m ~= length(zj))
            error('Locations (xj,zj) must have the same length.');
        end
        n = length(xk);
        if (n ~= length(zk))
            error('Locations (xk,zk) must have the same length.');
        end
    end
    %Autocovariance function parameters
    if(nargin < 7)
        nu = 2;
        sig2 = 0.1/8;
        rho = sqrt(2)/2;
        %B = 0.1; A = 4;
    end
    %Number of nodes in Gaussian quadrature
    if(nargin < 12)
        nNodes = 30; 
    end
    
    int_K = zeros(m,n);
 
    % for k=1:n
    %     [node, w] = gaussquad(nNodes, 0, zk(k));  %Gaussian quadrature;
    %     int_K(:,k) = sum(w(:,ones(m,1)) .* cov_evalfun(nu,sig2,rho,cov_type, xj, xk(k), zj, node)); 
    % end

    % % if(sum(xj)==sum(xk))  %symmetry holds: K is mxm 
    % %     t1 = xk;
    % %     t2 = zj;
    % %     int_K = zeros(m,m);
    % %     for j=1:(m-1)
    % %         for k=(j+1):m
    % %             %compute 1-D integral
    % %             if(t1(j)==0 && t1(k)==0 && t2(j)==0 && t2(k)==0)
    % %                 int_K(j,k) = 0;
    % %                 disp('Warning, K is singular, zero-case');
    % %             else
    % %                 int_K(j,k) = quad(@(x)cov_evalfun(nu,sig2,rho,cov_type,t1(j),t1(k),t2(j),x), 0,t2(k));
    % %             end
    % %         end
    % %     end
    % % 
    % %     %combine upper and lower triangular matrix parts
    % %     int_K = int_K + int_K';
    % % 
    % %     %fill-in the diagonal terms
    % %     for j=1:n
    % %         int_K(j,j) = quad(@(x)cov_evalfun(nu,sig2,rho,cov_type,0,0,t2(j),x), 0,t2(j));
    % %     end
    % % else    %symmetry does not hold: K is mxn
        int_K = zeros(m,n);
        %%% version 1
    %    if(version == 1)
            for j=1:m
                for k=1:n
                    %compute 1-D integral
                    if(xj(j)==0 && xk(k)==0 && zj(j)==0 && zk(k)==0)
                        int_K(j,k) = 0;
                    else
                        int_K(j,k) = quad(@(x)cov_evalfun(nu,sig2,rho,cov_type,xj(j),xk(k),zj(j),x), 0,zk(k));
                    end
                end
            end
    %     else    
    %         %%% version 2
    %         for k=1:n
    %             [node, w] = gaussquad(20, 0, zk(k));  %Gaussian quadrature with 7 nodes;
    %             int_K(:,k) = sum(w(:,ones(n,1)) .* cov_evalfun(nu,sig2,rho,cov_type, xj, xk(k), zj, node)); 
    %             %weight = w(:,ones(n,1));
    %             %fun = cov_evalfun(nu,sig2,rho,cov_type, xj, xk(k), zj, node);
    %             %int_K(1:m,k) = sum(weight .* fun);
    %         end
    %     end
    % end

elseif(strcmp(cov_type,'dblexp'))

    %%% Compute Double Exponential covariance
    
    %Create mesh based on vector elements with index (j,k)
    [xxk,xxj] = meshgrid(xk,xj);
    [zzk,zzj] = meshgrid(zk,zj);
    
    %Compute the integral of Double Exponential Covariance : EXACT solution
    int_K = sqrt(2*pi)*exp(-0.5*(xxj-xxk).^2).*( normcdf(zzj) - normcdf(zzj-zzk) );
    
elseif(strcmp(cov_type,'inverse'))

    %%% Compute Inverse Square covariance
    
    %Create mesh based on vector elements with index (j,k)
    [xxk,xxj] = meshgrid(xk,xj);
    [zzk,zzj] = meshgrid(zk,zj);
    
    %Compute the integral of Inverse Square Covariance 
    int_K = (xxj-xxk) .* (atan((zzj-zzk)./(xxj-xxk)) - atan(-zzk./(xxj-xxk)));

elseif(strcmp(cov_type,'uncorr'))
 
    %Compute the integral of Double Exponential Covariance : EXACT solution
    int_K = zeros(length(zj), length(zk));

elseif(strcmp(cov_type,'simpex') || strcmp(cov_type,'simple') || strcmp(cov_type,'simpdi'))
 
    %Create mesh based on vector elements with index (j,k)
    [xxk,xxj] = meshgrid(xk,xj);
    [zzk,zzj] = meshgrid(zk,zj);
    
    %Compute the no integral
    int_K = cov_evalfun(nu,sig2,rho,cov_type,xxj,xxk,zzj,zzk,1);

end

    if(verbose)
        toc;
    end

% % 
% %     % integration exaples and comparison
% %     %% version 1: integrate using Adaptive Simpson quadrature
% %     tic
% %     Q = quad(@(x)cov_evalfun(4,0.1,1,2,x),0,2)
% %     toc
% % 
% %     %integrate using Gauss quadrature
% %     tic
% %     [node, w] = gaussquad(8, 0, 2);
% %     W = sum(w .* cov_evalfun(4,0.1,1,2,node))
% %     toc
% % 
% %     %double integral over a rectangle
% %     tic
% %     zmin = @(z) z;
% %     zmax = @(z) z - 1;
% %     Q = quad2d(@(x,z)cov_evalfun(4,0.1,1,2,x),0,1,ymin,ymax) %NO
% %     toc
% % 
% %     %correct double integral:
% %     F = @(x,z)0.1.*(1+(x-z).^2).*besselk(2,4*sqrt(1+(x-z).^2));
% %     Q = dblquad(F,0,1,0,1)
% % 
% %     %Result:
% %     Q = dblquad(@(x,z)cov_evalfun(4,0.1,1,2,x,z), 0,1,0,1)
% %     %and
% %     Q = quad(@(x)cov_evalfun(1,4,0.1,1,2,x), 0,1)
% % 
% % 
% %     %double integral with rectangles:   0<E<1,  0<C<1,  zj=zk=1.
% %     dx=0.001;
% %     x=0:dx:1;
% %     z=x;
% %     tmp=[];
% %     for j=2:length(z),  tmp(j,:) = cov_evalfun(4,0.1,1,2,z(j) - 1+x(2:end)); end
% %     sum(sum(tmp))
% %     sum(sum(tmp))*dx^2
