function  int_K = integrate2_cov_f(x,z,n,nu,sig2,rho,cov_type,verbose,a,nNodes,TOL,version)

%%%
% integrate2_cov_f  computes the double numerical integral of variance-covariance 
%          matrix of Gaussian Random Field for specified location
%          integrated w.r.t. location (0->y_1),(0->y_2).
%          Covariance model: Matern, nu = smoothness,
%          sig2 = sill,
%          rho = range
% Method:  Version 1: Adaptive Simpson quardature (MATLAB built-in function) 
%          evaluated at every point
%            NOT vectorized
%          Version 2: Gaussian quadrature
%            NOT vectorized
%          Version 3: Gaussian quadrature
%            Partially vectorized to array (small speed-up)
%
% Input: x = x location: sky/latitude
%        y = y location: redshift
%        n = number of locations 
%        nu = order of modified Bessel function of the second kind
%        A = alpha parameter
%        B = beta parameter
%
% Output: 2-dim integral of Covariance matrix K
%%%

eps = 10^(-3); %0; %10^(-3);

if(nargin<12)
    version = 5;  %use Self-written Simpson's rule double integration: 
end %- it is much faster than dblquad and a bit faster than gaussquad.
%NOTE: version 3 has an ERROR in implementation, have not been fixed yet!

if(nargin<11)
    TOL = 10^-7;  %abs error tolerance for dblquad() in version1 and quad2d in diagonal section
end

if(verbose)
    disp('Integrating 2D covariance');
end


    if(verbose)
        tic;
    end
    
if(strcmp(cov_type,'matern'))
    
    if(nargin < 3)
        n = length(x);
        if (n ~= length(y))
            error('Number of (x,y) locations must match.');
        end
    end
    %Autocovariance function parameters
    if(nargin < 4)
        nu = 2;
        sig2 = 0.1/8;
        rho = sqrt(2)/2;
        %B = 0.1; A = 4;
    end
    %number of nodes in Gaussian quadrature or function evaluations for
    %Simpson's Rule
    if(nargin < 10)
        nNodes = 24; 
    end

    int_K = zeros(n,n);
    if (version == 0)
        %Double integration by adaptive MATLAB's function quad2d
        for j=1:(n-1)
            for k=(j+1):n
                %compute 1-D integral
                if(x(j)==0 && x(k)==0 && z(j)==0 && z(k)==0)
                    int_K(j,k) = 0;
                    disp('Warning, K is singular, zero-case');
                else
                    int_K(j,k) = quad2d(@(zj,zk)cov_evalfun(nu,sig2,rho,cov_type,x(j),x(k),zj,zk,1), 0,z(j), 0,z(k),'AbsTol',TOL);
                end
            end
        end
    elseif (version == 1)
        %Double integration by adaptive Simpson's quadrature, MATLAB's function
        for j=1:(n-1)
            for k=(j+1):n
                %compute 1-D integral
                if(x(j)==0 && x(k)==0 && z(j)==0 && z(k)==0)
                    int_K(j,k) = 0;
                    disp('Warning, K is singular, zero-case');
                else
                    int_K(j,k) = dblquad(@(x_,z_)cov_evalfun(nu,sig2,rho,cov_type,x(j),x(k),x_,z_), 0,z(j), 0,z(k),TOL);
                end
            end
        end
    %     %create vectors of lower triangular matrices of ALL combination of
    %     %locations for integration
    %     [z3, x2] = meshgrid(z,x);
    %     [x3, z2] = meshgrid(x,z);
    %     id = tril(x2,-1)==0; %identify upper-triangular index
    %     x2(id)=[]; %turn into vector by removing diagonal and upper-triangle entries
    %     z2(id)=[];
    %     x3(id)=[];
    %     z3(id)=[];
    % 
    %     tmp1 = zeros(length(x2(:)),1);
    %     for j=1:length(x2(:))
    %         disp( 'do lop j=') ; disp(j);
    %         tmp1(j) = dblquad(@(zj,zk)cov_evalfun(nu,sig2,rho,cov_type,x2(j),x3(j),zj,zk,1), 0,z2(j), 0,z3(j));
    %     end
    %     squareform(tmp1)

    elseif(version == 2)
        %version 2:  Double integration by Gaussian quadrature, 
        %dbl for-loop, no speed
        xj = x;
        xk = x;
        zj = z;
        zk = z;
        for j=1:(n-1)        
            [node2, w2] = gaussquad(nNodes, 0, zj(j));
            for k=(j+1):n
                [node1, w1] = gaussquad(nNodes, 0, zk(k));  %Gaussian quadrature with 20 nodes;
                [weight2 weight1] = meshgrid(w2,w1);
                fun = cov_evalfun(nu,sig2,rho,cov_type, xj(j), xk(k), node2, node1);
                int_K(j,k) = sum(sum(weight1 .* weight2 .* fun));
            end
        end
    elseif(version == 3)
        %This method DOES NOT WORK YET!!!!
        %version 3:  Double integration by Gaussian quadrature
        %single for-loop, evaluated as array (Nx25x(25*p)) - need more memory !!!
        xj = x;
        xk = x;
        zj = z;
        zk = z;
        for k=2:n   %for each column of covariance...  
            [node1, w1] = gaussquad(nNodes, 0, zk(k));  %Gaussian quadrature with 25 nodes
            [node2, w2] = gaussquad(nNodes, 0, zj(1:(k-1)) );   %result is matrix of nodes
            %create mesh of all nodes, for each fixed node1 to have all 'sets' of all nodes2
            %'set' j in [1:(k-1)] will give an entry j of intergated covariance
            [node1_g node2_g] = meshgrid(node1,node2(:));
            [w1_g] = meshgrid(w1,w2(:));
            %expand xj entries to match 'sets' of meshed node2
            tmp = repmat(xj(1:(k-1)),size(w2,1),1); 
            [~, xj_g] = meshgrid(node1,tmp(:));
            dont_mesh = 1; %indicator not to mesh locations insice cov_evalfun()
            fun = cov_evalfun(nu,sig2,rho,cov_type, xj_g, xk(k), node2_g, node1_g, dont_mesh);
            inner_int = sum(w1_g .* fun,2);  %sum across node1
            %fill-in upper-diagonal part with integrated covariance by column
            int_K(1:(k-1),k) = sum(w2 .* reshape(inner_int,size(w2,1),[]) );  %sum across node2
        end
    elseif(version == 4)
        %version 4:  Self-implemented Simpson's rule.
        xj = x;
        xk = x;
        zj = z;
        zk = z;
        %Assume: fixed j,k
        %set dimensions
        N=nNodes;
        M=nNodes;
        %get Simpson's weights
        W = dblsimpweights(M,N);
        for j=1:(n-1)
            for k=(j+1):n    
                %This method works!!! and fast.   =)
                %But, integrand must be smooth for this to work well.
                %do the same by hand-written Simpson's rule, 
                %unlike dblquad, method does NOT check against Tolerance
                %mesh y-locations after expanding for integration
                u = linspace(0,zk(k),M+1);
                v = linspace(0,zj(j),N+1);
                [v u] = meshgrid(v,u);
                %evaluate covariance function at u,v
                dont_mesh = 1;
                f = cov_evalfun(nu,sig2,rho,cov_type, xj(j), xk(k), v, u, dont_mesh);
                %compute integral
                int_K(j,k) = zj(j)*zk(k)/(9*M*N)*sum(sum(W.*f));
            end 
        end
    else
        %version 5:  Self-implemented Simpson's rule, Vectorized down to
        %1-for-loop
        xj = x(:);
        xk = xj;
        zj = z(:);
        zk = zj;
        %Assume: fixed j,k
        %set dimensions
        nNodes=24;
        N=nNodes;   %value of 24 gives 10^-5 accuracy to the double integral.
        M=nNodes;   %value of 20 gives 10^-4 accuracy 
        %get Simpson's weights
        W = dblsimpweights(M,N);    
        for k=2:n   %for each column of covariance...

            %This method works!  =)
            %But, integrand must be smooth for this to work well.
            %do the same by hand-written Simpson's rule, 
            %unlike dblquad, method does NOT check against Tolerance
            %mesh y-locations after expanding for integration
            u = linspace(0,zk(k),M+1);
            %generate j 'sets' of v for all j=1:(k-1) 
            v = zeros(k-1,N+1);
            for j=1:(k-1)
                v(j,:) = linspace(0,zj(j),N+1); %result is matrix of node_v
            end
            %create mesh of all nodes, for each fixed node1 to have all
            %'sets' of all node_v
            %'set' j in [1:(k-1)] will give an entry j of intergated covariance
            [v_g u_g] = meshgrid(v,u);
            v_g = reshape(v_g,[],N+1);
            u_g = reshape(u_g,[],N+1);
            W_g = repmat(W,k-1,1);
            %expand xj entries to match 'sets' of meshed node_v
            xj_g = reshape(repmat(xj(1:(k-1))',M+1,N+1) , [],N+1);
            %evaluate covariance function at u,v
            dont_mesh = 1; %indicator not to mesh locations insice cov_evalfun()
            fun = cov_evalfun(nu,sig2,rho,cov_type, xj_g, xk(k), v_g, u_g, dont_mesh);
            %compute integral
            inner_int = sum(W_g .* fun,2);  %sum across node_u
            %unstack 'sets' and
            %fill-in upper-diagonal part with integrated covariance by column
            int_K(1:(k-1),k) = sum(reshape(inner_int,M+1,[]))'.*zj(1:(k-1))*zk(k)/(9*M*N);  %sum across node_v
        end
    end

    %combine upper and lower triangular matrix parts
    int_K = int_K + int_K';

    %fill-in the diagonal terms using Simpson's method
    for j=1:n
        %quad2d is faster and more persise than dblquad for calculating the diagonal terms
        %int_K(j,j) = dblquad(@(yj,yk)cov_evalfun(nu,sig2,rho,cov_type,0,0,yj,yk,1), 0,y(j), 0,y(j),TOL);
        int_K(j,j) = quad2d(@(zj,zk)cov_evalfun(nu,sig2,rho,cov_type,0,0,zj,zk,1), 0,z(j), 0,z(j),'AbsTol', TOL);  
    end
    
    %%% Add correction term to the diagonal of matern-Matrix:
    int_K = int_K + eye(n)*eps;
    
elseif(strcmp(cov_type,'dblexp'))
 
    %Create mesh based on vector elements with index (j,k)
    [xxk,xxj] = meshgrid(x(:),x(:));
    [zzk,zzj] = meshgrid(z(:),z(:));
    
    %Compute the double integral of Double Exponential Covariance : EXACT solution
    int_K = sqrt(2*pi)*exp(-0.5*(xxj-xxk).^2).*( zzj.*normcdf(zzj) - (zzj-zzk).*normcdf(zzj-zzk) - zzk.*normcdf(-zzk) ...
           + sqrt(0.5/pi)*(exp(-0.5*zzj.^2)-1-exp(-0.5*(zzj-zzk).^2)+exp(-0.5*zzk.^2)));
     
elseif(strcmp(cov_type,'inverse'))

    %%% Compute double Integral of Inverse Square covariance
    
    %Create mesh based on vector elements with index (j,k)
    [xxk,xxj] = meshgrid(x(:),x(:));
    [zzk,zzj] = meshgrid(z(:),z(:));
    
    %Compute the integral of Inverse Square Covariance 
    int_K = 1./(xxj-xxk) .* ( -(zzj-zzk).*atan((zzj-zzk)./(xxj-xxk)) + zzk.*atan(zzk./(xxj-xxk)) + zzj.*atan(zzj./(xxj-xxk)) ) ...
            + 1/2 .* ( log((zzj-zzk).^2 + (xxj-xxk).^2) - log(zzk.^2 + (xxj-xxk).^2) - log(zzj.^2 + (xxj-xxk).^2)  + log( (xxj-xxk).^2) );

 %check the diagonal term:
    for j=1:n
         int_K(j,j) = pi/2;
    end
    
elseif(strcmp(cov_type,'uncorr'))
    %Compute the double integral of Double Exponential Covariance : EXACT solution
    int_K = diag(z.^2);
         
elseif(strcmp(cov_type,'simpex') || strcmp(cov_type,'simple') || strcmp(cov_type,'simpdi'))
    %Compute the no integral, simply evaluate the covariance
    int_K = cov_f(x(:),z(:),nu,sig2,rho,cov_type,verbose,a);    
end

    if(verbose)
        toc;
    end
    
    
% %check the diagonal term:
%     for j=1:n
%         [node w] = gaussquad(nNodes, 0, zj(j));
%         [weight2 weight1] = meshgrid(w,w);
%       %  [node1 node2] = meshgrid(node,node);
%         fun = cov_evalfun(nu,A,B,cov_type, 0, 0, node, node);
%         int_K(j,j) = sum(sum(weight1 .* weight2 .* fun));
%     end

if(verbose==2)
    %Check if int_K is positive (semi) definite
    E=sort(eig(int_K));
    S=sort(svd(int_K));
    disp('Check if double integral gives proper result in integrate2_cov_f');
    disp('Eigenvalues and Singular values must be positive and equal in Left-Right of 10 elements:');
    disp([E(1:10) S(1:10)]);                 %  <-  these must be equal!
    disp([E((end-9):end) S((end-9):end)]);   %  <-  these must be equal!
end

% 
%     % integration exaples and comparison
%     %integrate using Adaptive Simpson quadrature
%     tic
%     Q = quad(@(x)cov_evalfun(4,0.1,1,2,x),0,2)
%     toc
% 
%     %integrate using Gauss quadrature
%     tic
%     [node, w] = gaussquad(8, 0, 2);
%     W = sum(w .* cov_evalfun(4,0.1,1,2,node))
%     toc
% 
%     %double integral over a rectangle
%     tic
%     zmin = @(z) z;
%     zmax = @(z) z - 1;
%     Q = quad2d(@(x,z)cov_evalfun(4,0.1,1,2,x),0,1,zmin,zmax) %NO
%     toc
% 
%     %correct double integral:
%     F = @(x,z)0.1.*(1+(x-z).^2).*besselk(2,4*sqrt(1+(x-z).^2));
%     Q = dblquad(F,0,1,0,1)
% 
%     %Result:
%     Q = dblquad(@(x,z)cov_evalfun(4,0.1,1,2,x,z), 0,1,0,1)
%     %and
%     Q = quad(@(x)cov_evalfun(1,4,0.1,1,2,x), 0,1)
% 
% 
%     %double integral with rectangles:   0<E<1,  0<C<1,  zj=zk=1.
%     dx=0.01;
%     x=0:dx:1;
%     z=x;
%     
%     tmp=[];
%     for j=2:length(z),  tmp(j,:) = cov_evalfun(2,4,0.1,0,0,z(j) - 1+eps + x(2:end));  end   %NaN is introduced if you remove eps.
%     %tmp(isnan(tmp))=0.1/8;
%     sum(sum(tmp))
%     sum(sum(tmp))*dx^2
