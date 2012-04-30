function  s_tilde_result = integrate_f(z,f,x,x_locf,z_locf)
                 
%%%
% cov_grf  computes the vector (matrix) of 2D Gaussian Random Field 
%          integrated over one coordinate z from 0 up to a constant specified in z. 
%          Simple version of code generates integral up to every element in
%          z, similar to a cdf
% Input: z = vector (or matrix) of z location: redshift
%        f = vector (or matrix) of GRF evaluated at equidistant timepoints 
%        x = vector (or matrix) of x location, size must be same as z
%        x_locf = matrix grid of x locations of f
%        z_locf = matrix grid of z locations of f       
%        
% Output: vector (matrix) of integrated GRF
%%%
  
    m = size(z,1);
    if(nargin<3) %integrate f up to all z, for all x-locations

        %f is a vector
        if(isvector(f) && sqrt(length(f))==m) 
            dz = z(3)-z(2);
            f_eqM = reshape(f,m,m);
            s_tilde_eqM = cumsum(f_eqM,1)*dz;
            res = s_tilde_eqM(:);

        %f is a matrix
        elseif(size(f,1)==m) 
            dz = z(3)-z(2);
            res = cumsum(f,1)*dz;
        end
    else     %integrate f up to select (x,z)-locations

        %f is matrix evaluated on a grid
        %result is known up to grid, so interpolate between grid
        grid1 = x_locf(1,:); 
        grid1 = grid1(:); %coerse to column vector
        grid2 = z_locf(:,1);
        grid2 = grid2(:); 
        m1 = length(grid1);
        m2 = length(grid2);

        x=x(:)';        %coerse to row vector
        z=z(:)';
        m_long = length(x);

        f_gal = interp2(x_locf,z_locf,f,x(:),z(:));
        s_tilde_eqM = integrate_f(z_locf,f);

        %code is vectorized:            
        %find vertical location, compute distance from nearest vertical eq.location
        idz = sum(grid2(:,ones(m_long,1)) <= z(ones(m2,1),:),1);   %z_locf must be sorted
        if(isemptz(idz))
            idz = 1;
        end
        dzz = z' - grid2(idz);

        %snap to nearest grid horizontal location
        idx = sum(grid1(:,ones(m_long,1)) <= x(ones(m1,1),:),1);
        % Integrate f over dzz, using s_tilde_eqM as reference up to dzz part
        id = sub2ind(size(s_tilde_eqM), idz, idx);  %linear indexing
        s_tilde = s_tilde_eqM(id)' + f_gal.*dzz;  

        %coerse result to dimensions equivalent to the dimension of entry
        if(m ~= m_long)
            res = reshape(s_tilde,m,[]);
        else
            res = s_tilde;
        end        
    end
    s_tilde_result = res;
end