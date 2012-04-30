function [prior] = prior_y_pdf(y, y_grid, Pz_y, do_log) 
    
%####################################################################################################
% prior_y_eval   Function calculates the prior density at given redshift y,
%                using the numerical result of the prior Pz_y; vectorized.
%
% Input: y      = current value of redshift, want to get p(y)
%        y_grid = redshift-locations on a grid at which Pz is defined
%        Pz_y   = redshift posterior from file on a grid
%        log    = (T/F) return log-value of regular
%
% Output: prior
%####################################################################################################

N = length(y);
prior = zeros(N,1);
for j = 1:N
  idy = sum(y_grid <= y(j)); 
  if(idy == 0)
      prior(j) = 0;
  else
      prior(j) = Pz_y(j,idy);
  end
end

if(do_log==true)
  prior = log(prior);
end