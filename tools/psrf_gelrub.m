function [R,var_hat,W]=psrf_gelrub(chain)%,a,s,e)
%PSRF  Gelman-Rubin-Brooks psrf MCMC convergence diagnostic
% Gelman-Rubin-Brooks potential scale reduction factor MCMC convergence diagnostic
% Returns (length of total sequence (1-a)% interval)/
%         (average length of within sequence intervals)
%
% Reference:  Gelman,...,Rubin, Bayesian Data Analysis, 2004, p296-297 
%
% Parallel chains must be given in chain{1}, chain{2}, ...
%
% Date: March 11 2012

nchain = 0;
if iscell(chain)
  nchain = length(chain);
  [nsimu,npar1] = size(chain{1});
elseif isnumeric(chain)
  error('Chain must be a cell array with at least 2 chains.');
end

if nchain<2
  error('Chain must be a cell array with at least 2 chains.');
end

% if nargin < 4 | isempty(e); e = nsimu; end
% if nargin < 3 | isempty(s); s = fix(e/2)+1; end
% nsimu = e-s+1;
% 
% 
% if nargin < 2
%   a = 0.1;
% end
% p = [a/2,1-a/2];

mean_j = zeros(nchain,npar1);
s2_j   = zeros(nchain,npar1);
for j=1:nchain
    mean_chain = mean(chain{j},1);
    var_chain = var(chain{j});
    mean_j(j,:) = mean_chain;
    s2_j(j,:) = var_chain;    
end
mean_tot = mean(mean_j,1);
B        = nsimu/(nchain-1) * sum( (mean_j - mean_tot(ones(nchain,1),:) ).^2, 1);
W        = mean(s2_j,1);
var_hat  = (nsimu-1)/nsimu * W + 1/nsimu *B ;
V        = var_hat + B/(nsimu*nchain);
R        = V./W;

% disp(['B/n', num2str(B/nsimu)]);
% disp('************');
% disp(['W = ',num2str(W)]);
% disp('************');
% disp(['var_hat = ',num2str(var_hat)]);
% disp('************');
% disp(['V = ',num2str(V)]);
% disp('************');
% disp(['R = ',num2str(R)]);
% disp('************');

