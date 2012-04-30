function W = dblsimpweights(m,n)
% Produces the m by n matrix of weights for Simpson?s rule
% for double integrals
% Inputs: m -- number of intervals in the row direction.
% must be even.
% n -- number of intervals in the column direction.
% must be even.
% Output: W -- a (m+1)x(n+1) matrix of the weights
if (rem(m,2)~=0 || rem(n,2)~=0)
error('m and n must be even')
end

w = ones(m+1,1);
for i = 2:m
    if rem(i,2)==0
        w(i)=4;
    else
        w(i)=2;
    end
end
u = w;

w = ones(n+1,1);
for i = 2:n
    if rem(i,2)==0
        w(i)=4;
    else
        w(i)=2;
    end
end
v = w;

W = u*v';