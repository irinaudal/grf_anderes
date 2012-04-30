function res = trapezoid(f,u,du)

%This function will attempt to integrate the covariance function for all
%possible of (x,y) locations
N=length(u);

trap = zeros(N-1,1);

for j=1:(N-1)
    if(j==1)
        trap(j) = du*(f(1)+f(j+1))/2;
    else
        trap(j) = du*(f(1) + 2*sum(f(2:j)) + f(j+1))/2;
    end
end  

res = [u [0; trap]];