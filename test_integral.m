%function [] = test_integral(nu,sig2,rho,xmin,xmax,ymin,ymax,N,idj,idk)

%This function will attempt to integrate the covariance function for all
%possible of (x,y) locations using many numerical integration procedures

N = 10000;
xmin = 0;        %minimum value of x-grid
xmax = 4;        %maximum value of x-grid
zmin = 0.005;    %minimum value of z-grid (for now, same as posterior z-grid)
zmax = 4.005;    %maximum value of z-grid (for now, same as posterior z-grid)

maxTry = 10;  %should be even. total number of different integrals to evaluate and different points and compare. 

linsp=floor(linspace(1,N,maxTry));
for jj=1:maxTry
for ii=1:maxTry

%Define indices to evaluate integral
if(nargin<9)
    idj = linsp(jj);
    idk = linsp(ii);
end
  

addpath(['DERIVESTsuite']);


%Prepare equiditant locations
x_grid = linspace(xmin,xmax,N)';  %equidistant timepoints  %seq(1,4,length=m1)
z_grid = linspace(zmin,zmax,N)';  %equidistant timepoints  
%[zj_gridM zk_gridM] = meshgrid(z_grid,z_grid); 

xj = x_grid;
xk = x_grid;
zj = z_grid;
zk = z_grid;

dont_mesh = 1;

% f = cov_evalfun(nu,sig2,rho,xj(idj),xk(idk),yj(idj),yk,dont_mesh);
xx2 = (xj(idj)-xk(idk))^2;
f = exp(-0.5.*(xx2 + (zj(idj)-zk).^2));

u=zk;
du = zk(3) - zk(2);

% F = @(x)cov_evalfun(nu,sig2,rho,xj(idj),xk(idk),zj(idj),x);
F = @(x)( sqrt(2*pi) * exp(-0.5.*xx2) * (normcdf(x)-normcdf(x-zj(idj))) );

f = feval(F,u);

%try Reimann's sum
%res.rsum = [u cumsum(f)*du];

%try Trapezoid Rule
% res.trap = trapezoid(f,u,du);
res.trap = [u cumtrapz(u,f)];


%try adaptive Simpson's rule, from Matlab, up to tolerance
Q = zeros(N,1);
for i=1:N
    Q(i) = quad(F, 0,u(i),10^-8);
end
res.simp = [u Q];

%try adaptive Lobatto quadrature
Q = zeros(N,1);
for i=1:N
    Q(i) = quadl(F, 0,u(i),10^-8);
end
res.loba = [u Q];

%try adaptive Gauss-Kronrod quadrature    %DOES NOT WORK!!!!
% Q = zeros(N,1);
% errbnd = Q;
% for i=1:N
%     [Q(i),errbnd(i)] = quadgk(F, 0,u(i));
% end
% res.qakr = [u Q];

%try Gauss quadrature
Q = zeros(N,1);
for i=1:N
    [node, w] = gaussquad(10, 0, u(i));  %Gaussian quadrature with 10 nodes;
%     Q(i) = sum(w .* cov_evalfun(nu,sig2,rho, xj(idj), xk(idk), yj(idj), node,dont_mesh)); 
    Q(i) = sum(w .* feval(F,node)); 
end
res.gaus = [u Q];

%other approximations
% Flog = @(x)(-1)*log(cov_evalfun(nu,sig2,rho,xj(idj),xk(idk),zj(idj),x));
Flog = @(x)(-1)*log( sqrt(2*pi) * exp(-0.5.*xx2) * (normcdf(x)-normcdf(x-zj(idj))) );


%try Laplace approx
peak = zj(idj);
ddf = hessian(Flog,peak);
lim1 = (u - peak)*sqrt(ddf);
lim2 = -peak*sqrt(ddf);
Q = sqrt(2*pi/ddf).*feval(F, peak).*( normcdf(lim1)-normcdf(lim2) );
res.lapl = [u Q];


%plot results
figure(jj);
subplot(2,maxTry/2+1,1);
hold on;
plot(u,res.trap(:,2),'k');
plot(u(idk),res.trap(idk,2),'om');
plot(u,res.gaus(:,2),'b');
plot(u(idk),res.gaus(idk,2),'om');
%hold on;
plot(u,res.lapl(:,2),'r');
plot(u(idk),res.lapl(idk,2),'om');
plot(u,res.simp(:,2),'g');
plot(u(idk),res.simp(idk,2),'om');
plot(u,res.loba(:,2),'c');
plot(u(idk),res.loba(idk,2),'om');
title('Double integral');
xlabel('y');
ylabel('Double Integral of Cov=NegExpo')
hold off;
subplot(2,maxTry/2+1,2);
hold on;
plot(u,abs(res.trap(:,2)-res.simp(:,2)),'g');
plot(u(idk),abs(res.trap(idk,2)-res.simp(idk,2)),'om');
plot(u,abs(res.trap(:,2)-res.loba(:,2)),'c');
plot(u(idk),abs(res.trap(idk,2)-res.loba(idk,2)),'om');
%plot(u,abs(res.trap(:,2)-res.lapl(:,2)),'r');
%plot(u(idk),abs(res.trap(idk,2)-res.lapl(idk,2)),'om');
plot(u,abs(res.trap(:,2)-res.gaus(:,2)),'b');
plot(u(idk),abs(res.trap(idk,2)-res.gaus(idk,2)),'om');
title('Abs Er');
xlabel('y');
ylabel('Absolute Error |Trap-Othr.Method|')
hold off;

subplot(2,maxTry/2+1,ii+2);
hold on;
plot(u,abs(res.trap(:,2)-res.simp(:,2)),'g');
plot(u(idk),abs(res.trap(idk,2)-res.simp(idk,2)),'om');
plot(u,abs(res.trap(:,2)-res.loba(:,2)),'c');
plot(u(idk),abs(res.trap(idk,2)-res.loba(idk,2)),'om');
%plot(u,abs(res.trap(:,2)-res.lapl(:,2)),'r');
%plot(u(idk),abs(res.trap(idk,2)-res.lapl(idk,2)),'om');
plot(u,abs(res.trap(:,2)-res.gaus(:,2)),'b');
plot(u(idk),abs(res.trap(idk,2)-res.gaus(idk,2)),'om');
title('Abs Er');
xlabel('y');
ylabel('Absolute Error |Trap-Othr.Method|')
hold off;

end


figure_name = ['test_integral' num2str(jj) '.pdf'];
print('-painters', '-dpdf', '-r600', figure_name);
end