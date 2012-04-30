function [id_p,Pz,zgrid,zs,zPM,zMode,zMed,zRand,zVar] = posterior_redshift_access(N,Posterior_filename,BPZ_filename,res_dir,verbose,sd_z,type,zgrid) 
    
%####################################################################################################
% posterior_redshift_access   Function accesses the posterior file of redshift and
%                             computes posterior statistics
%
% Input: N    = Number of galaxies to read from file
%        Posterior_filename = filename where redshift posterior is evaluated       
%        BPZ_filename = filename where spectroscopic redshift is stored (pair file to one above)       
%        verbose = (T/F) display progress of program
%        zgrid = Grid of redshift for posterior       
%
% Output: Pz    = Posterior of redhsift
%         zgrid = Grid of redshift for posterior 
%         zs    = Spectroscopic redshift 
%         zPM   = Posterior mean
%         zMode = Posterior mode
%         zMed  = Posterior median
%         zRand = Random sample of posterior dedshift
%         zVar  = Posterior variance
%####################################################################################################

if(verbose)
    disp('Reading Redshift-Posterior file, Reasing BPZ file...');
end

if(nargin < 8 )   %grid values are not defined for posterior file
    zmin = 0.005;    %redshift minimum - according to Posterior file       
    zmax = 4.005;    %redshift maximum
    dz = 0.01;       %redshift increment
    zgrid = linspace(zmin,zmax,401)';    %Redshift grid for Posterior
end

%read Spectroscopic Redshift file
[~, bpz] = hdrload(BPZ_filename); 
%read posterior file
[~, Pz] = hdrload(Posterior_filename);    %[header Pz] = hdrload('mockPosterior1.txt');
id_p = randsample(length(Pz),N);       %select a random sample of posteriors of size Ng

   %new
   zs = bpz(:,10);
   id1 = find(zs>3);        id1 = id1(randsample(length(id1),N/4));
   id2 = find(zs>2 & zs<3); id2 = id2(randsample(length(id2),N/4));
   id3 = find(zs>1 & zs<2); id3 = id3(randsample(length(id3),N/4));
   id4 = find(zs>0 & zs<1); id4 = id4(randsample(length(id4),N/4));
   id = [ id1; id2; id3; id4 ];
   id_p = id(randsample(N,N));
   
zs = bpz(id_p,10);
Pz = Pz(id_p,:);        
Pz(:,1) = [];            %remove ID column of Posterior  
for i=1:N
    %normalize the discretized density
    Pz(i,:) = Pz(i,:)/sum(Pz(i,:));
end
    
if (strcmp(type,'normal'))
    %%%%%%%
    %Change the prior to discretized normal density:  z ~ N(mu_z, sd_z)
    %sample mean mu_z ~ U(zmin,zmax)     %(  suitable posterior is: N(mu=2,sd=0.775)  )
    mu_z = zmin + (zmax-zmin).*rand(N,1);
    %compute normal distributions with various means mu_z
    for i=1:N
        %compute normal density over grid with mu_z(i)
        Pz(i,:) = normpdf(zgrid,mu_z(i),sd_z);
        %normalize the discretized density
        Pz(i,:) = Pz(i,:)/sum(Pz(i,:));
        %generate zpectroscopic redshift as a random draw from the same prior
        zs(i) = round(10000*normrnd(mu_z(i),0.1))/10000;
    end
    %%%%%%%
end

%%%%%%%
% Plot a small sample of posteriors
fig = figure();
plot(zgrid,Pz(randsample(N,50),:));
xlabel('z'); ylabel('Pz');

title(['20 Prior distributions: N(zS,',num2str(sd_z),'^2)']);

export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/posterior.pdf'])
close(fig);
%%%%%%%

Nz = size(Pz,2);         %Number of data points in Posterior
%%%Calculate Posterior Mean
zPM = sum(Pz.*zgrid(:,ones(N,1))', 2); %sum over all columns (2) within each row
%%%Calculate Posterior Variance
zVar = sum((zgrid(:,ones(N,1))' - zPM(:,ones(1,Nz))).^2.*Pz, 2);
%%%Calculate Posterior Mode
zMode = zeros(N,1);       
for j=1:N
    zMode(j) = min(zgrid( Pz(j,:) == max(Pz(j,:)) ));  %select Ties by min()        
end
%%%Calculate Posterior Median
zMed = zeros(N,1);       
for j=1:N
    tmp = cumsum(Pz(j,:));
    idx = [find(tmp <= 0.5,1,'last'), find(tmp >= 0.5,1,'first')];
    if(length(idx)==1)
        idx = [idx idx];
    end
    zMed(j) = zgrid(idx(1)) + 0.5*(zgrid(idx(2)) - zgrid(idx(1)));         
end
%%%Sample from Posterior Distr
id = rand(N,1);
zRand = zeros(N,1); 
for j = 1:N
    zRand(j) = zgrid(find(cumsum(Pz(j,:))>=id(j),1));
end