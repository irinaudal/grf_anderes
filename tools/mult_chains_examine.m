function [rchain schain pchain saveR1 saveR2 saveR3] = mult_chains_examine(MAINchain, mu_f,z_obs,f_gridM, x_gridM,z_gridM, fixed_z,fixed_f, res_dir)
      
    if (~iscell(MAINchain))
        error('MAINchain must be a cell array.');
    end
      
    Nchains = length(MAINchain);
    N = length(z_obs);
    M = numel(f_gridM);
    
    rchain  = cell(Nchains);
    schain  = cell(Nchains);

      
    for jj=1:Nchains
        
        %load chain
        resmc = MAINchain{jj};
%        % Restrict the range of iterations of MCMC draws
%        resmc = resmc(examine_range,:);
                        
        % compute statistics
        disp('Checking MCMC draws...');
        z_true = z_obs'; 
        f_true = f_gridM(:)';
%    [delta delta_true d c P P_true xsi stat_foci foci_true stat_ex ex_true stats stats_true] = diagnose_grf(resmc, ...
        [stats stats_true statsP stats_trueP] = eval_stats(resmc,jj,mu_f, ...
            z_true, f_true, x_gridM(1,:),z_gridM(:,1), res_dir);
        
        %store all chains
        rchain{jj}    = resmc;
        schain{jj}    = stats;
        pchain{jj}    = statsP;

    end
    
    
    %Make Diagnostic plots (script)
    gelman_diag
    
%     id = 1:36;
%     if(fixed_z==false)
%         gew = gewekeplot(resmc(:,id),res_dir);
%     else    
%         gew = gewekeplot(resmc(:,(N)+(id)),res_dir);
%     end
%     
%     gew = gewekeplot(stats,res_dir);
    
end  %END examine MCMC draws: resmc and their statistics
