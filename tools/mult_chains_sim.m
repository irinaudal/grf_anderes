function [MAINchain]=mult_chains_sim(Nchains, niter, burnin, saveimage,storage, s_obs, x_obs, x_gridM, z_gridM, Pz, z_grid_Pz,... 
                            sd_e, mu_f, S11, nu,sig2,rho, cov_type, verbose,doplots,res_dir,...
                            fixed_z,zS,zPM,zMode,zMed,fixed_f,f_gridM,a)
%%% Generate multiple chains, return as a cell array                       
                        
    MAINchain = cell(Nchains);
    for jj=1:Nchains
        %run Gibbs sampler
        disp(['Running a simple MCMC of chain: ',num2str(jj),'...']);
        [ resmc ] = analyze_grf(niter, s_obs, x_obs, x_gridM, z_gridM, Pz, z_grid_Pz,... 
                            sd_e, mu_f, S11, nu,sig2,rho, cov_type, verbose,doplots,res_dir,...
                            fixed_z,zS,zPM,zMode,zMed,fixed_f,f_gridM,a);
        disp(['All done with MCMC of chain: ',num2str(jj),'. :)']);
               
        %remove burnin
        if(burnin>0)
            resmc(1:burnin,:) = [];
        end
        %store the mcmc results to file
        if(saveimage)
            save_chain_to_file(res_dir,resmc,f_gridM,storage,jj);
        end
            

        %store all chains
        MAINchain{jj} = resmc;
    end
end %END simulating MAINchain of MCMC draws
