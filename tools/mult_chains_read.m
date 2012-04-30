function [MAINchain] = mult_chains_read(Nchains,res_dir)
      
% Read multiple chains from folder, return cell array
    MAINchain = cell(Nchains);

    for jj=1:Nchains
        %run Gibbs sampler
        disp('Reading data result for MCMC test...');
        [ resmc ] = load([res_dir,'/mcmcdraws',num2str(jj),'.txt']);
        
        %store all chains
        MAINchain{jj} = resmc;
    end
end %END reading MAINchain of MCMC draws from previously stores file
