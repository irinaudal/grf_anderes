function []=save_chain_to_file(res_dir,resmc,f_gridM,storage,jj)
% Store all current results of the MCMC chain to a file
if (nargin<5)
    jj=false;
end
    if (jj)
        save([res_dir,'/mcmcdraws',num2str(jj),'.txt'], 'resmc', '-ASCII');
    else
        save([res_dir,'/mcmcdraws.txt'], 'resmc', '-ASCII');
    end
    save([res_dir,'/data.txt'], 'storage', '-ASCII');
    save([res_dir,'/true_f.txt'],'f_gridM','-ASCII');
    clear('storage');
end
