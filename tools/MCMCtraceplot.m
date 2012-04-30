function []=MCMCtraceplot(chain,res_dir)%,a,s,e)

if(nargin<2)
    res_dir = false;
end

[nsimu npar] = size(chain);

it = 0;
maxRow = 4;
if (npar<maxRow)
    maxRow = npar;
end

for k=1:npar
    if(gcd(k-1,maxRow)==maxRow)
        it = it+1;
        fig=figure();
        trackK = k;
    end
    %take subset of ids
    sub = (1:3)  + (k-trackK)*3;

    %%% select the range of points to plot, then:
    X = chain(:,k);
    
    % 1. make traceplot
    subplot(maxRow,3,sub(1));
    plot(X);
    xlabel('iter, t');
    ylabel('MCMC draws');
    title(['Traceplot ',num2str(k)]);

    % 2. make hist and fitted curve to that hist
    subplot(maxRow,3,sub(2));
    [Ct Bin]=hist(X,20);
    dx = Bin(2)-Bin(1);
    bar(Bin,Ct,'FaceColor',[.8 .8 1])
    hold on;
    [F,XI]=ksdensity(X);
    plot(XI,F*nsimu*dx,'--r','LineWidth',2);
    hold off;
    xlabel('MCMC draws');
    ylabel('Freq');
    title(['Histogram ',num2str(k)]);

    % 3. make ACF plot
    subplot(maxRow,3,sub(3));
    autocorr(X,99);
    xlabel('Lag');
    ylabel('ACF');
    title(['ACF ',num2str(k)]);

    if (sub(3)==maxRow*3 || k==npar)
        if(res_dir)
            screensz = get(0, 'MonitorPositions');
            set(fig, 'Position', [0 0 screensz(3) screensz(4) ] );
            export_fig(fig,'-painters', '-r600', '-q101', [res_dir,'/mcmc_trace_',num2str(k),'_meanALL.pdf'])
            close(fig);
        end
    end
end        
