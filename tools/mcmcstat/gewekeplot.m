function zz=gewekeplot(chain,res_dir)
%GEWEKEPLOT Plot Geweke's diagnostic
% gewekeplot(chain) plots Geweke's diagnostic for increasing number of 
% iterations. See geweke.m

% ML, 2002
% $Revision: 1.3 $  $Date: 2003/05/07 12:22:19 $

if(nargin<2)
    res_dir = false;
end

[nsimu,npar]=size(chain);

n  = 40;
e  = fix(nsimu/2);
l  = fix(e/n);
ii = 1:l:e;

z = zeros(length(ii),npar);
for i=1:length(ii)
  z(i,:) = geweke(chain(ii(i):end,:));
end

it = 0;
maxRow = 6;
if (npar<maxRow)
    maxRow = npar;
end
for i=1:npar
  if(gcd(i-1,maxRow)==maxRow)
      it = it+1;
      fig=figure();
  end
  %take subset of ids
  sub = i - (it-1)*maxRow;
  h = subplot(maxRow,1,sub);
  
  plot(ii,z(:,i));
  set(h,'XLim',[1 e]);
  ylabel(['par.',num2str(i)]);
%   if i==1
  if sub==1 
    title('Geweke diagnostics');
  end
%  if i==npar
  if sub==maxRow
    xlabel('first iteration used');
  end
  if (sub==maxRow || i==npar)
      if(res_dir)
          screensz = get(0, 'MonitorPositions');
          set(fig, 'Position', [0 0 screensz(3)/2 screensz(4) ] );
          print(fig,'-djpeg','-r300',[res_dir,'/geweke_',num2str(it),'.jpg']);
          close(fig);
      end
  end
end  

if nargout > 0
  zz=z;
end
