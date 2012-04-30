%%% Coupling from the past: Propp-Wilson Simulation
function [x] = proppwilson(N,u,p,T,verbose)

%%%
% Parameters:
%    N   = number of observations desired from the stationary distribution
%    set = set of possible stating states
%    P   = transition matrix probabilities
%    verbose = (T/F) for printing dialog of program
%%%

%%%%%%%% Begin function:
%%% Find all possible combinations of starting states
%combos = fullfact([n N]);
combos = ff2n(N)+1;
len = size(combos,1);

% prepare P matrix of correct size
if(all(size(p)==[N,N]))
    P=p;
else
%     P = zeros(N,N);
%     for j=1:(length(p))
%         P(j,:) = [p(j) , ones(1,N-1)*(1-p(j))/(N-1)];
%     end
     P = eye(N)*p(1) + (ones(N,N)-eye(N))*p(2)/(N-1);
end

disp(P);

%%% For each configuration, determine the result k steps forward
%%% If the result is the same state for all N values, stop, return that
%%% configuration.

disp(combos);

%%% Check each configuration
for t=1:T
    for j=2:(len-1)
        x = combos(j,:)';
        temp = u*(P^t);
        temp=temp(:);
        
        xnew = binornd(1,temp) + 1;
        
        if(verbose)
            disp([temp, xnew]);
        end

        % Check curent configuration
        if (all(xnew==1) || all(xnew==2) )
            if(verbose)
                disp('Yey - we found the solution!')
                disp(t);
            end
            return;
        else
            x = 0; %set faulty return to nothing
        end
    end
end


end
