% this script set up some artificial data to be used in demo1
% run this one first, then demo1_learn, then demo1_classify

clear all;

% number of states
Q = 10;
% observed alphabet spaces, 
% V(k) is domain of kth observed element(feature)
V = 2 * ones(1,1000);

% number of phases, set to 1 if using HMM
M = 5;

% number of independent "features" in each observation
K = length(V);

% generate some observation sequences
N = 1;     % number of idd observation sequences
T = 100;   % length of each sequence

NULL = -1; % used in the null trick

for (n=1:N)
	mobsq{n} = zeros(K+1,T);
	for k=1:K
		clear uniDist;
		uniDist = (1/(V(k)+1))*ones(1,V(k)+1);
		for (t=1:T)
            
            % construct random sequence for feature k
            tmp = sample_discrete(uniDist);
            if (tmp <= V(k))
                mobsq{n}(k,t) = tmp;
            else
                % testing the null trick
                mobsq{n}(k,t) = NULL;
            end
		end
    end
    uniDist = (1/(Q+1)) * ones(1,Q+1);
    for t=1:T
        
        % construct random state sequence
        tmp=sample_discrete(uniDist);
        if tmp <= Q 
            mobsq{n}(K+1,t)=tmp;
        else
            mobsq{n}(K+1,t)= NULL;
        end
    end
end
