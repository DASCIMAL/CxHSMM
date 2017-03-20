clear all;

% number of states
Q = 5;
% observed alphabet spaces, 
% V(k) is domain of kth observed element(feature)
V = [Q];
%V = [Q Q Q];

% number of phases, set to 1 if using HMM
M = 3;

% number of independent "features" in each observation
K = length(V);

% generate some observation sequences
N = 1;    % number of iid observation sequences
T = 16;   % length of each sequence

NULL = -1; % used in the null trick

for (n=1:N)
	mobsq{n} = zeros(K+1,T);
	for k=1:K
		clear uniDist;
		uniDist = (1/(V(k)+1)) * ones(1,V(k)+1);

		% assume there is no missing ovservation for now
		uniDist = (1/(V(k))) * ones(1,V(k));
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
        tmp = sample_discrete(uniDist);
        if tmp <= Q 
            mobsq{n}(K+1,t)=tmp;
        else
            mobsq{n}(K+1,t)= NULL;
        end

		  % make state unobserved for now
		  mobsq{n}(K+1,t) = NULL;
    end
end

% just pick one for experiment now
obsq = mobsq{1};


% initialize parameters, see init_cxhsmm for further information
[PI,A,P,D,B] = init_cxhsmm(Q,M,K,V,'uniform');

% make what we osbserve is x% chance of being the state
for (k=1:K)
	B{k} = ones(Q,Q);
	for i=1:Q
		B{k}(i,i) = 1000;
	end
	B{k} = mk_stochastic(B{k});
end

% compute the observation probability for obsq
% H(i,t) = Pr(y_t | x_t =i)
[H lnH] = compute_obprob(B,obsq);

% compute forward
% [alpha, phi] = forward_cxhsmm(PI,A,P,D,H);

% compute backward
% [beta] = backward_cxhsmm(PI,A,P,D,H,phi);
	
% compute 1/2TBN marginals
% [xi] = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi);

% testing viteri algorithm
[vtbsq prob] = viterbi_cxhsmm(PI,A,P,D,H);   % use normal scale
[lvtbsq lprob] = viterbi_cxhsmm(PI,A,P,D,H,'uselog');

% testing viterbi deconding based on smoothing distribution 
[smdsq] = smdecode_cxhsmm(PI,A,P,D,H);
