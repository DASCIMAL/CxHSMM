%
%[PI,A,P,D,B,loglik]=em_cxhsmm(Q,M,K,V,mobsq,numIter,threshold,...
%                              obsStateStr)
%
% The last argument is optional
% If they are supplied, then obsStateStr must be 'observe_state'
% and the last "feature" in each obs sequence is 
% treated as the observed state sequence
% (-1 is used to indicate no-observation)
%
% Last updated: Hung Bui 06/10/2005
%               DQ Phung 26/09/2005
%               VT Duong 12/05/2005

function [PI,A,P,D,B,loglik] = em_cxhsmm(Q,M,K,V,mobsq,numIter,threshold,obsStateStr)


clear PI; clear A; clear P; clear D; clear B; clear loglik; 

if (nargin~=7)&&(nargin~=8)
    fprintf('Wrong number of arguments in em_cxhsmm\n');
    exit;
elseif nargin==7
    % by default, state is not observed
    stateObserved= 0;
elseif strcmp(obsStateStr,'observe_state')
    stateObserved = 1;
else
    fprintf('Something wrong with observe_state argument in em_cxhsmm.\n');
    exit;
end
    
N = length(mobsq); 

% add a small prior to prevent division by zero
% alternatively, we need to do zero-checking before doing division.
dirichlet_priors = 1e-100;

%  initialize the parameters, we can choose 'random' or 'uniform'
%  see inside this routine for further information.

if stateObserved
    % artifically increase the number of features
    K = K+1;
    % the last feature is the state observation
    % so its domain size is Q (not including the null value)
    V = [V Q];
end

[PI,A,P,D,B] = init_cxhsmm(Q,M,K,V,'uniform'); 

if stateObserved
    % we have to modify the last component of B here
    % use a diagonal matrix as the observation model
    B{K}=eye(Q);
end
    

% first iteration
iter = 1;
ctime = cputime;  % current time stamp
[essPI,essA,essP,essD,essE,essB,loglik(iter)] = compute_ess_cxhsmm(PI,A,P,D,B,V,mobsq);

% get elapsed time
ctime = cputime - ctime;

ploglik = -inf;          % previous log likelihood
cloglik = loglik(iter);  % current log likelihood
disp(sprintf('\t + loglik = %.20f, iteration: %d -- took %.2f (s)',cloglik,iter,ctime));

while ((abs(cloglik - ploglik) > threshold) & (iter < numIter))
	if (iter > 1) % the first E-step has been done 
		ctime = cputime; 

		% comput expected sufficient statistics given newly updated parameters
		[essPI,essA,essP,essD,essE,essB,loglik(iter)] = compute_ess_cxhsmm(PI,A,P,D,B,V,mobsq);

		ctime = cputime - ctime;
		ploglik = loglik(iter-1);
		cloglik = loglik(iter);
		disp(sprintf('\t + loglik = %.20f, iteration: %d -- took %.2f (s)',cloglik,iter,ctime));
	end

	% reestimating initial state probability PI
	essPI = essPI + dirichlet_priors;
	PI = normalize(essPI);
	
	
	% reestimating transition probability A
	essA = essA + dirichlet_priors;
	% don't add dirichlet_priors for diagonal entries
	A = A .* (1 - eye(Q,Q));
	A = mk_stochastic(A); 

	
	% re-estimating phase initial distribution P
	essP = essP + dirichlet_priors;
	P = normalize(essP,1);

	% re-estimating duration model D
	essD = essD + dirichlet_priors;
	essE = essE + dirichlet_priors;
	D(:,:) = essD(:,:) ./ essE(:,:); 

	% re-estimating observation model
    
    % do not reestimate the last component
    % if state is observed
    if stateObserved
        tmpK = K-1;
    else
        tmpK = K;
    end
    
	for k=1:tmpK
		essB{k} = essB{k} + dirichlet_priors;
		B{k} = mk_stochastic(essB{k});	
	end

	% increase the loop count
	iter = iter + 1;
end

if stateObserved
    % do not return the last component of B
    tmpB=B(1:K-1);
    B=tmpB;
end
