
function [smdsq] = smdecode_cxhsmm(PI,A,P,D,H)

% performs the best sequence viterbi decoding using smoothing distribution
% It returns:  smdsq(t) is the best state at time t given the entire observation sequence
%              (as opposed to the normal Viterbi, which conditions only observation up to time t only)
%                     i.e., smdsq(t) = argmax_Xt P(Xt | Y)


Q = size(H,1);     % state space
T = size(H,2);     % observation length  
M = size(D,1);     % number of Coxian phases
BOOL = 2;

% compute forward
[alpha, phi] = forward_cxhsmm(PI,A,P,D,H);

% compute backward
[beta] = backward_cxhsmm(PI,A,P,D,H,phi);
	
% compute 1/2TBN marginals
[xi] = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi);

% allocate memory
smdsq = zeros(3,T);  % 3 here means place holders for x_t, m_t, e_t

for (t=1:T)
	% compute one-slice smoothing at time t
	if (t < T)
		tmp = xi(:,:,:,:,:,t);
   	onesm = reshape(sum(sum(tmp,3),2),[Q M BOOL]);  % 2 means TRUE or FALSE
		a = reshape(onesm,[1 Q * M * BOOL]);   % first expand to get maximum value 
		[b indx] = max(a);                  %
		[bi bn bk] = ind2sub(size(onesm),indx); % then map back to get its indices

		% get the state that returns best prob
		smdsq(:,t) = [bi bn bk]';
	else
		tmp = xi(:,:,:,:,:,t);
   	onesm = reshape(sum(sum(sum(tmp,5),4),1),[Q M]);  % 2 means TRUE or FALSE
		a = reshape(onesm,[1 Q * M]);        % first expand to get maximum value 
		[b indx] = max(a);                   %
		[bi bn] = ind2sub(size(onesm),indx); % then map back to get its indices

		% get the state that returns best prob
		smdsq(:,t) = [bi bn 1]';  % assume the last ending status is TRUE=1
	end
end
