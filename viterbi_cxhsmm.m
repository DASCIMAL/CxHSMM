function [vtbsq, prob] = viterbi_cxhsmm(PI,A,P,D,H,logStr);

% function [vtbsq] = viterbi_cxhsmm(PI,A,P,D,H)
% This function perform the viterbi decoding. It returns:
%        'vtbsq' is the best sequence decoded ['vtbsq' is short for 'viterbi sequence']
%        'prob"    is the probability for this decoded sequence, log(prob) if operating in log space
%
% The viterbi algorithm works as follows: (here x_t is a product state at time t) 
% We keep track of the best sequence upto time t via the variable 'del'
%     del(t,i) = max_{x_1:t-1} Pr(x_1:t-1, x_t = i, y_1:t)
%     del(t+1,i) = [max_j del(t,j) A(j->i)] * Pr(y_t+1 | i) 
% We need phi(t,i) to keep the value of j in above recursion that returns the maximization.
% When all is done, we perform a backtracking as follows:
%     prob = max_i del(T,i)
%     vtbsq(T) = argmax_i del(T,i)
%     vtbsq(t-1) = phi(t,vtbsq(t))  for t = T-1, ..., 1
%
% Note: The code is not written to optimize for speed but rather for [somewhat] clarity.
%       We can operate on logarithmic space to avoid numerical unstable. Use 'uselog' option
%       to enable this operation.
% 
% Last updated:   DQ. Phung  07/10/2005 
%                 HH. Bui    30/09/2005 
%                 VT. Duong  17/05/2005 

if (nargin == 5)
	useLogFlag = 0;    % operate on normal probability number 
elseif ((nargin == 6) && (strcmp(logStr,'uselog'))) 
	useLogFlag = 1;    % operate on logarithmic space
else
    fprintf('Something wrong in uselog argument in viterbi_cxhsmm\n');
    return;
end

Q = size(H,1);     % state space
T = size(H,2);     % observation length  
M = size(D,1);     % number of Coxian phases

% boolean values definitions
BOOL = 2; 
TRUE = 1;
FALSE = 2;

% allocate memory
if (useLogFlag) 
	del = -Inf * ones(Q,M,BOOL,T); % maximum probability at time t
else
	del = zeros(Q,M,BOOL,T);       % -Inf in log-domain, and zero in prob-domain   
end
phi = zeros(Q,M,BOOL,T,3);  % previous state that returns maximization


% at the first time slice t = 1
t = 1;
for i=1:Q,  % loop through state  node
	for m = 1:1:M, % loop through Coxian phase node,
			         % the end node is consider separately through k

	   % clqM, clqX, are not affected by {k,l}
		clqM = P(m,i);
		clqX = PI(i);

		k = TRUE; % Pr(e = TRUE | m,i)
		clqE = (m==M) * D(m,i);

		if (useLogFlag)
			del(i,m,k,t) = zlog(clqE) + zlog(clqM) + zlog(clqX) + zlog(H(i,t));
		else
			del(i,m,k,t) = clqE * clqM * clqX * H(i,t); 
		end
		phi(i,m,k,t,:) = 0;

		k = FALSE; % 1 - Pr(e = TRUE|m,i)
		clqE = 1 - (m==M) * D(m,i); 

		if (useLogFlag)
			del(i,m,k,t) = zlog(clqE) + zlog(clqM) + zlog(clqX) + zlog(H(i,t)); 
		else
			del(i,m,k,t) = clqE * clqM * clqX * H(i,t); 
		end
		phi(i,m,k,t,:) = 0;
	end % m  
end % i 

% at other time slices, the naming convension is this:
% 		for state node:    j  -> i
% 		for phase node:    n  -> m
%     for end node:      l  -> k
for t = 2:1:T  
	for i = 1:1:Q  			 % loop through state node
		for m = 1:1:M         % loop through Coxian phase node
			for k = 1:1:BOOL   % loop through end node

				% loop through the previous slice
				for j = 1:1:Q, 
					for n = 1:1:M,
						for l = 1:1:BOOL,

							% clqE
							clqE = (k==TRUE) * (m==M) * D(m,i) ...
							       + (k==FALSE) * (1 - (m==M) * D(m,i)); 

							% clqM;
							clqM = (l==TRUE) * P(m,i) ...
							       + (l==FALSE) * (n~=M) * ((m==n) * (1-D(n,j)) + ((n+1)==m) * D(n,j)); 

							% clqX:
							clqX = (l==TRUE) * A(j,i);
							clqX = clqX + (l==FALSE) * (j==i); 

							if (useLogFlag)
								currDel = del(j,n,l,t-1) + zlog(clqE) + zlog(clqM) + zlog(clqX) + zlog(H(i,t));
							else
								currDel = del(j,n,l,t-1) * clqE * clqM * clqX * H(i,t);
							end

							if (currDel > del(i,m,k,t))
								del(i,m,k,t) = currDel;
								phi(i,m,k,t,1) = j;
								phi(i,m,k,t,2) = n;
								phi(i,m,k,t,3) = l;
							end;
										
						end % l 
					end %  
				end % j 
			end % k 
		end % m 
	end % i
end % t


% now the backtracking 
vtbsq = zeros(3,T);

% at time T, sorry for a little trick to extract fast maximum index :-)
c = del(:,:,:,T);
a = reshape(c,[1 Q * M * BOOL]);   % first expand to get maximum value 
[b indx] = max(a);                 %
[bi bn bk] = ind2sub(size(c),indx);% then map back to get its indices

vtbsq(1,T) = bi;
vtbsq(2,T) = bn;
vtbsq(3,T) = bk;

% the probability for this sequence, it is the log(prob) if useLogFlag is set
prob = del(bi,bn,bk);  

% other time slices other than T 
for t = (T-1):-1:1,
   i = vtbsq(1,t+1);
	n = vtbsq(2,t+1);
	k = vtbsq(3,t+1);

	vtbsq(1,t) = phi(i,n,k,t+1,1);
	vtbsq(2,t) = phi(i,n,k,t+1,2);
	vtbsq(3,t) = phi(i,n,k,t+1,3);
end

% that's it
