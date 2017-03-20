function [alpha,phi] = forward_cxhsmm(PI,A,P,D,H)

% function [alpha,phi] = forward_cxhsmm(PI,A,P,D,H)
% This function performs the forward pass. It returns:
%      alpha(i,m,k,t) = Pr(x_t=i,m_t=m,e_t=k|y_1:t)
%      phi(t) = Pr(y_t | y_1:t-1)
% Last updated: DQ Phung 26/09/2005 
%					 VT Duong 12/05/2005

Q = size(H,1);     % state space
T = size(H,2);     % observation length  
M = size(D,1);     % number of Coxian phases

% boolean values definitions
BOOL = 2; 
TRUE = 1;
FALSE = 2;

alpha = zeros(Q,M,BOOL,T);  % scaled alpha
palpha = zeros(Q,M,BOOL,T); %  partially scaled alpha:
phi = zeros(1,T);           % scaling factor

% at first time slice
t = 1;
for i = 1:1:Q,     % loop over state variable
	for m = 1:1:M,  % loop over phase variable, end variable is considered separately

		clqM = P(m,i);
		clqX = PI(i);
		  
		k = TRUE; 
		clqE = (m==M) * D(m,i);
		palpha(i,m,k,t) = H(i,t) *  clqE * clqM * clqX;

		k = FALSE; 
		clqE = 1 - (m==M) * D(m,i);
		palpha(i,m,k,t) = H(i,t) * clqE * clqM * clqX;
	end; % m  
end;  % i 

phi(t) = sum(sum(sum(palpha(:,:,:,t))));
% prevent the case phi(t) is too small
if (phi(t) < eps) 
	phi(t) = eps; 
	warning('phi(t) is too small in forward_cxhsmm, replaced by eps');
end;
alpha(:,:,:,t) = palpha(:,:,:,t) ./ phi(t); 
		

% at other time slice, the convesion is
% state node:  j -> i
% phase node:  n -> m
% end node:    l -> k
for t = 2:1:T,
	for i = 1:1:Q, 
		for m = 1:1:M,
			
			tmpSum1 = 0;
			tmpSum2 = 0;
		
			% k = TRUE, this case only happens when m = M 
			clqE1 = (m==M) * D(m,i);
				
			% k = FALSE;
			clqE2 = (m~=M) + (m==M)*(1 - D(m,i));

			% loop through time slice previous slice: 
			% previous ending status is TRUE
			l = TRUE; 
			n = M;

			clqM = P(m,i);
			clqX = A(:,i);

			tmp = clqM .* clqX .* H(i,t);
			tmp = tmp .* alpha(:,n,l,t-1); 
			
			tran1 = tmp .* clqE1; 
			tran2 = tmp .* clqE2; 

			tmpSum1 = tmpSum1 + sum(tran1);  
			tmpSum2 = tmpSum2 + sum(tran2);  

			% previous ending status is FALSE
			l = FALSE; 
			j = i;
			clqX = 1;	

			for n = max(1,m-1):1:m,
				if (n == M)
					clqM = 1;
		 		else
					clqM = (n==m) * (1 - D(n,j)) + (n==m-1) * D(n,j); 
		 		end;

				tmp = clqM * clqX * H(i,t) * alpha(j,n,l,t-1);
				
				tran1 = clqE1 * tmp; 
				tran2 = clqE2 * tmp; 

				tmpSum1 = tmpSum1 + tran1;
				tmpSum2 = tmpSum2 + tran2;
		 	end

			% update partialled scaled alpha 
			palpha(i,m,TRUE,t) = tmpSum1;		
			palpha(i,m,FALSE,t) = tmpSum2;		

		end; % n 
	end; % i 

	phi(t) = sum(sum(sum(palpha(:,:,:,t))));
	% prevent division by zero
	if (phi(t) < eps)
		warning('phi(t) is too small in forward_cxhsmm, replaced by eps');
		phi(t) = eps;
	end;

	% now compute the scaled alpha
	alpha(:,:,:,t) = palpha(:,:,:,t) ./ phi(t); 
end % loop over t = 2:1:T
