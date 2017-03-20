function [beta] = backward_cxhsmm(PI,A,P,D,H,phi)

% function [beta] = backward_cxhsmm(PI,A,P,D,H,phi)
% It returns:
%   	beta(i,m,k,t) = Pr(y_t+1:T | x_t=i,m_t=m,e_t = k) / scaling factor
% Last updated: DQ PHUNG 26/09/2005 
%               VT Duong 12/05/2005 

Q = size(H,1);     % state space
T = size(H,2);     % observation length  
M = size(D,1);     % number of Coxian phases

% boolean values definitions
BOOL = 2;
TRUE = 1;
FALSE = 2;

% memory allocation
beta = zeros(Q,M,BOOL,T);

% notations:
%     e: e node
%

% --------------------------------------------------------------------
% This is a backward recursive procedure, starting at time t = T
% t = T
% --------------------------------------------------------------------
t = T;
beta(:,:,:,T) = ones(Q,M,BOOL);


% --------------------------------------------------------------------
% Starting the recursive procedure 
% t = T-1:1:1 
% --------------------------------------------------------------------
% NOTE: clqM; clqX; do not depend on eps_{t+1},e_{t+1}

for t = T-1:-1:1,
	for i = 1:1:Q, 
											
   	clear tmpSum;
		tmpSum = 0;
      % ----------------------------------------------------------------------------			  
		% CASE 1:  e = TRUE, 
		% happens only when the last phase has been reached
      % ------------------------------------------------------------------------------			  
		e = TRUE; 
		n = M;

		for nx_n = 1:1:M,
			% Cliques: clqM, clqX 
			clqM = P(nx_n,:); %P(nx_n,nx_i)
			clqX = A(i,:);		%A(i,nx_i)


			if ( nx_n == M)
				% -------------------------------------------
				% Case 1.1: nx_n = M 
				% -------------------------------------------
				nx_e = TRUE; 

				% Cliques:  clqE
				clear clqE;

				clqE = D(nx_n,:); % clqE(nx_n,nx_i) 

				clear tran;
				tran = clqM .* clqX .* clqE;
				tran = tran';
				tran = tran .* H(:,t+1) .* beta(:,nx_n,nx_e,t+1);
				tmpSum = tmpSum + sum(tran);
			end;

			% -------------------------------------------
			% Case 1.2: nx_e = FALSE 
			% -------------------------------------------
			nx_e = FALSE; 

			% Cliques:  clqE
		 	clear clqE;

			if ( nx_n ==  M)
				clqE = 1 - D(nx_n,:);
		 	else
				clqE = 1;
			end;

			clear tran;
			tran = clqM .* clqX .* clqE;
			tran = tran';
			tran = tran .* H(:,t+1) .* beta(:,nx_n,nx_e,t+1);
			tmpSum = tmpSum + sum(tran);

	 	end; % nx_n = 1:1:M

		beta(i,n,e,t) = tmpSum / phi(t+1);
		% END CASE 1


	
	   clear tmpSum;
		tmpSum = 0;
     	% ----------------------------------------------------------------------------			  
		% CASE 2:  e = FALSE, 
      % ------------------------------------------------------------------------------			  
		e = FALSE; 
		for n = 1:1:M,

			tmpSum = 0;
         
			nx_i = i;
				
			clqX = 1;

			for nx_n = n:1:min(n+1,M),
				if ( n == M )
					clqM = (nx_n==M);
		 		else
					clqM = (n==nx_n) * (1 - D(n,i)) + (n+1==nx_n) * D(n,i);
				end;

					
				% -------------------------------------------
				% Case 3.1 
				% -------------------------------------------
				nx_e = TRUE; 
				clqE = (nx_n==M) * D(nx_n,nx_i);
				tran = clqM * clqX * clqE;
				tmpSum = tmpSum + tran * H(nx_i,t+1) * beta(nx_i,nx_n,nx_e,t+1);

				% -------------------------------------------
				% Case 3.2 
				% -------------------------------------------
				nx_e = FALSE; 
				clqE = (nx_n==M) * (1 - D(nx_n,nx_i)) + (nx_n~=M); 
				tran = clqM * clqX * clqE;
				tmpSum = tmpSum + tran * H(nx_i,t+1) * beta(nx_i,nx_n,nx_e,t+1);
			end; % nx_n

			beta(i,n,e,t) = tmpSum / phi(t+1);
		end; % n 
		% END CASE 2 

	end; % i 
end; % t = T-1:1:1

