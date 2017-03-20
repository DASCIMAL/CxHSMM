% function [xi] = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi)
% Assume forward/backward passes have been done.
% This function returns the one-a-half TBN marginals:
%     xi(i,j,n,m,e,t) = Pr(x_{t-1}=i,x_t=j,m_{t-1}=n,m_t=m,e_{t-1}=e|Y)
%
% If call with parameters
% function [xi] = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi,lower,upper)
% then xi is computed only for t in the range [lower,upper]
% This feature is added to minimize memory usage
%
function [xi] = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi,lower,upper)
% 
% Last updated: 
%           Hung Bui   14/10/05, modify to minimize memory usage
%           DQ Phung   26/09/05
%           VT Duong   12/05/05	


if (nargin~=8)&&(nargin~=10)
    fprintf('Wrong number of arguments in fwback_cxhsmm\n');
    return;
elseif nargin==8
    % by default, lower=1 and upper=T
    lower = 1;
    upper = size(H,2); % observation length
end

Q = size(H,1);     % state space
T = size(H,2); % number of elements of xi to be returned
M = size(D,1);     % number of Coxian phases

BOOL = 2; 
TRUE = 1;
FALSE = 2;

xiLength = upper-lower+1;
% allocate memory, xi(i,j,n,m,e,t) = Pr(x{t-1:t}=[i,j],m{t-1:t}=[n,m],e{t-1}=e|Y)
xi = zeros(Q,Q,M,M,BOOL,xiLength);
fprintf('2TBN size (in Mbytes): %d\n',round(4*prod(size(xi))/1000000));



% note that at t = T: xi = 0
% we use following convension:
%    state variable:  i -> j
%    phase variable:  n -> m 
%    end variable:    k -> l 

for t=lower:min(upper,T-1)
    tmpIndex = t - lower + 1;
    
	for i = 1:1:Q

		% current ending status is TRUE, this happens only when n = M 
   	k = TRUE; 
		n = M;
		
		for m = 1:1:M,
			tmpSum = 0;
			clqM = P(m,:); % clqM(m,j)
			clqX = A(i,:); % A(i,j)

			if (m == M)
				l = TRUE;
				clear clqE;

				clqE(1,:) = D(m,:); % clqE(m,j) 
				tran = clqM .* clqX .* clqE;
				tran = tran';
				tmpSum = tmpSum + tran .* H(:,t+1) .* beta(:,m,l,t+1) .* alpha(i,n,k,t);
		
			end;

			l = FALSE; 
			if (m == M)
				clqE(1,:) = 1 - D(m,:); %D(m,j)
			else
				clqE = 1;
			end;

			tran = clqM .* clqX .* clqE; 
			tran = tran';
			tmpSum = tmpSum + tran .* H(:,t+1) .* beta(:,m,l,t+1) .* alpha(i,n,k,t);

			xi(i,:,n,m,k,tmpIndex) = tmpSum ./ phi(t+1); 

	 	end; % m 

		% current ending status is FALSE
		k = FALSE; 
		for n = 1:1:M,

			j = i;
			clqX = 1;

			for m = n:1:min(n+1,M),
				tmpSum = 0;
				if (n == M)
					clqM = (m==M);
		 		else
					clqM = (n==m) * (1 - D(n,i)) + (n+1==m) * D(n,i);
				end;

				l = TRUE; 
				clqE = (m==M) * D(m,j);
				tran = clqM * clqX * clqE;
				tmpSum = tmpSum + tran * H(j,t+1) * beta(j,m,l,t+1) * alpha(i,n,k,t);

				l = FALSE; 
				clqE = (m==M) * (1 - D(m,j)) + (m~=M);
                
				tran = clqM * clqX * clqE;
				tmpSum = tmpSum + tran * H(j,t+1) * beta(j,m,l,t+1) * alpha(i,n,k,t);  

				% update xi
				xi(i,j,n,m,k,tmpIndex) = tmpSum  / phi(t+1);
			end % m
		end % n 
	end % i 
end % loop over tmpIndex
 	
