% function[essPI,essA,essP,essD,essE,essB,loglik] = ess_cxhsmm(PI,A,P,D,B,V,obsq,maxXiLength)
%
% Compute the expected sufficient statistics given current parameter and one
% observation sequence. 
%
% The last parameter maxXiLength is optional. This controls the maximum length
% of the 2-TBN array xi to be computed in each call to fwdback_cxhsmm.
% If not given, its default value is length(obsq).
% This is used to miminize the memory usage. If length(obsq) is too large
% and you are getting an out of memory error, then supply a smaller value
% for maxXiLength would fix the problem, at the cost of matlab running abit
% slower due to multiple calls to the same function. Always use the largest
% possible value for maxXiLength that your computer memory can allow for.
%
function[essPI,essA,essP,essD,essE,essB,loglik] = ess_cxhsmm(PI,A,P,D,B,V,obsq,maxXiLength)

% Last updated:  Hung Bui 05/10/2005
%                DQ Phung 26/09/2005 
%					  VT Duong 12/05/2005 						
%====================================================================

% xi(z{t},z{t+1},exp{t},x{t},x{t+1},m{t},m{t+1},e{t},t)
% gamma(z{t},x{t},m{t},exp{t},e{t},t)

% Note that xi(:,:,:,:,:,:,T) = 0;

if (nargin~=7)&&(nargin~=8)
    fprintf('Wrong number of arguments in ess_cxhsmm.\n');
    return;
elseif nargin==7
    % by default, maxXiLength is the length of the observation sequence
    maxXiLength = size(obsq,2);
end

T = size(obsq,2);
Q = size(A,1);
M = size(D,1);
K = length(V);

TRUE = 1;
FALSE = 2;


% compute the observation probability ONCE at every instance of time
% H(i,t) = Pr(y_t | x_t =i)

% remove last argument if dont want to do the H scale trick
[H lnHScale] = compute_obprob(B,obsq, 'scale');

% forward pass 
[alpha,phi] = forward_cxhsmm(PI,A,P,D,H);

% backward pass
[beta] = backward_cxhsmm(PI,A,P,D,H,phi);


% loglikelihood for this obsq
loglik = sum(log(phi)) - lnHScale;

% compute the one-slice smoothing marginal
% gamma(i,m,e,t) = Pr(x_t=i,m_t=m,e_t=e | Y)
gamma = alpha .* beta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the ESS that does not require xi first

%
% expected sufficient statistics for PI
%
clear tmp; tmp = gamma(:,:,:,1);
essPI = sum(sum(tmp,3),2);

%
% essE
%
essE = zeros(M,Q);

%tmp = gamma(i,n,FALSE,1:T-1);
%essE(n,i) = essE(n,i)+sum(tmp);
clear tmp;
tmp = gamma(:,1:M-1,FALSE,1:T-1);
essE(1:M-1,:) = essE(1:M-1,:) + (sum(tmp,4))';


% tmp = gamma(i,M,:,:);
% essE(M,i) = essE(M,i)+sum(sum(tmp));
clear tmp;
tmp = gamma(:,M,:,:);
essE(M,:) = essE(M,:) + (sum(sum(tmp,4),3))';


%
%essB
%
essB = cell(1,K);
for k=1:K
    essB{k} = zeros(Q,V(k));
end

% from gamma to marginal of x
% stateMarg(i,t) = sum(sum(gamma(i,:,:,t)))
%

% Is squeeze still correct when no of state
% of no of time slice is 1?
stateMarg = squeeze(sum(sum(gamma,2),3));

for i=1:Q,
    for (k=1:K)
        for v=1:V(k)
            %clear tmpB;
            
            %tmpB(1) = 0;
            tmpB = zeros(1,T);
            %for t = 1:1:T,
                %clear tmpBt;
                %tmpBt = gamma(i,:,:,t);
                
                %%tmpB(t) = sum(sum(tmpBt)) * delta(obsq(k,t),v);
                %tmpB(t) = stateMarg(i,t) * (obsq(k,t)==v);
               
            %end;
            
            % vectorized version
            tmpB = stateMarg(i,:) .* (obsq(k,:)==v);

            essB{k}(i,v) = sum(tmpB);
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now onto the ESS that does require xi
%

% initialization for essA
essA = zeros(Q,Q);


% initialization for essP
essP = zeros(M,Q);
clear tmp1; tmp1 = sum(gamma(:,:,:,1),3);
essP = essP + tmp1';

% initialization for essD
essD = zeros(M,Q);
%         tmp = gamma(i,M,TRUE,:);
%         essD(M,i) = essD(M,i)+sum(tmp);
clear tmp;
tmp = gamma(:,M,TRUE,:);
essD(M,:) = essD(M,:) + (sum(tmp,4))';



% dont do this to conserve memory
% compute the smoothing 2TBN marginals
%[xi] = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi);


for tmpIndex = 1:maxXiLength:T
    lower = tmpIndex;
    
    % somehow MATLAB took too long to compute this min
    %upper = min((tmpIndex + maxXiLength - 1),T);
    upper = tmpIndex + maxXiLength - 1;
    if (upper > T)
        upper = T;
    end
    
    fprintf('Processing time %d to time %d\n',lower,upper);
    
    clear xi;
    xi = fwdback_cxhsmm(PI,A,P,D,H,alpha,beta,phi,lower,upper);


    % expected sufficient statistics for A
    clear tmp; tmp = xi(:,:,:,:,TRUE,:);
    essA = essA + sum(sum(sum(sum(tmp,6),5),4),3);

    
    %-----------------------------------------------------------
    % essP: initial phase
    %-----------------------------------------------------------
    %xi(i,nx_i,n,nx_n,k,t)
    % gamma(i,n,k,t);

    %
    % for i = 1:1:Q,
    % 	for n = 1:1:M,
    % 		clear tmp1;
    % 		tmp1 = gamma(i,n,:,1);
    % 		tmp1 = sum(tmp1);
    %
    %
    % 		clear tmp2;
    % 		tmp2 = xi(:,i,:,n,TRUE,:);
    % 		tmp2 = sum(sum(sum(sum(sum(sum(tmp2))))));
    %
    % 		essP(n,i) = tmp1 + tmp2;
    % 	end;
    % end;

    % moved this out of the loop over tmpIndex
    % vectorize version
    % tmp1(i,n) is the same as tmp1 above
    %clear tmp1; tmp1 = sum(gamma(:,:,:,1),3);


    % tmp2(i,n) is the same as tmp2 above
    clear tmp2; tmp2 = squeeze(sum(sum(sum(xi(:,:,:,:,TRUE,:),6),3),1));
    essP = essP + tmp2';



    %-----------------------------------------------------------
    % essD; essE
    %-----------------------------------------------------------

    for i = 1:1:Q,
        for n = 1:M-1,
            % xi(i,nx_i,n,nx_n,k,t)
            clear tmp;

            % dont know how to vectorize this
            % is there another way to compute the same thing?
            tmp = xi(i,i,n,n+1,FALSE,:);
            essD(n,i) = essD(n,i)+sum(tmp);

            % move this out of the loop over tmpIndex
%             clear tmp;
%             %tmp = xi(i,i,n,n,FALSE,:);
%             tmp = gamma(i,n,FALSE,1:T-1);
%             essE(n,i) = essE(n,i)+sum(tmp);
        end;

    end;

end %for tmpIndex = 1:maxXiLength:T

