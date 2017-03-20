function[PI,A,P,D,B] = init_cxhsmm(Q,M,K,V,initStr) 

% function[PI,A,P,D,B] = init_cxhsmm(Q,M,K,V) 
% Input
%    Q: number of states 
%    M: number of phases
%    K: number of observables, i.e., dimension of the observation vector,
%        we can get K from parameter V, but I include here for clarity
%    V: a column 1 x K vector, where V(k) is the number of alphabets for
%        kth observation elements
% It returns:
%    PI: the initial probability of states 
%    A: a zero-diag-entry transition probability of states (A_ii = 0)
%    P: the intial probability of entering a Coxian phase 
%        P(m,i) : prob. of entering phase m given state i  
%    D: the terminating probability of phases. 
%        D(m,i): prob. of phase m 'go to ends' given current state i
%    B: is a cell array of size K
%        B{k} is an Vk x Q  observation matrix
%        B{k}(v,i) = Pr(y_t^k = v | x_t = i)
%
% Last updated: DQ Phung   26/09/2005
%               VT Duong   24/05/2005  

% initialization method, should be passed from calling function

if (nargin ~= 5)
    printf('Something wrong with number of arguments in init_cxhsmm.\n');
    exit;
end

initmethod = initStr;  % CURRENTLY USE 'random'

% initialization for PI 
if (strcmp(initmethod,'uniform'))
	PI = ones(1,Q);
elseif (strcmp(initmethod,'random'))
	PI = rand(1,Q);
else
    printf('I dont know this init method in init_cxhsmm.\n');
    exit;
end
PI = normalize(PI);

% initialization for A
if (strcmp(initmethod,'uniform'))
	A = ones(Q,Q);
elseif (strcmp(initmethod,'random'))
	A = rand(Q,Q);
end

% ensure that the diagonal entries are zero, since this is
% a semi-Markov chain.
for i = 1:1:Q,
	A(i,i) = 0;
end;
A = mk_stochastic(A);


% entering phase initial probability 
if (strcmp(initmethod,'uniform'))
	P = ones(M,Q);
elseif (strcmp(initmethod,'random'))
	P = rand(M,Q);
end
P = normalize(P,1);

% ending phase probability, note that we don't need to normalize!
if (strcmp(initmethod,'uniform'))
		  D = 0.5 * ones(M,Q);  % 1/2 chance going to end
elseif (strcmp(initmethod,'random'))
	D = rand(M,Q);
end

% observation probability matrices
for k=1:K
	if (strcmp(initmethod,'uniform'))
		B{k} = ones(Q,V(k));
	elseif (strcmp(initmethod,'random'))
		B{k} = rand(Q,V(k));			  	
	end
	B{k} = mk_stochastic(B{k});
end

% that's all for initialization
