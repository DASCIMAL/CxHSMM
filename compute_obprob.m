function [H lnHScale] = compute_obprob(B,obsq,scaleStr);


% function [H lnHScale] = compute_obprob(B,obsq,scaleStr);
% To generalize and speed up the computation, we compute ONCE the observation
% probability at every time slice. It returns:
% 		H is an Q x T matrix where 
%     H(i,t) = Pr(y_t|x_t=i)=\prod_k=1^K Pr(y_t^k | x_t =i)

% Last updated:  Hung Bui 30/09/05
%                DQ Phung 26/09/05

if (nargin == 2)
    scaleFlag = 0;
elseif (nargin == 3)&&(strcmp(scaleStr,'scale'))
    scaleFlag = 1;
else
    fprintf('Something wrong in the scale argument in compute_obprob\n');
    exit;
end

K = length(B);    % observation vector cardinality 
Q = size(B{K},1); % state space
T = size(obsq,2); % observation length
NULL = -1;        % the special 'null symbol used in the "null trick"

% allocate memory
H = ones(Q,T);
lnHScale = 0;

for t=1:T
	for k=1:K
        % here is the NULL trick
        % ignore all values equal to NULL
        if ~(obsq(k,t) == NULL)
            H(:,t) = H(:,t) .* B{k}(:,obsq(k,t));
            
            if (scaleFlag)
                % rescale H immediately, make all columns of H sum to 1
                tmp = 1/sum(H(:,t));
                H(:,t) = H(:,t) .* tmp; 
                lnHScale = lnHScale + log(tmp);
            end
        end
	end
end

% that's it!
