% This function return, at each time t,
% a vector of state sequence, sorted in descending order
% according to its marginal probability
% e.g.
% rankedState(1,t) is the best state at time t
% rankedState(1:i,t) is the set of i best states at time t
%
% the second output argument simply returns the
% sorted probabilities
% 
% Last updated:  Hung Bui 05/10/2005

function [rankedState probs] = smstate_decode(PI,A,P,D,H)

% forward pass 
[alpha,phi] = forward_cxhsmm(PI,A,P,D,H);

% backward pass
[beta] = backward_cxhsmm(PI,A,P,D,H,phi);

% compute the one-slice smoothing marginal
% gamma(i,m,e,t) = Pr(x_t=i,m_t=m,e_t=e | Y)
gamma = alpha .* beta;

stateMarg = squeeze(sum(sum(gamma,2),3));

[probs rankedState] = sort(stateMarg,1,'descend');
