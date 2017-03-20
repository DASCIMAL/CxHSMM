%classify
obsq = mobsq{1}; % pick one sequence

% must do this first
[H lnH] = compute_obprob(B,obsq,'scale');

% there are two methods for infering the hidden states
% below is a simple argmax P(x_t | obs) for each t

[rankedState probs] = smstate_decode(PI,A,P,D,H);
smoothedLabels=rankedState(1,:);

% 
% Here is the viterbi decode which returns
% the most likely sequence, rather than a sequence
% of most likely states as above
%
% argmax P(x_{1:T} | obs)
%
% If you're not sure which method to use,
% then try the viterbi first

[lvtbsq lprob] = viterbi_cxhsmm(PI,A,P,D,H,'uselog');
viterbiLabels = lvtbsq(1,:);
fprintf('Look at smoothedLabels and viterbiLabels for results.\n');
