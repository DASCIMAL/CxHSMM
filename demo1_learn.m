 % testing EM algorithm
    % Q : number of hidden states
    % M : number of Coxian phases
    % K : number of "features" in each observation
    % V : sizes for all features
    % mobsq : the set of obs sequences
    % 20 : max number of EM iterations
    % 1e-30 : loglik difference threshold
    % 'observe_state' : observation including direct state value
    % (in that case the last row in each of mobsq{i} is the state sequence
    %

[PI,A,P,D,B,loglik] = em_cxhsmm(Q,M,K,V,mobsq,5,1e-30,'observe_state');

