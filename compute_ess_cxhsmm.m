% function [essPI,essA,essP,essD,essE,essB,loglik] = compute_ess_cxhsmm(PI,A,P,D,B,V,mobsq)
% This function compute the ESS FOR ALL OBSERVATION SEQUENCES
%
% There are N iid obsq contained in variable mobsq. It basically iterate through each obsq
% and add up the ESS for each parameters.
function [essPI,essA,essP,essD,essE,essB,loglik] = compute_ess_cxhsmm(PI,A,P,D,B,V,mobsq)

% Last updated: 
%       Hung Bui 14/10/2005
%       DQ PHUNG 26/09/2005	
%       VT Duong 12/05/2005	

% housekeeping thing
clear essPI; clear essA; clear essP; clear essD; clear essE; clear essB; clear loglik;

% number of iid observation sequences
N = length(mobsq); 
% observation dim
K = length(V);

% control maximum number of time slices allocated
% in the 2TBN matrix "xi"
% use the largest number that your RAM can handle
%
% if you don't care about memory usage, just comment out this line
% note that the memory needed for each 2TBN segment will be
% M^2 * Q^2 * 2 * maxLength2TBN * 4 bytes
% where M is no phases, Q is no of states
%
% if you want to know roughly the largest memory block your machine
% can allocate, then try
% x=feature('DumpMem');
% 
maxLength2TBN = 2;

% compute the ESS for the first sequence separately
if exist('maxLength2TBN')
    [essPI,essA,essP,essD,essE,essB,loglik] = ess_cxhsmm(PI,A,P,D,B,V,mobsq{1},maxLength2TBN);
else
    [essPI,essA,essP,essD,essE,essB,loglik] = ess_cxhsmm(PI,A,P,D,B,V,mobsq{1});
end


% now loop through the rest
for n = 2:N

	clear tmpPI; clear tmpA; clear tmpP; clear tmpD; clear tmpE; clear tmpB; clear tmplog;

	% collect ESS from obs sequence nth
    if exist('maxLength2TBN')
        [tmpPI,tmpA,tmpP,tmpD,tmpE,tmpB,tmplog] = ess_cxhsmm(PI,A,P,D,B,V,mobsq{n},maxLength2TBN);
    else
        [tmpPI,tmpA,tmpP,tmpD,tmpE,tmpB,tmplog] = ess_cxhsmm(PI,A,P,D,B,V,mobsq{n});
    end
	
	essPI = essPI + tmpPI; % add up for PI
	essA = essA + tmpA;    % add up for A
	essP = essP + tmpP;    % add up for P
	essD = essD + tmpD;    % add up for D
	essE = essE + tmpE;    % add up for E

	% add up for B
	for (k=1:K)
		essB{k} = essB{k} + tmpB{k};
	end

	% update the log-likelihood
	loglik = loglik + tmplog;
end

% that's it!
