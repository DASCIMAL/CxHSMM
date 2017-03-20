function [M, c] = normalize(A, dim)


% ACKNOWLEDGEMENT: THIS NORMALIZE FUNCTION WAS DOWNLOADED FROM THE WEB: 
% www.mit.edu/afs/sipb.mit.edu/user/arolfe/FullBNT/KPMtools


% NORMALISE Make the entries of a (multidimensional) array sum to 1
% [M, c] = normalise(A, dim)
% c is the normalizing constant

% If dim is specified, we normalise the specified dimension only,
% otherwise we normalise the whole array.


if nargin < 2
  c = sum(A(:));
  % Set any zeros to one before dividing
  % This is valid, since c=0 => all i. A(i)=0 => the answer should be 0/1=0
  d = c + (c==0);
  M = A / d;
else % Keith Battocchi
  s=sum(A,dim);
  L=size(A,dim);
  d=length(size(A));
  v=ones(d,1);
  v(dim)=L;
  %c=repmat(s,v);
  c=repmat(s,v');
  M=A./c;
end
