function X = trirnd(a,c,m,n)
% TRIRND samples from the symmetric triangle distribution.
% Inputs:
%   a, the minimum of the support of the distribution.
%   c, the maximum of the support of the distribution.
%   m, the number of rows in the sample matrix.
%   n, the number of columns in the sample matrix.
%
% Output:
%   X, an m x n matrix of independent samples.
rng(1);
% define distribution
pd = makedist('Triangular','a',a,'b',(a+c)/2,'c',c);

% generate samples
X = random(pd,m,n);

end

