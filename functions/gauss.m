function X = gauss(a,b,m,n)
% TRIRND samples from the symmetric triangle distribution.
% Inputs:
%   a, the minimum of the support of the distribution.
%   b, the maximum of the support of the distribution.
%   m, the number of rows in the sample matrix.
%   n, the number of columns in the sample matrix.
%
% Output:
%   X, an m x n matrix of independent samples.
rng(1)
% define distribution
pd = makedist('Normal','mu',(a+b)/2,'sigma',(b-(b+a)/2)/2.58);
% generate samples
X = random(pd,m,n);

end
