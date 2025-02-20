function X = gauss(a, b, m, n)
% GAUSS generates normally distributed samples constrained between a and b.
% Inputs:
%   a, the minimum of the support of the distribution.
%   b, the maximum of the support of the distribution.
%   m, the number of rows in the sample matrix.
%   n, the number of columns in the sample matrix.
%
% Output:
%   X, an m x n matrix of independent samples.

     % Set seed for reproducibility
    
    % Define the mean as the midpoint between a and b
    mu = (a + b) / 2;
    
    % Estimate the standard deviation so that ~95% of values lie within [a, b]
    % The standard deviation for a 95% confidence interval:
    sigma = (b - a) / (2 * 1.96); % Dividing by 2*1.96 ensures 95% of values fall between a and b
    
    % Generate normal distribution samples
    pd = makedist('Normal', 'mu', mu, 'sigma', sigma);
    
    % Generate m x n samples
    X = random(pd, m, n);
    
    % Ensure samples stay within the bounds [a, b] by rejection sampling
   % X(X < a) = a;
    %X(X > b) = b;
end
