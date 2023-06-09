function fullStairs(tspan,u,varargin)
% fullStairs draws a stairstep plot of u vs. tspan, much like stairs.m. 
% The difference is that the number of stairs (i.e., horizontal lines) 
% equals the length of the input vectors.
%
% Inputs:
%   tspan, a vector of length K+1.
%   u, a K x M matrix.
%
% Optional arguments (line/marker specs, etc.) are allowed, provided they
% follow the syntax accepted by stairs.m.
%
% This function is designd to plot the input signals to a discrete-time 
% dynamical system in M Monte Carlo runs. fullStairs assumes that the input
% signal u is defined such that
%
%   u(t,:) = u(1,:) for tspan(1) <= t <= tspan(2)
%   u(t,:) = u(2,:) for tspan(2) < t <= tspan(3)
%        .
%        .
%        .
%   u(t,:) = u(K,:) for tspan(K) < t <= tspan(K+1).
%
% fullStairs modifies the behavior of stairs.m so that a final stair is 
% drawn at height u(K,m) corresponding to t in (tspan(K),tspan(K+1)], for
% each m = 1, ..., M.

% Check dimensions.
if ~isvector(tspan)
    error('tspan must be a vector.')
end
nt = length(tspan);
nu = size(u,1);
if nt ~= nu + 1
    error('length(tspan) must equal size(u,1) + 1.')
end

% Input augmentation.
u = [u;u(end,:)];

% Stairstep plot.
stairs(tspan,u,varargin{:})

end

