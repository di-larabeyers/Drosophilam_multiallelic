function cost = CostSparseShahRidge(w, X, gens, eps, lambda)
% Cost function for sparse trajectories with Shahshahani weighting
% and ridge regularization.
%
% X    : allele frequency trajectory (n x (T+1))
% gens : generation numbers corresponding to columns of X
% w    : parameter vector of length n(n+1)/2 - 1
% eps  : stabilization constant
% lambda : ridge penalty strength

if nargin < 4
    eps = 1e-8;     % default epsilon
end

if nargin < 5
    lambda = 0;     % default lambda
end

if size(X,2) ~= length(gens)
    error("Length of gens must match number of columns in X");
end

T = size(X,2)-1;

v = [0; w];
F = SMatVec(v);

S = 0;

for t = 1:T

    x_obs = X(:,t);
    xp_obs = X(:,t+1);

    delta_T = gens(t+1) - gens(t);

    y = x_obs;

    % propagate model through missing generations
    for s = 1:delta_T
        V = diag(y) - y*y';
        D = V*F*y/(1 + y'*F*y);
        y = y + D;
    end

    % predicted state at next observed timepoint
    Y = y;

    current_diff = xp_obs - Y;

    S = S + sum(current_diff.^2 ./ max(x_obs, eps));

end

cost = S + lambda*(w'*w);

end