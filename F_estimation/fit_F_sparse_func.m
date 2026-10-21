function [F, Cost_opt, F_fminunc, Cost_opt_fminunc] = fit_F_sparse_func(traj_dat, gen_dat, n)
%Estimates matrix of fitness effects, F, from trajectory and generation
%data. It can be applied to trajectories where data is not available for
%every generation (i.e. for sparse trajectories)


% Minimisation procedure
% Use fitness effects of zero as starting values
L=n*(n+1)/2;
v0=zeros(L-1,1);

% Set optimization options
fminunc_options = optimoptions('fminunc', ...
        'MaxFunctionEvaluations', 100000, ...
        'TolFun', 1e-9, ...  % Function tolerance
        'TolX', 1e-9, ...    % Step tolerance
        'MaxIter', 100000, ... % Maximum number of iterations
        'Display', 'off');  % Do not display iteration information

% Minimize CostSparseFn(v, traj_dat, gen_dat)
f=@(v)CostSparseFn(v, traj_dat, gen_dat);
[v_opt, Cost_opt]=fminunc(f, v0, fminunc_options);

F=SMatVec([0;v_opt]);

F_fminunc = F;
Cost_opt_fminunc = Cost_opt;


% test for non-sensical fintess values
has_value_smaller_neg_one = sum(v_opt < -1.0);

if has_value_smaller_neg_one
    ps_options = optimoptions('patternsearch', ...
                'MaxFunctionEvaluations', 100000, ...
                'TolFun', 1e-9, ...  % Function tolerance
                'TolX', 1e-9, ...    % Step tolerance
                'MaxIter', 100000, ... % Maximum number of iterations
                'Display', 'off');  % Display iteration information
    % Contraints
    lb =-ones(L-1,1);           % Lower bound of -1  on all elements of v
    ub=Inf(L-1,1);
    % Minimize CostSparseFn(v, traj_dat, gen_dat)
    f=@(v)CostSparseFn(v, traj_dat, gen_dat);
    [v_opt, Cost_opt]=patternsearch(f, v0, [], [], [], [], lb, ub, ps_options);

    F=SMatVec([0;v_opt]);

end


end