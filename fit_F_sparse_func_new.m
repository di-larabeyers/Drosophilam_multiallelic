%%%%%%%%% updated fit_F_sparse_func.m


function [F, Cost_opt] = fit_F_sparse_func_new(traj_dat, gen_dat, n)% Estimates matrix of fitness effects, F, from trajectory and generation
% data. It can be applied to sparse trajectories where data is not
% available for every generation.
%
% Uses:
%   - MultiStart
%   - fmincon (SQP)
%   - Shahshahani-weighted cost function with ridge regularization
%     (CostSparseShahRidge)

%% Parameters

L = n*(n+1)/2;

v_init = zeros(L-1,1);

lambda = 6.1e-3;
eps = 1e-6;

n_starts = 80;
s_start = 0.1;

%% fmincon options

fmincon_options = optimoptions('fmincon', ...
    'Algorithm','sqp', ...
    'MaxFunctionEvaluations',10000, ...
    'TolFun',1e-7, ...
    'TolX',1e-7, ...
    'MaxIter',1000, ...
    'Display','off');

lb = -0.999 * ones(L-1,1);
ub = Inf(L-1,1);

solver = 'fmincon';
solver_options = fmincon_options;

%% MultiStart object

ms = MultiStart( ...
    'UseParallel', false, ...
    'Display', 'off');

%% Cost function

f = @(v) CostSparseShahRidge(v, traj_dat, gen_dat, eps, lambda);

%% Generate start points

start_points = [v_init'; ...
                s_start * randn(n_starts, L-1)];

% Enforce lower bound
start_points = max(start_points, -0.999);

sp = CustomStartPointSet(start_points);

%% Optimization problem

problem = createOptimProblem( ...
    solver, ...
    'objective', f, ...
    'x0', v_init, ...
    'lb', lb, ...
    'ub', ub, ...
    'options', solver_options);

%% Run optimization

[v_opt, Cost_opt] = run(ms, problem, sp);

%% Convert vector to fitness matrix

F = SMatVec([0; v_opt]);


end