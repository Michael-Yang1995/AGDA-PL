clear
clc

%% dimension of data
n = 1000;
m = 500;

%% generate A, x* and y0
rng(400);
mu = zeros(m,1);
sigma = eye(m);
A = mvnrnd(mu, sigma, n);
x_star = random('normal', 0, 1, m,1);
y0 = A*x_star + random('normal', 0, 0.01, n,1);
F = eye(n,n);

lambda = 3;

%% Setup

param.iter = 200;
param.tol = 5*10^(-11);
param.x_0 = random('normal',0,1,m,1);
param.y_0 = random('normal',0,1,n,1);
param.stepsize = [0.00031,0.051];


opt = Optimizer1(A, F, lambda, y0, param, 'AGDA');
g_opt = obj_func(A, F, lambda, y0, opt.x, opt.y);

%% Stepsize of SGDA
%param.iter = 200;
%param.stepsize = [0.0005,0.00008];
%SGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SGDA', opt.x, opt.y, g_opt, 1);

%% Stepsize of AGDA
param.iter = 150;
param.stepsize = [0.00031,0.060];
AGDA_result_all = Optimizer(A, F, lambda, y0, param, 'AGDA', opt.x, opt.y, g_opt, 1);%%
%AGDA_result_all_2 = Optimizer(A, F, lambda, y0, param, 'AGDA', opt.x, opt.y, g_opt, 2);

%[0.00031,0.0485]
%[0.00031,0.060]


%%  Stepsize of EG two time scale
param.iter = 100;
param.stepsize = [0.00020, 0.069];   %0.00017
EG2_result_all = EG(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 1);
%EG2_result_all_2 = Optimizer(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 2);

%%  Stepsize of EG 
param.iter = 100;
param.stepsize = [0.00017, 0.00017];   %0.00017
EG_result_all = EG(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 1);
%EG_result_all_2 = Optimizer(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 2);


%% Stepsize of Stoc-AGDA and backtracking parameters
param.iter = 200;
param.iter = param.iter*n;
param.stepsize = [0.00029,0.04];
SAGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 1);
%SAGDA_result_all_2 = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 2);

%[0.00029,0.04]

%% Stepsize and iters of SVRG 
rng(400);
param_svrg.iter = 150;
param_svrg.x_0 = param.x_0;
param_svrg.y_0 = param.y_0;
param_svrg.stepsize = [0.000001110, 0.000195];
param_svrg.N = 500;
param_svrg.T = 1;
SVRG_result_all = SVRG(A, F, lambda, y0, param_svrg, opt.x, opt.y, g_opt, 1);
%SVRG_result_all_2 = SVRG(A, F, lambda, y0, param_svrg, opt.x, opt.y, g_opt, 2);

% N =200, [0.00000189, 0.000390]  0.03 50400foc
% N = 1000, [0.000000591, 0.000300]  0.15 50000foc
% N = 350, [0.00000149, 0.000190] 0.0169  51300foc
% N = 500, [0.000001110, 0.000195] 0.0088 potential 50000foc
%% plot errors
MS = 16;
LW = 4;
LegFont = 16;
iter = 150;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_AGDA_errors = loglog(0:iter-1, AGDA_result_all.errors(1:iter), 'b-', 'LineWidth', LW);
h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all.errors(1:iter), 'r-', 'LineWidth', LW);
h_EG_errors = loglog(2*(0:(75-1)), EG2_result_all.errors(1:75), 'c-', 'LineWidth', LW);
h_SVRG_errors = loglog(SVRG_result_all.foc(1:101)/n, SVRG_result_all.errors(1:101), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_AGDA_errors,h_SAGDA_errors,h_EG_errors, h_SVRG_errors], 'AGDA', 'Stoc-AGDA', 'EG', 'AGDA-SVRG', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 30)
ylabel('$\Vert x_t-x^*\Vert^2+\Vert y_t-y^*\Vert^2$','Interpreter','latex', 'fontsize', 30)


%% plot
MS = 16;
LW = 2;
LegFont = 12;
iter = 150;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_AGDA_errors = loglog(0:iter-1, AGDA_result_all_2.errors(1:iter), 'b-', 'LineWidth', LW);
h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all_2.errors(1:iter), 'r-', 'LineWidth', LW);
h_SVRG_errors = loglog(SVRG_result_all_2.foc(1:101)/n, SVRG_result_all_2.errors(1:101), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_AGDA_errors,h_SAGDA_errors, h_SVRG_errors], 'AGDA', 'Stoc-AGDA','AGDA-SVRG', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 24)
ylabel('$P_t$','Interpreter', 'latex', 'fontsize', 24)