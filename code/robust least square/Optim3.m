clear
clc

%% dimension of data
n = 1000;
m = 500;

%% generate A, x* and y0
rng(400);
mu = zeros(m,1);
sigma = zeros(m,m);
for i=1:m
    for j=1:m
        sigma(i,j) = 2^(-abs(i-j)/10);
    end
end
A = mvnrnd(mu, sigma, n);
x_star = random('normal', 0, 1, m,1);
y0 = A*x_star + random('normal', 0, 0.01, n,1);

p = floor(n*0.95);  % rank = n*0.9
U = randn(n, n);  
v = rand(p,1)*1.6+0.2;   % simga has eigenvalues between (0.2,1.8)
v = sqrt(v);
sigma = diag([v; zeros(n-p,1)]);
U = orth(U);
F = transpose(U)*sigma*U; 


lambda = 1.5;

%% Setup

param.iter = 200;
param.tol = 8*10^(-10);
param.x_0 = random('normal',0,1,m,1);
param.y_0 = random('normal',0,1,n,1);
param.stepsize = [0.000031,0.0051];
 

opt = Optimizer1(A, F, lambda, y0, param, 'AGDA');
g_opt = obj_func(F*A, F, lambda, F*y0, opt.x, opt.y);

%% Stepsize of SGDA
%param.iter = 200;
%param.stepsize = [0.0005,0.00008];
%SGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SGDA', opt.x, opt.y, g_opt, 1);

%% Stepsize of AGDA
param.iter = 5000;
param.stepsize = [0.0000328,0.0178];  %[0.0000328,0.0178];
AGDA_result_all = Optimizer(A, F, lambda, y0, param, 'AGDA', opt.x, opt.y, g_opt, 1);%
%AGDA_result_all_2 = Optimizer(A, F, lambda, y0, param, 'AGDA', opt.x, opt.y, g_opt, 2);

%%  Stepsize of EG two time scale
param.iter = 2500;
param.stepsize = [0.000015, 0.5500];   %0.00017
EG2_result_all = EG(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 1);
%EG2_result_all_2 = Optimizer(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 2);

%%  Stepsize of EG 
param.iter = 2500;
param.stepsize =  [0.0000328,0.0178];   %0.000015
EG_result_all = EG(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 1);
%EG_result_all_2 = Optimizer(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 2);



%% Stepsize of Stoc-AGDA and backtracking parameters
param.iter = 5000;
param.iter = param.iter*n;
param.stepsize = [0.00003,0.01];   %[0.00003,0.01];
SAGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 1);
%SAGDA_result_all_2 = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 2);





%% Stepsize and iters of SVRG 

% find the opt point for SVRG
param_svrg.iter = 10000;
param_svrg.x_0 = param.x_0;
param_svrg.y_0 = param.y_0;
param_svrg.stepsize =[0.0000008000, 0.0003100];
param_svrg.N = 1000;
param_svrg.T = 1;
SVRG_result_all_opt = SVRG(A, F, lambda, y0, param_svrg, opt.x, opt.y, g_opt, 1);



param_svrg.iter = 5100;
param_svrg.x_0 = param.x_0;
param_svrg.y_0 = param.y_0;
param_svrg.stepsize = [0.0000008000, 0.0003100];
param_svrg.N = 1000;
param_svrg.T = 1;
SVRG_result_all = SVRG(A, F, lambda, y0, param_svrg, SVRG_result_all_opt.x,  SVRG_result_all_opt.y, g_opt, 1);
%SVRG_result_all_2 = SVRG(A, F, lambda, y0, param_svrg, opt.x, opt.y, g_opt, 2);


%N=800, [0.0000000728, 0.0003680]
%N =1000,  [0.0000000547, 0.0002960]  627 at 750000foc
%N = 500, [0.0000001360, 0.000600]  931
%N = 350 [0.0000001684, 0.000843]  734 at 750600foc








%% plot errors
MS = 16;
LW = 4;
LegFont = 16;
iter = 5000;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_AGDA_errors = loglog(0:iter-1, AGDA_result_all.errors(1:iter), 'b-', 'LineWidth', LW);
h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all.errors(1:iter), 'r-', 'LineWidth', LW);
h_EG_errors = loglog(2*(0:(2500-1)), EG2_result_all.errors(1:2500), 'c-', 'LineWidth', LW);
h_SVRG_errors = loglog(SVRG_result_all.foc(1:2500)/n, SVRG_result_all.errors(1:2500), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_AGDA_errors,h_SAGDA_errors, h_EG_errors h_SVRG_errors], 'AGDA', 'Stoc-AGDA', 'EG','AGDA-SVRG', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 30)
ylabel('$\Vert x_t-x^*\Vert^2+\Vert y_t-y^*\Vert^2$','Interpreter','latex', 'fontsize', 30)


%% plot
MS = 16;
LW = 2;
LegFont = 12;
iter = 5000;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_AGDA_errors = loglog(0:iter-1, AGDA_result_all_2.errors(1:iter), 'b-', 'LineWidth', LW);
h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all_2.errors(1:iter), 'r-', 'LineWidth', LW);
h_SVRG_errors = loglog(SVRG_result_all_2.foc(1:2500)/n, SVRG_result_all_2.errors(1:2500), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_AGDA_errors,h_SAGDA_errors, h_SVRG_errors], 'AGDA', 'Stoc-AGDA','AGDA-SVRG', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 24)
ylabel('$P_t$','Interpreter', 'latex', 'fontsize', 24)