clear
clc

%% dimension of data
n = 546;
m = 8;

%% generate A, x* and y0
A = readtable('qsar_aquatic_toxicity.csv');
A = A{:,:};
y0 = A(:,end);
A = A(:,1:end-1);
F = eye(n,n);

mean_A = mean(A);
var_A = var(A);  
A = (A-mean(A))./sqrt(var_A); 

lambda = 2;

%% Setup

param.iter = 200;
param.tol = 10^(-12);
param.x_0 = random('normal',0,1,m,1);
param.y_0 = random('normal',0,1,n,1);
param.stepsize = [0.0004,0.048];


opt = Optimizer1(A, F, lambda, y0, param, 'AGDA');
g_opt = obj_func(A, F, lambda, y0, opt.x, opt.y);

%% Stepsize of SGDA
%param.iter = 200;
%param.stepsize = [0.0005,0.00008];
%SGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SGDA', opt.x, opt.y, g_opt, 1);

%% Stepsize of AGDA
param.iter = 200;
param.stepsize = [0.000448,0.060];
AGDA_result_all = Optimizer(A, F, lambda, y0, param, 'AGDA', opt.x, opt.y, g_opt, 1);
%AGDA_result_all_2 = Optimizer(A, F, lambda, y0, param, 'AGDA', opt.x, opt.y, g_opt, 2);

%[0.000448,0.059]


%%  Stepsize of EG two time scale
param.iter = 100;
param.stepsize = [0.00026,0.05];   %0.00025
EG2_result_all = EG(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 1);
%EG2_result_all_2 = Optimizer(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 2);


%%  Stepsize of EG
param.iter = 100;
param.stepsize = [0.00025,0.00025];   %0.00025
EG_result_all = EG(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 1);
%EG_result_all_2 = Optimizer(A, F, lambda, y0, param, opt.x, opt.y, g_opt, 2);



%% Stepsize of Stoc-AGDA and backtracking parameters
param.iter = 200;
param.iter = param.iter*n;
param.stepsize = [0.0002,0.038];
SAGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 1);
%SAGDA_result_all_2 = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 2);


%%
rng(400);
param_svrg.iter = 300;
param_svrg.x_0 = param.x_0;
param_svrg.y_0 = param.y_0;
param_svrg.stepsize = [0.00000620, 0.00082];
param_svrg.N = 350;
param_svrg.T = 1;
SVRG_result_all_opt = SVRG(A, F, lambda, y0, param_svrg, opt.x, opt.y, g_opt, 1);
g_opt = obj_func(A, F, lambda, y0, SVRG_result_all_opt.x, SVRG_result_all_opt.y);



%% Stepsize and iters of SVRG





rng(400);
param_svrg.iter = 180;
param_svrg.x_0 = param.x_0;
param_svrg.y_0 = param.y_0;
param_svrg.stepsize = [0.00000450, 0.00052];
param_svrg.N = 350;
param_svrg.T = 1;
SVRG_result_all = SVRG(A, F, lambda, y0, param_svrg, SVRG_result_all_opt.x, SVRG_result_all_opt.y, g_opt, 1);
%SVRG_result_all_2 = SVRG(A, F, lambda, y0, param_svrg, SVRG_result_all_opt.x, SVRG_result_all_opt.y, g_opt, 2);
% N = 350, [0.00000610, 0.00082] 

%% plot errors
MS = 16;
LW = 4;
LegFont = 16;
iter = 100;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_AGDA_errors = loglog(0:iter-1, AGDA_result_all.errors(1:iter), 'b-', 'LineWidth', LW);
h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all.errors(1:iter), 'r-', 'LineWidth', LW);
h_EG_errors = loglog(2*(0:(50-1)), EG2_result_all.errors(1:50), 'c-', 'LineWidth', LW);
h_SVRG_errors = loglog(SVRG_result_all.foc(1:61)/546, SVRG_result_all.errors(1:61), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_AGDA_errors,h_SAGDA_errors, h_EG_errors, h_SVRG_errors], 'AGDA', 'Stoc-AGDA','EG', 'AGDA-SVRG', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 30)
ylabel('$\Vert x_t-x^*\Vert^2+\Vert y_t-y^*\Vert^2$','Interpreter','latex', 'fontsize', 30)


%% plot
MS = 16;
LW = 2;
LegFont = 12;
iter = 100;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_AGDA_errors = loglog(0:iter-1, AGDA_result_all_2.errors(1:iter), 'b-', 'LineWidth', LW);
h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all_2.errors(1:iter), 'r-', 'LineWidth', LW);
h_SVRG_errors = loglog(SVRG_result_all_2.foc(1:61)/546, SVRG_result_all_2.errors(1:61), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_AGDA_errors,h_SAGDA_errors, h_SVRG_errors], 'AGDA', 'Stoc-AGDA','AGDA-SVRG', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 24)
ylabel('$P_t$','Interpreter', 'latex', 'fontsize', 24)


    
    
    
    