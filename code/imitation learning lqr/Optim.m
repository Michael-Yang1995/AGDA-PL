%% Setup


rng(30);
param.iter = 500;
param.Q = 50*eye(n);
param.R = 50*eye(m);
param.K = zeros(m,n);
d = 100; %number of sample of initial points
param.rollout = 1500;
mu = zeros(n,1);
sig = eye(n);
param.X = mvnrnd(mu, sig, d);


param.bd1 = bd1;
param.bd2 = bd2;

%%

param.stepsize = [0.00045,55];
param.iter = 1500;
output0 = PolicyGrad1(A, B, Qe, Re, Q0, R0,lambda, param);

%% compute saddle point value
[Ke, ~, ~] = dlqr(A, B, Qe, Re); % compute optimal control Ke corresponding to Qe, Re
sigma_e = 0;
for t = 1:d
    x = param.X(t,:).';
    for i=1:param.rollout
        sigma_e = sigma_e + x*x.';
        u = -Ke*x;
        x = A*x + B*u;
    end
end
sigma_e = sigma_e/d;

sigma = 0;
for t = 1:d
    x = param.X(t,:).';
    for i=1:param.rollout
        sigma = sigma + x*x.';
        u = -output0.K*x;
        x = A*x + B*u;  
    end
end
sigma = sigma/d;

opt = trace(sigma*output0.Q) + trace(output0.K*sigma*output0.K.'*output0.R) -...
      (trace(sigma_e*output0.Q) + trace(Ke*sigma_e*Ke.'*output0.R)) -...
       lambda*(norm(output0.Q, 'fro')^2+norm(output0.R,'fro')^2);


%% policy gradient
param.iter = 500;
param.stepsize = [0.000220,0.40];
output = PolicyGrad(A, B, Qe, Re, Q0, R0, lambda, param);

%[0.0028,0.55]
%[0.000055,0.40]
%[0.000220,0.40] n = 30, m = 20
%[0.000330,0.15] n = 20, m = 10



%% SVRG
rng(30);
param_SVRG.iter = 250  ;
param_SVRG.rollout = param.rollout;
param_SVRG.Q = param.Q; 
param_SVRG.R = param.R;
param_SVRG.K = param.K;
param_SVRG.X = param.X;
param_SVRG.bd1 = bd1;
param_SVRG.bd2 = bd2;
param_SVRG.stepsize = [0.0000099,0.00315];
param_SVRG.N = 100;
param_SVRG.T = 1;
output3 = SVRG(A, B, Qe, Re,Q0, R0,lambda, param_SVRG);

%[0.00030,0.008]
% [0.0000097,0.0032] n = 30, m = 20
% [0.0000200,0.0032] n = 20, m = 10
%% plot errors
MS = 16;
LW = 2;
LegFont = 12;
iter = 500;


figure
axes('XScale', 'linear', 'YScale', 'log');
hold on
h_PolicyGrad_errors = loglog(0:iter-1, output.errors(1:iter), 'b-', 'LineWidth', LW);
%h_SAGDA_errors = loglog(0:iter-1, SAGDA_result_all.errors(1:iter), 'r-', 'LineWidth', LW);
h_SVRG_errors = loglog(output3.foc(1:250)/d, output3.errors(1:250), 'g-', 'LineWidth', LW);
hold off
set(gca, 'fontsize', 18)
l = legend([ h_PolicyGrad_errors, h_SVRG_errors], 'AGDA','VR-AGDA', 'Location','northeast');
l.FontSize = LegFont;
 %title('(b) Convergence of deterministic GDA', 'fontsize', 18)
xlabel('#grad/n', 'fontsize', 24)
ylabel('$\Vert K_t-K^*\Vert^2+\Vert Q_t-Q^*\Vert^2+\Vert R_t-R^*\Vert^2$','Interpreter','latex', 'fontsize', 15)
