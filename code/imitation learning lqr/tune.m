


stepsize_x = 0.0001:0.0001:0.001;
stepsize_y = 10:10:100;
final_error_AGDA = zeros(100,100);
final_error_SAGDA = zeros(100,100);

% 0.0009, 100


param1.iter = 30;
param1.Q = 50*eye(n);
param1.R = 50*eye(m);
param1.K = zeros(m,n);
d = 500; %number of sample of initial points
param1.rollout = 1500;
param1.X = param.X;


param1.bd1 = bd1;
param1.bd2 = bd2;

%opt = Optimizer1(A, F, lambda, y0, param1, 'AGDA');
%g_opt = obj_func(A, F, lambda, y0, opt.x, opt.y);

%%


for k = 1:10  %x stepsize
for m = 1:10  %y stepsize
    param1.stepsize = [stepsize_x(k),stepsize_y(m)];
    AGDA_result_all1 =  PolicyGrad(A, B, Qe, Re, Q0, R0, lambda, param1);
    %SAGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 1);
    %final_error_AGDA(k,m) = norm(AGDA_g - g_opt);  
    %final_error_SGDA(k,m) = SGDA_result_all.errors(param.iter);
    final_error_AGDA(k,m) = AGDA_result_all1.errors(end);
    %final_error_SAGDA(k,m) =  SAGDA_result_all.errors(param.iter);
end
end

%min_SGDA = min(final_error_SGDA(:));
min_AGDA = min(final_error_AGDA(:));
%[row1,col1] = find(final_error_SGDA==min_SGDA);
[row2,col2] = find(final_error_AGDA==min_AGDA);
        
