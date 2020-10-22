


stepsize_x = 0.000031:0.000001:0.000050;
stepsize_y = 0.0091:0.0001:0.0110;
%final_error_SGDA = zeros(100,100);
final_error_AGDA = zeros(100,100);
final_error_SAGDA = zeros(100,100);


param1.iter = 100;
param1.tol = 10^(-10);
param1.x_0 = random('normal',0,1,m,1);
param1.y_0 = random('normal',0,1,n,1);
param1.stepsize = [0.00025,0.048];

%opt = Optimizer1(A, F, lambda, y0, param1, 'AGDA');
%g_opt = obj_func(A, F, lambda, y0, opt.x, opt.y);

%%


for k = 1:20  %x stepsize
for m = 1:20  %y stepsize
    param1.stepsize = [stepsize_x(k),stepsize_y(m)];
    %SGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SGDA', opt.x, opt.y, g_opt, 1);
    AGDA_result_all1 = Optimizer(A, F, lambda, y0, param1, 'AGDA', opt.x, opt.y, g_opt, 1);
    %SAGDA_result_all = Optimizer(A, F, lambda, y0, param, 'SAGDA', opt.x, opt.y, g_opt, 1);
    %SGDA_g = max_oracle(A, F, lambda, y0, SGDA_result_all.x);
    %AGDA_g = max_oracle(A, F, lambda, y0, AGDA_result_all.x);
    %final_error_SGDA(k,m) = norm(SGDA_g - g_opt);
    %final_error_AGDA(k,m) = norm(AGDA_g - g_opt);  
    %final_error_SGDA(k,m) = SGDA_result_all.errors(param.iter);
    final_error_AGDA(k,m) = AGDA_result_all1.errors(param1.iter);
    %final_error_SAGDA(k,m) =  SAGDA_result_all.errors(param.iter);
end
end

%min_SGDA = min(final_error_SGDA(:));
min_AGDA = min(final_error_AGDA(:));
%[row1,col1] = find(final_error_SGDA==min_SGDA);
[row2,col2] = find(final_error_AGDA==min_AGDA);
        
