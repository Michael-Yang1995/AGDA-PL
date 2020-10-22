function Result = SVRG(A, F, lambda, y0, param, x_opt, y_opt, g_opt, error_type)



 %% Parameters
    x_0 = param.x_0;
    y_0 = param.y_0;
    iter = param.iter;
    N = param.N;    
    T = param.T;
    tau = param.stepsize;
    [n,~] = size(y_0);
    

    %% Initialization
    x = x_0;
    y = y_0;  
    %errors = zeros(1,iter*N*T+1);
    errors = zeros(1,iter*T+1);
    FA = F*A;
    Fy0 = F*y0;
    %foc = zeros(1, iter*N*T+1);
    foc = zeros(1,iter*T+1);
    [lx,ly] = size(FA);
    pos = 1;
   
    
    %% find error for x_0;
    if error_type == 2 
        max = max_oracle(FA, F, lambda, Fy0, x, y);
        errors(1) = max - g_opt + max - obj_func(FA, F, lambda, Fy0, x, y);
    elseif error_type == 1
        errors(1) = (norm(x - x_opt))^2 + (norm(y - y_opt))^2;
    end
     
    %% iterations for svrg
    for k = 1:iter
        for t = 1:T
            % store gradients of snapshots
            snap_grad = grad_func(FA, F, lambda, Fy0, x, y, 0, 'det');
            snap_x = x;
            snap_y = y;
            for j = 1:N
                % compute new x and y using stochastic gradient
                sample = randsample(lx, 1);
                grad_x = 2*n*transpose(FA(sample,:))*(FA(sample,:)*x-F(sample,:)*y) - n*(2*transpose(FA(sample,:))*(FA(sample,:)*snap_x-F(sample,:)*snap_y)) + snap_grad.g1;
                x = x - tau(1)*grad_x;
                sample = randsample(lx, 1);
                grad_y = 2*n*transpose(F(sample,:))*(F(sample,:)*y-FA(sample,:)*x)-2*n*lambda*transpose(F(sample,:))*(F(sample,:)*y - Fy0(sample)) - n*(2*transpose(F(sample,:))*(F(sample,:)*snap_y-FA(sample,:)*snap_x)-2*lambda*transpose(F(sample,:))*(F(sample,:)*snap_y - Fy0(sample))) + snap_grad.g2;
                y = y + tau(2)*grad_y;
                
                % compute number of foc and error 
                
                %pos = (k-1)*T*N + (T-1)*N + j;
                %pos = pos+1;
                %if error_type == 2
                %    max = max_oracle(FA, F, lambda, Fy0, x, y);
                %    errors(pos+1) = max - g_opt + max - obj_func(FA, F, lambda, Fy0, x, y);
                %else 
                %    errors(pos+1) = (norm(x - x_opt))^2 + (norm(y - y_opt))^2;
                %end
                %foc(pos+1) = foc(pos) + 1;
                %if j == 1
                %    foc(pos+1) = foc(pos+1) + lx;
                %end
                
                if j == N
                    pos = pos +1;
                    foc(pos) = foc(pos-1) + lx +N;
                    if error_type == 2
                       max = max_oracle(FA, F, lambda, Fy0, x, y);
                       errors(pos) = max - g_opt + max - obj_func(FA, F, lambda, Fy0, x, y);
                    else 
                       errors(pos) = (norm(x - x_opt))^2 + (norm(y - y_opt))^2;
                    end
                end
               
                
            end
           
        end
    end
    
    
 %% Outputs
    
    Result.x = x;
    Result.y = y;
    Result.errors = errors;
    Result.foc = foc;


end