function Result = Optimizer(A, F, lambda, y0, param, method, x_opt, y_opt, g_opt, error_type)
    %A is a matrix in the least square, F is the matrix associated with M norm
    %lambda is the constant for regularizer
    %y0 is observed vector in least square
    %param includes initial points and tolerance
    %x_opt, y_opt, g_opt are optimal values used to measure errors
    %error_type = 1 or 2. 
    
    %% Parameters
    x_0 = param.x_0;
    y_0 = param.y_0;
    iter = param.iter;    
    [n,~] = size(y_0); 

    %% Initialization
    x = x_0;
    y = y_0;  
    if strcmp(method, "SAGDA") == 1
        iter1 = ceil(iter/n) + 2;
    else 
        iter1 = iter;
    end
    errors = zeros(1,iter1);
    FA = F*A;
    Fy0 = F*y0;
    
    
    
    %% find error for x_0;
    if error_type == 2 
        max = max_oracle(FA, F, lambda, Fy0, x, y);
        errors(1) = max - g_opt + max - obj_func(FA, F, lambda, Fy0, x, y);
    elseif error_type == 1
        errors(1) = (norm(x - x_opt))^2 + (norm(y - y_opt))^2;
    end
    
    %% Iterations for deterministic algorithm
    if strcmp(method, "AGDA") == 1
        for i = 1:(iter-1)
            %% Update
            outputs = Optimizer_Sub(FA, F, lambda, Fy0, x, y, param, method);
            x = outputs.x;
            y = outputs.y;

        
            %% Record errors and print out
            if error_type == 2
                max = max_oracle(FA, F, lambda, Fy0, x, y);
                errors(i+1) = max - g_opt + max - obj_func(FA, F, lambda, Fy0, x, y);
            else 
                errors(i+1) = (norm(x - x_opt))^2 + (norm(y - y_opt))^2;
            end
            if mod(iter,10) == 0
                fprintf('The error is %e \n', errors(i+1));
            end
        end
        
    elseif strcmp(method, "SAGDA") == 1
        for i = 1:(iter-1)
            %% Update
            outputs = Optimizer_Sub(FA, F, lambda, Fy0, x, y, param, method);
            x = outputs.x;
            y = outputs.y;

        
            %% Record errors
            if mod(i, n) == 0
               if error_type == 2
                   max = max_oracle(FA, F, lambda, Fy0, x, y);
                   errors(i/n+1) = max - g_opt + max - obj_func(FA, F, lambda, Fy0, x, y);
               else 
                   errors(i/n+1) = (norm(x - x_opt))^2 + (norm(y - y_opt))^2;
               end
               fprintf('The error is %e \n', errors(i/n+1))
            end
            %if mod(i,100) == 0
            %    fprintf('The error is %e \n', errors(i+1));
            %end
        end
    end
    
    
    
    
    %% Outputs
    
    Result.x = x;
    Result.y = y;
    Result.errors = errors;
end
