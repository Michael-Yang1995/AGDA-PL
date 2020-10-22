function Result = Optimizer1(A, F, lambda, y0, param, method)
    %% Parameters
    x_0 = param.x_0;
    y_0 = param.y_0;
    tol = param.tol;
    iter_max = 50000;
    
    
    %% Initialization
    x = x_0;
    y = y_0;
    norm_square = 10;
    iter = 0;
    
    FA = F*A;
    Fy0 = F*y0;

    
    %% Iterations for deterministic algorithm
    if strcmp(method, "SGDA")+strcmp(method, "AGDA") == 1
        while norm_square > tol
            iter = iter +1;
            %% Update
            x_old = x;
            y_old = y;
            outputs = Optimizer_Sub(FA, F, lambda, Fy0, x, y, param, method);
            x = outputs.x;
            y = outputs.y;
            dif = norm(x-x_old,2) + norm(y-y_old,2);
            g = grad_func(FA, F, lambda, Fy0, x, y, 0,'det');
            norm_square = norm(g.g1)+norm(g.g2);
           
            
            if mod(iter,10) == 0
                fprintf('The difference is %e \n', dif);
                fprintf('The norm square is %e \n', norm_square);
            end
            %if dif > 1000
            %    fprintf('It diverges \n');
            %    break
            %end
            if iter > iter_max
                fprintf('It does not converge within 30000 iters \n')
                break
            end
            if norm_square < tol
                fprintf('The total iteration is %e \n', iter);
            end
        end
    end
    
    
    
    
    %% Outputs
    
    Result.x = x;
    Result.y = y;
end
