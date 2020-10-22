function Result = PolicyGrad(A, B, Q0, R0, Qe, Re, lambda, param)
    %% initialization
    [lx,ly] = size(B);
    X = param.X;  % X is samples, each row is a sample
    [n_sample, ~] = size(X);
    iter = param.iter;
    tau = param.stepsize;
    K = param.K;
    Q = param.Q;
    R = param.R;
    bd1 = param.bd1;
    bd2 = param.bd2;
    rollout = param.rollout;
    errors = zeros(1,iter);
    
    %% compute Ke and correlation matrix sigma_e
    [Ke, Pe, e] = dlqr(A, B, Qe, Re); % compute optimal control Ke corresponding to Qe, Re
    sigma_e = 0;
    for t = 1:n_sample
        x = X(t,:).';
        for i=1:rollout
            sigma_e = sigma_e + x*x.';
            u = -Ke*x;
            x = A*x + B*u;
        end
    end
    sigma_e = sigma_e/n_sample;
    
    
    % compute initial error
    errors(1) = norm(K-Ke, 'fro')^2 + norm(Q - Qe, 'fro')^2 + norm(R - Re, 'fro')^2;
    
    
    
    %% iterations
    for i=1:(iter-1)
        
        % update
        g = grad_func(A, B, Q, R, Q0, R0, lambda, K, Ke, sigma_e, X, rollout, 1, 'det');
        fprintf('The norm square of K is %e \n', norm(g.g1,'fro'));
        K = K - tau(1)*g.g1;
        g = grad_func(A, B, Q, R, Q0, R0, lambda, K, Ke, sigma_e, X, rollout, 2, 'det');
        Q = Q + tau(2)*g.g2;
        %Q = nearPD(Q, bd1, bd2);
        R = R + tau(2)*g.g3;
        %R = nearPD(R, bd1, bd2);
        if mod(iter,1) == 0
            fprintf('The norm square of Q is %e \n', norm(g.g2,'fro'));
            fprintf('The norm square of R is %e \n', norm(g.g3,'fro'));
        end
        
        %compute error
        i
        errors(i+1) = norm(K-Ke, 'fro')^2 + norm(Q - Qe, 'fro')^2 + norm(R - Re, 'fro')^2;
        
    end 
    
    
    
    Result.K = K;
    Result.Q = Q;
    Result.R = R;
    Result.errors = errors;
    
    
    
    


end



