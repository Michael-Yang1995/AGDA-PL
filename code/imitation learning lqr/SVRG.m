function Result = SVRG(A, B, Qe, Re, Q0, R0, lambda, param)
 %% initialization
    [lx,ly] = size(B);
    X = param.X;  % X is samples, each row is a sample
    [n_sample, ~] = size(X);
    iter = param.iter;
    tau = param.stepsize;
    K = param.K;
    Q = param.Q;
    R = param.R;
    N = param.N;    
    T = param.T;
    bd1 = param.bd1;
    bd2 = param.bd2;
    rollout = param.rollout;
    errors = zeros(1,iter);
    foc = zeros(1, iter);
    
    %% compute Ke and correlation matrix sigma_e
    [Ke, ~, ~] = dlqr(A, B, Qe, Re); % compute optimal control Ke corresponding to Qe, Re
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
    
    
    errors(1) = norm(K-Ke, 'fro')^2 + norm(Q - Qe, 'fro')^2 + norm(R - Re, 'fro')^2;
    
    
    %% iterations
    
    
    for k = 1:iter
        for t = 1:T
            % store gradients of snapshots
            snap_grad = grad_func(A, B, Q, R,  Q0, R0, lambda, K, Ke, sigma_e, X, rollout, 0, 'det');
            snap_K = K;
            snap_Q = Q;
            snap_R = R;
            for j = 1:N
                % compute new x and y using stochastic gradient
                X_sample = X(randsample(n_sample, 1),:);
                g_sample = grad_func(A, B, Q, R, Q0, R0, lambda, K, Ke, sigma_e, X_sample, rollout, 1, 'stoc');
                g_snap_sample = grad_func(A, B, snap_Q, snap_R, Q0, R0, lambda, snap_K, Ke, sigma_e, X_sample, rollout, 1, 'stoc');
                g_K = g_sample.g1 - g_snap_sample.g1 + snap_grad.g1;
                K = K - tau(1)*g_K;
                
                X_sample = X(randsample(n_sample, 1),:);
                g_sample = grad_func(A, B, Q, R, Q0, R0, lambda, K, Ke, sigma_e, X_sample, rollout, 2, 'stoc');
                g_snap_sample = grad_func(A, B, snap_Q, snap_R, Q0, R0,  lambda, snap_K, Ke, sigma_e, X_sample, rollout, 2, 'stoc');
                g_Q = g_sample.g2 - g_snap_sample.g2 + snap_grad.g2;
                Q = Q + tau(2)*g_Q;
                %Q = nearPD(Q, bd1, bd2);
                g_R = g_sample.g3 - g_snap_sample.g3 + snap_grad.g3;
                R = R + tau(2)*g_R;
                %R = nearPD(R, bd1, bd2);
                
                if j == N
                    fprintf('The norm square of K is %e \n', norm(g_K,'fro'))
                    fprintf('The norm square of Q is %e \n', norm(g_Q,'fro'));
                    fprintf('The norm square of R is %e \n', norm(g_R,'fro'));
                    %compute error
                    k
                    errors(k+1) = norm(K-Ke, 'fro')^2 + norm(Q - Qe, 'fro')^2 + norm(R - Re, 'fro')^2;
                    foc(k+1) = foc(k)+N+n_sample;
                end
                
            end
           
        end
    end
    
    
    
    
    Result.K = K;
    Result.Q = Q;
    Result.R = R;
    Result.errors = errors;
    Result.foc = foc;


end