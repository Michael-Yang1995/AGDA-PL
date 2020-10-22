function g = grad_func(A, B, Q, R, Q0, R0, lambda, K, Ke, sigma_e, X, rollout, part, stoc)
    [lx,ly] = size(B);
    [n_sample, ~] = size(X);
    g1 = zeros(ly, lx);
    g2 = zeros(lx, lx);
    g3 = zeros(ly, ly);
    
    
    iter_max = 3000;
    
    %compute sigma
    sigma = 0;
    for t = 1:n_sample
       x = X(t,:).';
       for m = 1:rollout
           sigma = sigma + x*x.';
           u = -K*x;
           x = A*x + B*u;
       end
    end
    sigma = sigma/n_sample;
    
    %compute sigma_e if we use stochastic gradient
    if stoc == "stoc"
        sigma_e = 0;
        for t = 1:n_sample
           x = X(t,:).';
           for m = 1:rollout
               sigma_e = sigma_e + x*x.';
               u = -K*x;
               x = A*x + B*u;
           end
        end
        sigma_e = sigma/n_sample;
    end
    
    
    % solve value matrix P
    P = zeros(lx, lx);
    dif = 10;
    iter = 1;
    while dif > 10^(-8)
         iter = iter+1;
         P_old = P;
         P = Q + K.'*R*K+(A-B*K).'*P*(A-B*K);
         dif = norm(P-P_old,'fro');
         if iter > iter_max
             fprintf('It does not converge within 3000 iters \n')
             break
         end
    end 
    
    
    
    if stoc == "det"
      if part == 0
          g1 = 2*((R+B.'*P*B)*K-B.'*P*A)*sigma;
          g2 = sigma-sigma_e-2*lambda*(Q-Q0);
          g3 = K*sigma*K.' - Ke*sigma_e*Ke.' - 2*lambda*(R-R0);
      elseif part == 1
          g1 = 2*((R+B.'*P*B)*K-B.'*P*A)*sigma;
      elseif part == 2
          g2 = sigma-sigma_e-2*lambda*(Q-Q0);
          g3 = K*sigma*K.' - Ke*sigma_e*Ke.' - 2*lambda*(R-R0);
      end
    elseif stoc == "stoc"
      if part == 0
          g1 = 2*((R+B.'*P*B)*K-B.'*P*A)*sigma;
          g2 = sigma-sigma_e-2*lambda*(Q-Q0);
          g3 = K*sigma*K.' - Ke*sigma_e*Ke.' - 2*lambda*(R-R0);
      elseif part == 1
          g1 = 2*((R+B.'*P*B)*K-B.'*P*A)*sigma;
      elseif part == 2
          g2 = sigma-sigma_e-2*lambda*(Q-Q0);
          g3 = K*sigma*K.' - Ke*sigma_e*Ke.' - 2*lambda*(R-R0);
      end
    end
      
    g.g1 = g1;
    g.g2 = g2;
    g.g3 = g3;
end