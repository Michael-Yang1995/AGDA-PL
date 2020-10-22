function max = max_oracle(FA, F, lambda, Fy0, x, y)
    norm_square = 10;
    tau = 0.5;
    max_iter = 5000;
    iter = 0;
    while norm_square > 10^(-11)
        g = grad_func(FA, F, lambda, Fy0, x, y, 2, 'det');
        y = y + tau*g.g2;
        norm_square = norm(g.g2);
        iter = iter +1;
        if iter > max_iter
              fprintf('Max oracle does not converge within 5000 iters \n')
              break
        end
    end
    max = obj_func(FA, F, lambda, Fy0, x, y);
    %iter
end


%stepsize for data 1 is 0.5
%stepsize for data 2 is 0.25
%stepsize for data 3 is 1