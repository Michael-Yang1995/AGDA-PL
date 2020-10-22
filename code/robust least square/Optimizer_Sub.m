function out = Optimizer_Sub(FA, F, lambda, Fy0, x, y, param, method)
    if method == "SGDA"
        tau = param.stepsize;
        g = grad_func(FA, F, lambda, Fy0, x, y, 0, 'det');
        x = x - tau(1) * g.g1;
        y = y + tau(2) * g.g2;
    elseif method == "AGDA"
        tau = param.stepsize;
        g = grad_func(FA, F, lambda, Fy0, x, y, 1, 'det');
        x = x - tau(1) * g.g1;
        g = grad_func(FA, F, lambda, Fy0, x, y, 2,  'det');
        y = y + tau(2) * g.g2;
    elseif method == "SAGDA"
        tau = param.stepsize;
        g = grad_func(FA, F, lambda, Fy0, x, y, 1,'stoc');
        x = x - tau(1) * g.g1;
        g = grad_func(FA, F, lambda, Fy0, x, y, 2,  'stoc');
        y = y + tau(2) * g.g2;
    end
    
  
    
    
    
    out.x = x;
    out.y = y;
end