function g = grad_func(FA, F, lambda, Fy0, x, y, part, stoc)
    [lx,ly] = size(FA);
    g1 = zeros(ly, 1);
    g2 = zeros(lx, 1);
    
    if stoc == "det"
      if part == 0
          g1 = 2*transpose(FA)*(FA*x-F*y);
          g2 = 2*F.'*(F*y-FA*x)-2*lambda*F.'*(F*y - Fy0);
      elseif part == 1
          g1 = 2*transpose(FA)*(FA*x-F*y);
      elseif part == 2
          g2 = 2*F.'*(F*y-FA*x)-2*lambda*F.'*(F*y - Fy0);
      end
    elseif stoc == "stoc"
      if part == 0
          sample = randsample(lx,1);
          g1 = 2*transpose(FA(sample,:))*(FA(sample,:)*x-F(sample,:)*y);
          g2 = 2*transpose(F(sample,:))*(F(sample,:)*y-FA(sample,:)*x)-2*lambda*transpose(F(sample,:))*(F(sample,:)*y - Fy0(sample));
      elseif part == 1
          sample = randsample(lx,1);
          g1 = 2*transpose(FA(sample,:))*(FA(sample,:)*x-F(sample,:)*y);
      elseif part == 2
          sample = randsample(lx,1);
          g2 = 2*transpose(F(sample,:))*(F(sample,:)*y-FA(sample,:)*x)-2*lambda*transpose(F(sample,:))*(F(sample,:)*y - Fy0(sample));
      end
    end
      
    g.g1 = g1;
    g.g2 = g2;
end