
function obj = obj_func(FA, F, lambda, Fy0, x, y)
    obj = (norm(FA*x-F*y,2))^2 - lambda*(norm(F*y-Fy0,2))^2;
end