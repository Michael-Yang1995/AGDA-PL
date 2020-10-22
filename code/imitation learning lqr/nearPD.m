function B = nearPD(A, bd1, bd2)
    A = .5*(A+A.');
    [V,D] = eig(A);
    d = diag(D);
    d = max(d, bd1);
    d = min(d, bd2);
    B = V*diag(d)*V.';
end