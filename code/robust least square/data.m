clear
clc

%% dimension of data
n = 100;
m = 50;

%% generate A, x* and y0
rng(400);
lambda = 3;
F = eye(n,n);
A = random('normal',0,1,n,m);
x_star = random('normal', 0, 1, m,1);
y0 = A*x_star + random('normal', 0, 0.01, n,1);


%% define M-norm and lambda
p = floor(n*0.95);  % rank = n*0.9
U = randn(n, n);  
v = rand(p,1)/10+0.9;   % simga has eigenvalues between (0.9,1.1)
sigma = diag([v; zeros(n-p,1)]);
U = orth(U);
F = transpose(U)*sigma*U;  % F = U*sigma*U^T
% M = F^tF is rank deficit and it has eigenvalues between (0.9^2, 1.1^2)


lambda = 3;