clear
clc

%% dimension of data
n = 546;
m = 8;

%% generate A, x* and y0
A = readtable('qsar_aquatic_toxicity.csv');
A = A{:,:};
y0 = A(:,end);
A = A(:,1:end-1);
F = eye(n,n);

mean_A = mean(A);
var_A = var(A);  

A = (A-mean(A))./sqrt(var_A); 





lambda = 2;

