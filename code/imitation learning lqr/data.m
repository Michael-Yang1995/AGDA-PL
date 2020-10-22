 clear
 clc
% 
n=20;
m=10;
 %% generate A and B
 rng(30);
 B = randn(n,m);
 U = randn(n,n);
 U = orth(abs(U));
 v = rand(n,1)*0.5;
 A = U*diag(v)/(U); % A has eigenvalues less than 0.5
% 


%% generate Qe and Re
bd1 = 0.2; % feasible set for Q and R
bd2 = 200;

%U = abs(randn(n,n));  % generate Qe
%Qe = nearPD(U, 0.8, 1.2);
% 
 U = orth(randn(n,n));
v = rand(n,1)*1.2+5; % Qe has eigvalues between 0.8^2 and 1.2^2
Qe = transpose(U)*diag(v)*U;
% 
% %U = abs(randn(m,m));  % generate Re
% %Re =  nearPD(U, 0.8, 1.2);
 U = orth(randn(m,m));
 v = rand(m,1)*1.2 +5; % Re has eigvalues between 0.8^2 and 1.2^2
 Re = transpose(U)*diag(v)*U;



%% Q0 and R0
Q0 = Qe;   %eye(n,n);
R0 = Re;   %eye(m,m);


lambda = 1;


