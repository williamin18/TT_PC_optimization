n = 4;
d = 10;
N = n*ones(d,1);
r = 3;
n_samples = 2000;
n_test_samples = 100;


X_true = TTrand(N,1);
A = 3*rand(d,n,n_samples);
A_test = 3*rand(d,n,n_test_samples);



b = multi_r1_times_TT(A,X_true);
b_test = multi_r1_times_TT(A_test,X_true);

x = TTrand(N,r);
% x = X_true;
%x{1} = x{1}-0.05*rand(n,r);
%x = TT_GD(A,b,x,16,r,1e-3,20,A_test,b_test)
x = TT_ALS(A,b,x,1,r,1e-3,50,A_test,b_test)