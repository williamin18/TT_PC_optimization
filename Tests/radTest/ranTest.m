n = 4;
d = 10;
N = n*ones(d,1);
r = 3;
n_samples = 2000;
n_test_samples = 100;


% X_true = TTrandComplex(N,3);
X_true = TTrand(N,3);
X_true = TTorthogonalizeLR(X_true);
X_true{d} = X_true{d}/norm( X_true{d},'fro');

A = cell(d,1);
A_test = cell(d,1);
for i =1: d
    A{i} = randn(n_samples,n);
    A_test{i} = rand(n_test_samples,n);
end


b = multi_r1_times_TT(A,X_true);
b_test = multi_r1_times_TT(A_test,X_true);

x = TTrand(N,3);
x = TTorthogonalizeLR(x);
x{d} = x{d}/norm( x{d},'fro');
% x = X_true;
% x{1} = x{1}-0.05*rand(n,r);
tic
[x,training_err,test_err,epoch] = TT_Newton_GD(A,b,x,r,1e-3,2000,A_test,b_test,0)
% % x = TT_simple_GD(A,b,x,1,r,1e-3,50,A_test,b_test)
%x = TT_ALS_old(A,b,x,n_samples,r,1e-3,2000,A_test,b_test)
toc
TTnorm(TTaxby(1,x,-1,X_true))/TTnorm(X_true)