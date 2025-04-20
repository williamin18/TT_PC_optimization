n = 4;
d = 10;
N = n*ones(d,1);
r = 3;
n_samples = 2000;
n_test_samples = 100;



X_true = TTrand(N,3);
A = cell(d,1);
A_test = cell(d,1);
for i =1: d
    A{i} = zeros(n,n_samples);
    idx = randi(n,1,n_samples);
    idx = idx + n*(0:n_samples-1);
    A{i}(idx) = 1;
    A{i} = A{i}';

    A_test{i} = zeros(n,n_test_samples);
    idx = randi(n,1,n_test_samples);
    idx = idx + n*(0:n_test_samples-1);
    A_test{i}(idx) = 1;
    A_test{i} = A_test{i}';
end



x_samples = multi_r1_times_TT(A,X_true);
x_test_sample =  multi_r1_times_TT(A_test,X_true);

x = TTrand(N,r);


x = TT_Riemannian_completion(A,x_samples,x,r,1e-4,2000,A_test,x_test_sample,0);
TTnorm(TTaxby(1,x,-1,X_true))/TTnorm(X_true)