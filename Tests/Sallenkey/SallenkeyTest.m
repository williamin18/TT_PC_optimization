load('Tests/Sallenkey/samples.mat')

n = 4;
d = 7;
N = n*ones(d,1);
r = 3;

n_data = length(b_samples);
n_training = 500;

A_train = A_samples(:,:,1:n_training);
b_train = b_samples(1:n_training);

A_test = A_samples(:,:,n_training+1:n_data);
b_test = b_samples(n_training+1:n_data);

x = TTrand(N,r);
% x = TT_simple_GD_adam(A_train,b_train,x,1,r,1e-3,50,A_test,b_test);

tic
x = TT_ALS(A_train,b_train,x,n_training,r,1e-3,1000,A_test,b_test);
toc