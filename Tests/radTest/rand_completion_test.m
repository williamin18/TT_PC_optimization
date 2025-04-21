n = 4;
d = 10;
N = n*ones(d,1);
r = 3;
n_samples = 3600;
n_test_samples = 100;



x_true = TTrand(N,3);
x_true = TTorthogonalizeLR(x_true);
x_true{d} = x_true{d}/norm( x_true{d},'fro');

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


x = TTrand(N,3);
x = TTorthogonalizeLR(x);
x{d} = x{d}/norm( x{d},'fro');




% comparsion_samples = load('Completion_comparsion_samples.mat');
A = comparsion_samples.A;
A_test = comparsion_samples.A_test;
x_true = comparsion_samples.x_true;
% x = comparsion_samples.x;


x_sample = multi_r1_times_TT(A,x_true);
x_test_sample =  multi_r1_times_TT(A_test,x_true);



% load('Completion_comparsion_samples.mat');

[x,training_err,test_err,epoch] = TT_Riemannian_completion(A,x_sample,x,r,1e-4,2000,A_test,x_test_sample,0);
TTnorm(TTaxby(1,x,-1,x_true))/TTnorm(x_true)