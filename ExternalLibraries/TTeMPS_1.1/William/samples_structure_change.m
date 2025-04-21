clear variables
load('TT_riemannian_samples.mat')

d = A.order;
r = A.rank;
m = A.size;
n_samples = length(A_Omega);
n_test_samples = length(A_Gamma);

x_true = cell(d,1);
x =  cell(d,1);
for i = 1:d
    x_true{i} = reshape(A{i},[r(i)*m(i) r(i+1)]);
    x{i} = reshape(X0{i},[r(i)*m(i) r(i+1)]);
end


x_sample = A_Omega;
x_test_sample = A_Gamma;

A =  cell(d,1);
A_test = cell(d,1);
for i =1: d
    A{i} = zeros(m(i),n_samples);
    idx = Omega(:,i)';
    idx = idx + m(i)*(0:n_samples-1);
    A{i}(idx) = 1;
    A{i} = A{i}';

    A_test{i} = zeros(m(i),n_test_samples);
    idx = Gamma(:,i)';
    idx = idx + m(i)*(0:n_test_samples-1);
    A_test{i}(idx) = 1;
    A_test{i} = A_test{i}';
end