function [b_predict,PC_coefficients,training_err,test_err] = ...
    pc_collocation_tensor_completion(sample_indices,b_train,y,A_predict,order,polynomial,method)
%y is TT-estimation of outputs, b is sample outputs,sample_indices are indices of b,
% A is polynomial samples(not used on training)
switch method
    case "TT-Riemannian"
        f = @TT_Riemannian_completion;
    otherwise
        err('Unsupported method')
end

switch polynomial
    case "Hermite"
        poly = @genHermite;
    case "Legendre"
        poly = @genLegendre;
    otherwise
        err('Unsupported polynomial')
end
    

%% initialization
[n_samples,d] = size(sample_indices);


n_train = round(0.9*n_samples);
training_indices = sample_indices(1:n_train,:);
training_samples = b_train(1:n_train,:);
test_indices = sample_indices(n_train+1:end,:);
test_samples = b_train(n_train+1:end,:);


training_sample_selector = cell(d,1);
test_sample_selector = cell(d,1);
for i = 1:d
    training_sample_selector{i} = zeros(order+1,n_train);
    idx = training_indices(:,i);
    idx = idx + (order+1)*(0:n_train-1);
    training_sample_selector{i}(idx) = 1;
    training_sample_selector{i} = training_sample_selector{i}';

    test_sample_selector{i} = zeros(order+1,n_samples-n_train);
    idx = test_indices(:,i);
    idx = idx + (order+1)*(0:n_samples-n_train-1);
    test_sample_selector{i}(idx) = 1;
    test_sample_selector{i} = training_sample_selector{i}';
end

%% TT outs estimation
[~,n_b] = size(b_train);
TT_outs = cell(n_b,1);
training_err = zeros(n_b,1);
test_err = zeros(n_b,1);
n_iterations = zeros(n_b,1);


for i = 1:n_b
    [y,training_err(i),test_err(i),n_iterations(i)] = f(training_sample_selector,training_samples(:,i),y,4,1e-4,2000,test_sample_selector,test_samples(:,i),0.1);
    [training_err(i) test_err(i) n_iterations(i)] 
    TT_outs{i} = y;
end

%% compute PC coefficients
PC_coefficients = cell(n_b,1);
polynomial_roots =  genPolynomialRoots(order+1,polynomial);
A_r1 = zeros(order+1,order+1);
for i = 1:order+1
    A_r1(:,i) = poly(polynomial_roots,i-1);
end
A_r1_inverse = A_r1^-1;
for i = 1:n_b
    PC_coefficients{i} = TT_outs{i};
    [~,m,r] = TTsizes(x);
    for j = 1:d
        Wj = PC_coefficients{i}{j};
        Wj = reshape(Wj,[r(j),m(j),r(j+1)]);
        Wj = reshape(permute(Wj, [2 1 3]),m(i),[]);
        Wj = A_r1_inverse*Wj;
        Wj = permute(reshape(Wj,m(i),r(i),r(i+1)),[2 1 3]);
        PC_coefficients{i}{j} = Wj;
    end
end

%% predict outputs
[n_predict_samples,~] = size(A_predict);
b_predict = zeros(n_predict_samples,n_b);
for i = 1:n_b
    b_predict(:,i) = multi_r1_times_TT(predict_samples,PC_coefficients{i});
end
end

