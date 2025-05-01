function [b_predict,PC_coefficients,training_err,test_err] = ...
    pc_collocation_tensor_completion(sample_indices,b_train,x,A_predict,order,polynomial,method)
%PC_COLLOCATION_TENSOR Summary of this function goes here
%   Detailed explanation goes here
switch method
    case "TT-Riemannian"
        f = @TT_Riemannian_completion;
    otherwise
        err('Unsupported method')
end
    


[n_samples,d] = size(sample_indices);



n_train = round(0.9*n_samples);
training_indices = sample_indices(1:n_train,:);
training_samples = b_train(1:n_train,:);
test_indices = sample_indices(n_train+1:end,:);
test_samples = b_train(n_train+1:end,:);

%preconditioning
for i = 1:d
    for j = 1:order
        training_samples{i}(:,j+1) = training_samples{i}(:,j+1)*lambda1^j;
        test_samples{i}(:,j+1) = test_samples{i}(:,j+1)*lambda1^j;
        predict_samples{i}(:,j+1) = predict_samples{i}(:,j+1)*lambda1^j;
    end
end


PC_coefficients = cell(n_b,1);
training_err = zeros(n_b,1);
test_err = zeros(n_b,1);
n_iterations = zeros(n_b,1);

for i = 1:n_b
    [x,training_err(i),test_err(i),n_iterations(i)] = f(training_samples,training_out(:,i),x,4,1e-4,2000,test_samples,test_out(:,i),0.1);
    [training_err(i) test_err(i) n_iterations(i)] 
    PC_coefficients{i} = x;
end



[n_predict_samples,~] = size(A_predict);
b_predict = zeros(n_predict_samples,n_b);
for i = 1:n_b
    b_predict(:,i) = multi_r1_times_TT(predict_samples,PC_coefficients{i});
end
end

