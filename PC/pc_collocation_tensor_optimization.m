function [b_predict,PC_coefficients] = pc_collocation_tensor_optimization(A_train,b_train,x,A_predict,order,polynomial,method)
%PC_COLLOCATION_TENSOR Summary of this function goes here
%   Detailed explanation goes here
switch method
    case "TT-ALS"
        f = @TT_ALS;
    case "TT-Newton"
        f = @TT_GD;
    otherwise
        err('Unsupported optimization type')
end


[n_samples,~] = size(A_train);
n_train = round(0.9*n_samples);
training_samples = genPolynomialSamplesTensor(A_train(1:n_train,:),order,polynomial);
test_samples = genPolynomialSamplesTensor(A_train(n_train+1:end,:),order,polynomial);
training_out = b_train(1:n_train,:);
test_out =  b_train(n_train+1:end,:);

[~,n_b] = size(b_train);
PC_coefficients = cell(n_b,1);
for i = 1:n_b
    x = f(training_samples,training_out(:,i),x,4,1e-4,2000,test_samples,test_out(:,i));
    PC_coefficients{i} = x;
end

predict_samples = genPolynomialSamplesTensor(A_predict,order,polynomial);
for i = 1:n_b
    b_predict(i,:) = multi_r1_times_TT(predict_samples,PC_coefficients{i});
end
end

