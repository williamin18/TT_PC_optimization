clear variables
load('Tests/e2_double_tlnet/e2_data.mat')
[~,d] = size(training_samples);

training_samples = training_samples(1:660,:);
vouts_train = vouts_train(1:660,:);
m = 3;
r = 8;
x = cell(d,1);
for i = 1:d
    x{i} = zeros(m+1,1);
end


tic
 [vouts_GD,PC_coefficients_GD,training_err,test_err,n_iterations] = ...
     pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-SGD',0.3,0.2,r,40,10);
norm(vouts_GD-vouts,"fro")/norm(vouts,"fro")
toc