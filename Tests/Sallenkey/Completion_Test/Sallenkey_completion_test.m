clear variables
clc

load('Tests/Sallenkey/Completion_Test/sample_completion_vouts.mat');
load('Tests/Sallenkey/Completion_Test/sample_xi.mat');
out_data = load('Tests/Sallenkey/Sallenkey_8par.mat');
out_samples = out_data.samples;
out_true = out_data.vouts;

load('Tests/Sallenkey/Completion_Test/r1_samples.mat');
% out_samples = abs(out_samples);
% vouts_train = abs(vouts_train);






[~,d] = size(sample_xi);
[n_samples,n_outs] = size(vouts_train);
m = 3;
N = (m+1)*ones(d,1);
r = 3;

y_init = cell(n_outs,1);
for i = 1:n_outs
    y_init{i} = formRank1Tensor(vout_ref(i),vouts_r1(:,i),m,d);
end


% [vouts_predicted,PC_coefficients,training_err,test_err] = ...
%     pc_collocation_tensor_completion(sample_idx,vouts_train(:,1),y_init(1),out_samples,m,'Hermite','TT-Riemannian');

tic
[vouts_predicted,PC_coefficients,training_err,test_err] = ...
    pc_collocation_tensor_completion(sample_idx,vouts_train,y_init,out_samples,m,'Hermite','TT-Riemannian');
toc
norm(vouts_predicted-out_true)/norm(out_true)