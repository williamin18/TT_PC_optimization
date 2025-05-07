clear variables
clc

load('Tests/Sallenkey/Completion_Test/sample_completion_vouts.mat');
load('Tests/Sallenkey/Completion_Test/sample_xi.mat');
out_data = load('Tests/Sallenkey/Sallenkey_8par.mat');
out_samples = out_data.samples;
out_true = out_data.vouts;

out_samples = abs(out_samples);
vouts_train = abs(vouts_train);

[~,d] = size(sample_xi);
m = 3;
N = (m+1)*ones(d,1);
r = 3;
y = TTrand(N,r);
y{d} = y{d}/norm(y{d},'fro');

[vouts_predicted,PC_coefficients,training_err,test_err] = ...
    pc_collocation_tensor_completion(sample_indices,vouts_train,y,out_samples,m,'Hermite','TT-Riemannian');