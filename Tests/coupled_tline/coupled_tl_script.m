clear variables

load('Tests/coupled_tline/coupled_tl_28par_10std.mat')

[n_samples,d] = size(training_samples);
training_samples = training_samples(1:660,:);
vouts_train = vouts_train(1:660,:);
m = 3;

% tic
% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts,"fro")/norm(vouts,"fro")
% toc

N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,r);
x{1}(1) = vouts_train(1,1);
for i = 1:d
    x{i}(1)=1;
end

% x = cell(d,1);
% for i = 1:d
%     x{i} = zeros(m+1,1);
% end
tic
 [vouts_TT,PC_coefficients,training_err,test_err,n_iterations] = ...
     pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-Newton',0.3,0.2,3);
toc
norm(vouts_TT-vouts,'fro')/norm(vouts,'fro')

