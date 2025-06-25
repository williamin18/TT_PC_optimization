clear variables

load('Tests/TlineNet/TlineNet_29par.mat')
training_samples = training_samples(1:667,:);
vouts_train = vouts_train(1:667,:);
[n_samples,d] = size(training_samples);
m = 3;

% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts)/norm(vouts);

N = (m+1)*ones(d,1);
r = 3;
% x = TTrand(N,r);
% x = TTorthogonalizeLR(x);
% x{d} = x{d}/norm( x{d},'fro');
% x{1}(1) = vouts_train(1,1);
% for i = 1:d
%     x{i}(1)=1;
% end

x = cell(d,1);
for i = 1:d
    x{i} = zeros(4,1);
end

tic
 % [vouts_TT,PC_coefficients,training_err,test_err,n_iterations] = ...
 %     pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-ALS',0.3,0.2,3);


 [vouts_TT,PC_coefficients,training_err,test_err,n_iterations] = ...
     pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-SGD',0.3,0.2,8);
toc
norm(vouts_TT-vouts,'fro')/norm(vouts,'fro')