clear variables

load('Tests/Sallenkey/Sallenkey_8par.mat')

[n_samples,d] = size(training_samples);
m = 3;

% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts)/norm(vouts);

N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,r);

x{1}(1) = vouts_train(1,1);
for i = 1:d
    x{i}(1)=1;
end

vouts_TT = pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-Newton',0.3);
norm(vouts_TT-vouts)/norm(vouts)