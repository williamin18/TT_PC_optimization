clear variables

load('Tests/Sallenkey/Sallenkey_8par.mat')
m = 3;
[n_samples,d] = size(training_samples);

% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts)/norm(vouts);
N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,3);
vouts_TT = pc_collocation_tensor_optimization(training_samples,vouts_train,x,samples,m,'Hermite','TT-ALS');
norm(vouts_TT-vouts)/norm(vouts)