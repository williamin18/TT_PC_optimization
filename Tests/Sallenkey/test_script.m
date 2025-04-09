clear variables

load('Tests/Sallenkey/Sallenkey_8par.mat')
m = 3;
vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
norm(vouts_total-vouts)/norm(vouts);

vouts_TT = pc_collocation_tensor_optimization(training_samples,vouts_train,samples,m,'Hermite','TT-Newton');