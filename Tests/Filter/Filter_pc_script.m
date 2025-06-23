clear variables
freq = linspace(0.5e9,14e9,101);
training_data = load('Tests/Filter/data_training.mat');
mc_data = load('Tests/Filter/data_training.mat');
S21_train = (abs(training_data.S21));
training_samples = (training_data.samples);
S21_mc = (abs(mc_data.S21));
mc_samples = (mc_data.samples);

s = log(S11_train./(1-S21_train));
s_mc = log(S21_mc./(1-S21_mc));
pdf_freq_idx = 80;
f = figure(5);
Hmc = histogram( s(:,pdf_freq_idx) ,50,'Normalization','pdf', 'DisplayStyle','bar', 'FaceColor',[0.7 0.7 0.7]);

[n_samples,d] = size(training_samples);
m = 4;

% vouts_total = pc_collocation_total(training_samples,vouts_train,samples,m ,'Hermite');
% norm(vouts_total-vouts)/norm(vouts);

N = (m+1)*ones(d,1);
r = 3;
x = TTrand(N,r);
x = TTorthogonalizeLR(x);
x{d} = x{d}/norm( x{d},'fro');
x{1}(1) = S11_train(1,1);
for i = 1:d
    x{i}(1)=1;
end

tic
% [y_out,PC_coefficients] = pc_collocation_total(training_samples,S21_train,mc_data.samples,m,'Legendre');
[y_out,PC_coefficients,training_err,test_err,n_iterations] = ...
    pc_collocation_tensor_optimization(training_samples,s,x,mc_samples,m,'Hermite','TT-Newton',0.3,0.2,r);
toc
