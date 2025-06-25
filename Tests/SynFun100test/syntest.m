n_train = 1000;
xi_test = 2*rand(10000,100)-1;
xi_train = 2*lhsdesign(n_train,100)-1;

x = xi_test/2+1.5;
x(:,20) = xi_test(:,20)+2;

y = xi_train/2+1.5;
y(:,20) = xi_train(:,20)+2;

% out_test = zeros(10000,1);
% out_train = zeros(n_train,1);
% for i = 1:1e4
%     out_test(i) = SynFun100d(x(i,:)');
% end
% for i = 1:n_train
%     out_train(i) = SynFun100d(y(i,:)');
% end

out_test = uq_many_inputs_model(x);
out_train = uq_many_inputs_model(y);


d = 100;
m = 3;


% s = TTrand(N,r);
% s{1}(1) = out_train(1);
% for i = 1:d
%     s{i}(1)=1;
% end

s = cell(d,1);
for i = 1:d
    s{i} = zeros(4,1);
end

out_predict = pc_collocation_tensor_optimization(xi_train,out_train,s,xi_test,m,'Hermite','TT-SGD',0.3,0.2,20);

N = (m+1)*ones(d,1);
r = 3;
s = TTrand(N,r);
s{1}(1) = out_train(1);
for i = 1:d
    s{i}(1)=1;
end

out_predict2 = pc_collocation_tensor_optimization(xi_train,out_train,s,xi_test,m,'Hermite','TT-Newton',0.3,0.2,3);


norm(out_predict-out_test,'fro')/norm(out_test,'fro')
norm(out_predict2-out_test,'fro')/norm(out_test,'fro')

[mean(out_predict),mean(out_predict2), mean(out_test)]
[std(out_predict),std(out_predict2), std(out_test)]
figure(1)
histogram(out_test)
