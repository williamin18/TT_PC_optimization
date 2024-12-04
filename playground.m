close all

xi = 2*rand(100,10000)-1;

x = xi/2+1.5;
x(20,:) = xi(20,:)+2;

y = uq_many_inputs_model(x');
mean(y)
std(y)
figure(1)
histogram(y)