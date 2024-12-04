xi = 2*rand(100,10000)-1;
y = zeros(10000,1);

for i = 1:1e4
    y(i) = SynFun100d(xi(:,i));
end

mean(y)
std(y)
figure(1)
histogram(y)