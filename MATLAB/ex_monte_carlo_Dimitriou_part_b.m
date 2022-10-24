% MONTE CARLO - DIMITRIOU ELEFTHERIOS - part - b - no for loop
clear all
clc

tic
rng('default')

sum=0;
n=10000000; % n = 10^4, 10^5, 10^6

x = -1 + 2*rand(n,1);
y = -1 + 2*rand(n,1);
z = -1 + 2*rand(n,1);
r = x.^2 + y.^2 + z.^2;

for i=1:n
    if r(i) <= 1
    sum = sum + 1;
    end
end


per = sum/n;
PI = 6*per;
st_d = 6*sqrt(var(r))/sqrt(n);
fprintf('The approximation of pi is %.10f with n = %d \n', PI, n)
fprintf('The standard deviation is s = %.10f \n', st_d)
interval = [PI - st_d, PI + st_d];
fprintf('The real value of pi belongs to the interval [%f %f]', interval)
fprintf('\n The real value of pi is pi = %.10f \n', pi)
toc   