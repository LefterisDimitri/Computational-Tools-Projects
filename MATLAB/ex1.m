% MONTE CARLO - DIMITRIOU ELEFTHERIOS
clear all
clc

tic
rng('default')

sum=0;
for n=1:100000 % n = 10^4, 10^5, 10^6, 10^7, 10^8 
    x(n) = -1 + 2*rand(1);
    y(n) = -1 + 2*rand(1);
    z(n) = -1 + 2*rand(1);
    r(n) = x(n)^2 + y(n)^2 + z(n)^2;
    if r(n) <= 1
        sum = sum +1;
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
% scatter3(x,y,z)
% hold on
% 
% [x,y,z]=sphere;
% Sphere3D=surf(x,y,z);
toc  
