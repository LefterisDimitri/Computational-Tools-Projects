% RANDOM VARIABLES - DIMITRIOU ELEFTHRIOS
clear all
clc

rng('default')
sum = 0;   
sum1 = 0;

figure(1)
subplot(2,1,1)
for m=1:1000
    measure(m) = randn(1);
    plot(m,measure(m),'b.')
    hold on
    sum = sum + measure(m);
end
title('Generated numbers from normal distribution with mu = 0 and sigma = 1','Fontsize',10)
xlabel('m','Interpreter','latex','Fontsize',11)
ylabel('Measurment','Interpreter','latex','Fontsize',11)

figure(2)
histogram(measure)
title('Histogram of measurments','Fontsize',10)

mean_measure_0 = sum/m;
mean_measure_1 = mean(measure);

for m=1:1000
    sum1 = sum1 + (measure(m) - mean_measure_0)^2;
end
st_d_measure_0 = sqrt(sum1/(m-1));
st_d_measure_1 = sqrt(var(measure));
fprintf('analytically mean value = %f and mean value from MATLAB = %f \n', mean_measure_0, mean_measure_1)
fprintf('analytically standard deviation = %f and standard deviation from MATLAB = %f \n', st_d_measure_0, st_d_measure_1)

km = 1000; % km = m

%rng('default')
figure(1)
subplot(2,1,2)
for m=1:1000
    measure_1 = randn(km,1);
    measure_mean_exper(m) = mean(measure_1);
    plot(m,measure_mean_exper(m),'r.')
    hold on
end
mean_data = mean(measure_mean_exper);
st_d_data = sqrt(var(measure_mean_exper));
fprintf('Mean value of mean values from MATLAB = %f \n', mean_data)
fprintf('Theoretically standard deviation of mean value of mean values = %f standard deviation of mean values from MATLAB = %f \n',st_d_measure_0/sqrt(m), st_d_data)

title('Measurments of #km experiments','Fontsize',10)
xlabel('m','Interpreter','latex','Fontsize',11)
ylabel('Measurment','Interpreter','latex','Fontsize',11)

