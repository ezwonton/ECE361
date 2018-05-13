clear, clc, warning('off'), close all

data = sort(xlsread('HW6_data_shankar_Spring(1).xlsx', 1, 'BZ:BZ'));
[p, i] = ksdensity(data, 'NumPoints', 200);
figure;
grid on;
hold on;
plot(i, p, 'k*');

% Normal Stuff
n = fitdist(data, 'Normal');
norm = pdf('Normal', i, n.mu, n.sigma);
plot(i, norm, '--');
[h,p,stats] = chi2gof(data, 'CDF', n)

% Laplace stuff
a = mean(data);
b = sqrt(var(data)/2);

lap_pdf = (1/(2*b))*exp(-abs(i-a)/b);
lap_cdf = (i>a).*(1-0.5*exp(-(i-a)/b)) + (i<a).*(0.5*exp(-(i-a)/b));
plot(i, lap_pdf, '.');
legend("data, (ksDensity)", "LaPlace", "Norm")

[Y,E] = discretize(data,10);
bin = [unique(Y) arrayfun(@(x) length(find(Y == x)), unique(Y))];
bin_cdf = bin(:,2)/length(data);
