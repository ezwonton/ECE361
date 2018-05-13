%% Eric Wan - ezw23@drexel.edu
clear, clc, warning('off'), close all

%% Question 10
%{
You are given a data set consisting of 200 entries. Determine the best fit
using chi square testing. Check for Nakagami, gamma, Weibull if the data
set is completely positive. If data set contains negative values, test for
normal or Laplacian. Laplacian is not a built in pdf in Matlab.  To get the
parameters of the Laplacian, simply find the mean and variance as indicated
in HW#4. [You will see that a file named HW4_data_shankar_Spring. You will
see your name (last   name only) at the top of the column. You are required
to use the data in that column].
%}

data = xlsread("HW6_data_shankar_Spring(1).xlsx", 1, "CF:CF");
datasort = sort(data);

[p, i] = ksdensity(datasort, 'NumPoints', 200);
figure;
grid on;
hold on;
plot(i, p, 'k*');

w = fitdist(datasort, 'Weibull'); % Weibull Probability Distribution
wei = pdf('Weibull', i, w.A, w.B);
[h,p,stats] = chi2gof(datasort,'CDF',w)
plot(i, wei, '--');

n = fitdist(datasort, 'Nakagami'); % Nakagami Probability Distribution
nak = pdf('Nakagami', i, n.mu, n.omega);
[h,p,stats] = chi2gof(datasort,'CDF',n)
plot(i, nak);

g = fitdist(datasort, 'Gamma'); % Gamma Probability Distribution
gam = pdf('Gamma', i, g.a, g.b);
[h,p,stats] = chi2gof(datasort,'CDF',g)
plot(i, gam, '.');

xlabel("input data");
ylabel("estimated pdf");
axis([0 16 0 0.3])
legend("data, (ksDensity)", "Weibull", "Nakagami (Best Fit)", "Gamma")

