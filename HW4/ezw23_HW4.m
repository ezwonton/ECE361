%% Eric Wan - ezw23@drexel.edu
clear, clc, warning('on'), close all

%% Question 12
%{
You are given a data set consisting of 200 entries. Check if the data fits
a Nakagami, gamma, Weibull if the data set is completely positive. If data
set contains negative values, test if it normal or Laplacian. Laplacian is
not a built in pdf in Matlab.

To get the parameters of the Laplacian,
simply find the mean and variance as indicated above. Once you have the
mean, you may substitute the values of a and b in f(x) to get the pdf. For
the other densities, the Matlab command is  pd = fitdist(data,distname)

Using the estimated parameters, you may get the theoretical density using
the Matlab command pd = pdf(distname,x,parameters)

Using the estimated parameters, plot the estimated densities. Plot the data
pdf using ksdensity and in each case, estimate the MSE. Determine the best
fit based on the lowest value of MSE.
%}

data = xlsread("HW4_data_shankar_Spring.xlsx", 1, "CF:CF");
datasort = sort(data);
syms x;
z = data(data < 0);
disp(z);
[p, i] = ksdensity(data, 'NumPoints', 200); % pdf = probability density, i = increment over range of data
w = fitdist(data, 'Weibull'); % Weibull Probability Distribution
wei = pdf('Weibull', i, w.A, w.B);
%wei = (w.b/w.B)*(x^(w.B-1)/w.A)*exp(-(x^w.B)/w.A);
msew = (1/length(data)*sum((p-wei).^2));

n = fitdist(data, 'Nakagami'); % Nakagami Probability Distribution
nak = pdf('Nakagami', i, n.mu, n.omega);
%nak = 2*(n.mu/n.omega)^(n.mu)*(x^(2*n.mu-1)/gamma(n.mu))*exp(-n.mu*x^2/n.omega);
msen = (1/length(data)*sum((p-nak).^2));

g = fitdist(data, 'Gamma'); % Gamma Probability Distribution
gam = pdf('Gamma', i, g.a, g.b);
%gam = (1/(g.b^(g.a)*gamma(g.a)))*(x^(-1)/g.a)*exp(-x/g.b);
mseg = (1/length(data)*sum((p-gam).^2));

fprintf("Weibull MSE: %f; Parameters: A = %f, B = %f\n", msew, w.A, w.B);
fprintf("Nakagami MSE: %f; Parameters: mu = %f, omega = %f\n", msen, n.mu, n.omega);
fprintf("Gamma MSE: %f; Parameters: a = %f, b = %f\n", mseg, g.a, g.b);

figure;
grid on;
hold on;
plot(i, p);
plot(i, wei);
plot(i, nak);
plot(i, gam);
axis([0 12 0 0.3])
legend("data (ksDensity)", "Weibull MSE: 0.000458", "Nakagami MSE: 0.000456", "Gamma MSE: 0.000527")


