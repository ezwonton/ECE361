%% Eric Wan - ezw23@drexel.edu
clear, clc, warning('off'), close all

%% Problem 7
%{
Data collected from a machine vision lab to see the efficiency of an object
recognition system is given. The 200 data points given are collected as
follows: Two receivers are mounted on the vehicle receiving the
backscattered signal from the observation region. The first set of 100
points is from receiver # 1 and the second set of 100 points is from
receiver # 2.  In this experiment, the interest is to see how the
performance of the receiver could be improved. To accomplish this, three
algorithms are explored, namely the arithmetic mean, geometric mean and the
maximum. The performance is characterized by the performance index:
n = mean/std

Obtain the performance indices for the raw data, arithmetic mean, geometric
mean and the maximum.
If X and Y represent the two outputs, 
V = Arithmetic mean = (X + Y) / 2
W = Geometric mean = sqrt(XY)
Z = Maximum = max(X, Y)
%}

data = xlsread("HW7_data_shankar_Spring.xlsx", 1, "CF:CF");
[p, i] = ksdensity(data, 'NumPoints', 200);
figure;
grid on;
hold on;
plot(i, p, 'k--');

data_N = mean(data)/std(data);
arith = ([data; 0] + [0; data]) ./ 2;
[ap, ai] = ksdensity(arith, 'NumPoints', 201);
plot(ai, ap, 'r--');
geo = sqrt([data; 0].*[0; data]);
[gp, gi] = ksdensity(geo, 'NumPoints', 201);
plot(gi, gp, 'g');
max = max([[data; 0] [0; data]]);
