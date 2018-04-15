%% Eric Wan
% recitation example set week2

clear, clc, close all

X = [ 2.6667, 3.4766, 3.7187, 4.0810, 9.4389, 5.0332, 1.0028, 4.2891, 5.2970,5.1346,4.6081 2.4199, 1.4691, 1.1571, 4.2389];
% no target was present
Y = [6.6541,7.0132,7.6296,6.0201,7.7718,1.0703,10.2948,3.9876,5.3993,14.1822, 5.0649, 5.1110 4.2848 4.6071];
% target was present

thr = 3;
Nf = X > thr;
count = sum(X > thr); % 10

Nc = Y > thr;
sum(Y > thr); %13


%{
thr 3
10/15 false alarm
1/14 miss rate
11/29 error rate
13/23 ppv 

ksdensity function


%}




