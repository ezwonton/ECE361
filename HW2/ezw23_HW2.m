%% Eric Wan
% ezw23@drexel.edu
clc, clear, close all

%% Question 12
X=[1.6978 7.7368 8.7591 1.9001 0.7486 0.2048 6.1228, 1.6926 2.2690 13.2931];
Y=[11.4976,4.5624,17.0159,5.1402,8.9841,8.9526,2.0939,6.1407];
X = sort(X);
Y = sort(Y);

for i = [3, 5, 6]
    thr = i; % threshold
    Nx = numel(X); % # of elements in X (No target)
    Ny = numel(Y); % # of elements in Y (W/ target)
    N = Nx + Ny; % # of total elements

    Nf = X > thr; % samples X > threshold
    count = sum(X > thr); % # of samples X > threshold
    pTa = Nx/N; % a priori probability for X (no target was present)
    pF = count/Nx; % probability of false alarm (X > thr / count X)

    Nc = Y > thr; % samples Y- > threshold
    count = sum(Y > thr); % # of samples Y > threshold
    pTp = Ny/N; % a priori probabilty for Y (target was present)
    pD = count/Ny; % probability of detection (Y > thr / count Y)
    pM = 1 - pD;

    pE = pM * pTp + pF * pTa; % prob error = prob miss * a priori prob present + prob fail * apriori prob not present
    % a priori prob, false alarm, miss, error rate
    % thresh of 3, 5, 6
    sprintf("a priori probabilty for X (target was not present) of threshold (%f): %f", thr, pTa)
    sprintf("a priori probabilty for Y (target was present) of threshold (%f): %f", thr, pTp)
    sprintf("probabilty for a false alarm (target was not present but detects as such) of threshold (%f): %f", thr, pF)
    sprintf("probabilty for a miss (target was present but not detects) of threshold (%f): %f", thr, pM)
    sprintf("probabilty of error of threshold (%f): %f", thr, pE)
end
[x, xt] = ksdensity(X);
[y, yt] = ksdensity(Y);
figure;
hold on;
plot(xt, x);
plot(yt, y);



%% Question 13
data = xlsread("HW2_data_shankar_Spring.xlsx", "hw5data_", "CE:CE");
data1 = data(1:40,:); % target is present
data1 = sort(data1);
data2 = data(41:70,:); % target is not present
data2 = sort(data2);
[p, pt] = ksdensity(data1);
[np, npt]= ksdensity(data2);
%{
figure;
hold on;
hist(p);
hist(np);
%}
figure;
hold on;
plot(pt, p);
plot(npt, np);
area(pt(pt>4.1 & pt<pt(end)), p(pt>4.1 & pt<pt(end)));



thr = 4.7; % threshold
Nx = numel(data2); % # of elements in X (No target)
Ny = numel(data1); % # of elements in Y (W/ target)
N = Nx + Ny; % # of total elements

Nf = X > thr; % samples X > threshold
count = sum(X > thr); % # of samples X > threshold
pTa = Nx/N; % a priori probability for X (no target was present)
pF = count/Nx; % probability of false alarm (X > thr / count X)

Nc = Y > thr; % samples Y- > threshold
count = sum(Y > thr); % # of samples Y > threshold
pTp = Ny/N; % a priori probabilty for Y (target was present)
pD = count/Ny; % probability of detection (Y > thr / count Y)
pM = 1 - pD;

pE = pM * pTp + pF * pTa; % prob error = prob miss * a priori prob present + prob fail * apriori prob not present
% a priori prob, false alarm, miss, sensitivity, specificity, positive predictive value, accuracy
sprintf("a priori probabilty for X (target was not present) of threshold (%f): %f", thr, pTa)
sprintf("a priori probabilty for Y (target was present) of threshold (%f): %f", thr, pTp)
sprintf("probabilty for a false alarm (target was not present but detects as such) of threshold (%f): %f", thr, pF)
sprintf("probabilty for a miss (target was present but not detects) of threshold (%f): %f", thr, pM)


% ^^sensitivity, specificity, positive predictive value, accuracy
