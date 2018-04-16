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
title("Question 12")
xlabel("Values");
ylabel("Counts");
legend("Target Absent", "Target Present");

%% Question 13
data = xlsread("HW2_data_shankar_Spring.xlsx", "hw5data_", "CE:CE");
data1 = data(1:40); % target is not present
data1 = sort(data1);
data2 = data(41:70); % target is present
data2 = sort(data2);
thr = (median(data1) + median(data2))/2;

figure;
hold on;
grid on;
histogram(data1, 6);
histogram(data2, 6);
line(thr, 10, 'Color', 'red', 'LineWidth', 3);
title("Question 13: WAN");
xlabel("Values");
ylabel("Counts");
legend("Target Absent", "Target Present");

Nx = numel(data1); % # of elements in d2 (No target)
Ny = numel(data2); % # of elements in d1 (W/ target)
N = Nx + Ny; % # of total elements

Nf = sum(data1 > thr); % samples d2 > threshold
pTa = Nx/N; % a priori probability for d2 (no target was present)
pF = Nf/Nx; % probability of false alarm (d2 > thr / count d2)

Nc = sum(data2 > thr); % samples d1 > threshold
pTp = Ny/N; % a priori probabilty for d1 (target was present)
pD = Nc/Ny; % probability of detection (d1 > thr / count d1)
pM = 1 - pD; % probabilty of miss

pSens = Nf/Ny; % # of correct positive detections / # of positive cases
pSpec = Nc/Nx; % # of correct negative detections / # of negative cases
PPV = Nc / (Nc + Nf); % # of correct positive detections / # of positive detections
pACC = (Nc + (Nx - Nf)) / N; % # of correct detections / # of cases

pE = pM * pTp + pF * pTa; % prob error = prob miss * a priori prob present + prob fail * apriori prob not present
% a priori prob, false alarm, miss, sensitivity, specificity, positive predictive value, accuracy
sprintf("a priori probabilty for X of threshold (%f): %f", thr, pTa)
sprintf("a priori probabilty for Y of threshold (%f): %f", thr, pTp)
sprintf("probabilty for a false alarm of threshold (%f): %f", thr, pF)
sprintf("probabilty for a miss of threshold (%f): %f", thr, pM)
sprintf("sensitivity of threshold (%f): %f", thr, pSens)
sprintf("specificity of threshold (%f): %f", thr, pSpec)
sprintf("positive predictive value of threshold (%f): %f", thr, PPV)
sprintf("accuracy of threshold (%f): %f", thr, pACC)

fprintf('------Confusion Matrix for Histogram w/ Threshold %f------\n', thr);
fprintf('%-15s\t %-15s\t %-15s\t %-15s\n', 'Data', 'Target', 'Target', 'Total')
fprintf('%-15s\t %-15s\t %-15s\t %-15s\n', 'Collected', 'Detected', 'Not Detected', 'Count')
fprintf('%-15s\t %-15d\t %-15d\t %-15d\n', 'Target Absent',Nf, Nx - Nf, Nx)
fprintf('%-15s\t %-15d\t %-15d\t %-15d\n', 'Target Present', Nc, Ny - Nc, Ny)
fprintf('%-15s\t %-15d\t %-15d\t %-15d\n', 'Total Count', Nf + Nc, Ny - Nf + Ny - Nc, N)

pt = 0:0.0001:10;
npt = 0:0.0001:25;
p = ksdensity(data1, pt);
np = ksdensity(data2, npt);
thr = sum((pt((abs(p - np(1:numel(p))) < 0.00001)) + npt((abs(p - np(1:numel(p))) < 0.00001)))/2)/2;
figure;
hold on;
plot(pt, p); % absent
plot(npt, np); % present
area(pt(pt>thr & pt<pt(end)), p(pt>thr & pt<pt(end))); % miss
area(npt(npt<thr & npt>npt(1)), np(npt<thr & npt>npt(1))); % false alarm
title("Question 13: WAN");
xlabel("Values");
ylabel("Est. Density");
legend("Target Absent", "Target Present", "False Alarm", "Miss");

Nx = numel(data1); % # of elements in d2 (No target)
Ny = numel(data2); % # of elements in d1 (W/ target)
N = Nx + Ny; % # of total elements

Nf = sum(data1 > thr); % samples d2 > threshold
pTa = Nx/N; % a priori probability for d2 (no target was present)
pF = Nf/Nx; % probability of false alarm (d2 > thr / count d2)

Nc = sum(data2 > thr); % samples d1 > threshold
pTp = Ny/N; % a priori probabilty for d1 (target was present)
pD = Nc/Ny; % probability of detection (d1 > thr / count d1)
pM = 1 - pD; % probabilty of miss

pSens = Nf/Ny; % # of correct positive detections / # of positive cases
pSpec = Nc/Nx; % # of correct negative detections / # of negative cases
PPV = Nc / (Nc + Nf); % # of correct positive detections / # of positive detections
pACC = (Nc + (Nx - Nf)) / N; % # of correct detections / # of cases

pE = pM * pTp + pF * pTa; % prob error = prob miss * a priori prob present + prob fail * apriori prob not present
% a priori prob, false alarm, miss, sensitivity, specificity, positive predictive value, accuracy
sprintf("a priori probabilty for X of threshold (%f): %f", thr, pTa)
sprintf("a priori probabilty for Y of threshold (%f): %f", thr, pTp)
sprintf("probabilty for a false alarm of threshold (%f): %f", thr, pF)
sprintf("probabilty for a miss of threshold (%f): %f", thr, pM)
sprintf("sensitivity of threshold (%f): %f", thr, pSens)
sprintf("specificity of threshold (%f): %f", thr, pSpec)
sprintf("positive predictive value of threshold (%f): %f", thr, PPV)
sprintf("accuracy of threshold (%f): %f", thr, pACC)

fprintf('------Confusion Matrix for kDensity w/ Threshold %f------\n', thr);
fprintf('%-15s\t %-15s\t %-15s\t %-15s\n', 'Data', 'Target', 'Target', 'Total')
fprintf('%-15s\t %-15s\t %-15s\t %-15s\n', 'Collected', 'Detected', 'Not Detected', 'Count')
fprintf('%-15s\t %-15d\t %-15d\t %-15d\n', 'Target Absent',Nf, Nx - Nf, Nx)
fprintf('%-15s\t %-15d\t %-15d\t %-15d\n', 'Target Present', Nc, Ny - Nc, Ny)
fprintf('%-15s\t %-15d\t %-15d\t %-15d\n', 'Total Count', Nf + Nc, Ny - Nf + Ny - Nc, N)

