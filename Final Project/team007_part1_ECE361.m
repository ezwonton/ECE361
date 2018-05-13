%% Eric Wan - ezw23@drexel.edu
clear, clc, warning('off'), close all

%{
Fs=1000;
x=linspace(0,20,Fs*20);
Ray=pdf('Rayleigh',x,1);
Ric=pdf('rician',x,3,1);
[~, index] = min(abs(Ray./Ric-1));
index
%}

%% Plotting densities of Target Absent and Target Present
data = xlsread("shankar_project_spring#7.xls", 1, "A:A");
% data = xlsread("PATH TO --> shankar_project_spring#7.xls", 1, "A:A"); %
% ^^^ may want to change this so taht it runs locally on your machine^^^

no_target = data(1:70); % first 70 - no_target present
target = data(71:100); % last 30  - target present
[pn, n] = ksdensity(no_target); % pn = probability density, n = increment over range of data
[pt, t] = ksdensity(target); % pt = probability density, t = increment over range of data
% Figure
figure;
hold on;
grid on;
plot(n, pn, 'r');
plot(t, pt, 'k--');
title("Estimated Densities (Intersect) - Team 7")
xlabel("Input Data");
ylabel("Estimated PDF");
axis([0 9.5 0 0.45])
%Finding intersection
fnX = @(x) interp1(n, pn, x) - interp1(t, pt, x); % interp over no_target - interp over target
thr = fzero(fnX, 2.5); % finds the intersection after my guess

% Plotting Prob Miss
x = 0:0.0001:thr; % x values: 0 -> thr
p = polyfit(t, pt, 20); % fitting target present ksdensity curve over x inc
p1 = polyval(p, x); % getting values of target present over x inc
area(x, p1, 'FaceColor', 'b'); % plotting area of target present up til threshold

% Plotting Prob Fail
x = thr:0.0001:n(100); % x values: thr -> end of target absent ksdensity
p = polyfit(n, pn, 20); % fitting target present ksdensity curve over x inc
p1 = polyval(p, x); % getting values of target present over x inc
area(x, p1, 'FaceColor', 'g');% plotting area of target present after threshold

% Plotting the Threshold
p = polyfit(n, pn, 20);
p1 = polyval(p, thr);
plot(thr, p1, "mo")

legend("Target Absent", "Target Present", "PM", "PF", "Thr = 2.7185"); % Legend

%% Finding Thr Optimum ROC Curve
gs0 = [zeros(70, 1) no_target];
gs1 = [ones(30, 1) target];
gs = sortrows([gs0; gs1], 2, 'descend');
counts = [0 0];
for i = 1:length(gs)
    counts(i+1, :) = [sum(gs(1:i) == 1) sum(gs(1:i) == 0)];
end
prob = [counts(:, 1)/length(target) counts(:, 2)/length(no_target)];
dist = sqrt(prob(:, 2).^2 + (1-prob(:, 1)).^2);
[M, I] = min(dist);
opt_pt = prob(I, :);
opt_PD = opt_pt(1, 1);
opt_PF = opt_pt(1, 2);
thr = gs(I, :);
figure;
hold on;
grid on;
plot(prob(:, 2), prob(:, 1), 'k');
plot(opt_PF, opt_PD, 'c*');
title("Data - Team 7")
xlabel("PF");
ylabel("PD");
legend({"ROC", "Opt.Pt. [PF, PD]=[0.15, 0.80]. Thr. Opt=3.7254"}, 'Location', 'southeast');

%% Plotting Densities of Threshold Optimum
figure;
hold on;
grid on;
plot(n, pn, 'r');
plot(t, pt, 'k--');
title("Estimated Densities (Optimum) - Team 7")
xlabel("Input Data");
ylabel("Estimated PDF");
axis([0 9.5 0 0.45])

x = 0:0.0001:thr(2);
p = polyfit(t, pt, 20);
p1 = polyval(p, x);
area(x,p1, 'FaceColor', 'b');

x = thr(2):0.0001:7;
p = polyfit(n, pn, 20);
p1 = polyval(p, x);
area(x,p1, 'FaceColor', 'g');

p = polyfit(n, pn, 20);
p1 = polyval(p, thr(2));
plot(thr(2), p1, "c*")
legend("Target Absent", "Target Present", "PM", "PF", "Thr = 2.5368");
