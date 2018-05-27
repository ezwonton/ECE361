%% Eric Wan - ezw23@drexel.edu
clear, clc, warning('off'), close all

%% Plotting densities of Target Absent and Target Present
data = xlsread("shankar_project_spring#7.xls", 1, "A:A");
% data = xlsread("PATH TO --> shankar_project_spring#7.xls", 1, "A:A"); %

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

legend("Target Absent", "Target Present"); % Legend

%% Chi Squared Testing

% TARGET PRESENT DENSITIES:
% rayleigh
ray_pd = fitdist(target, 'rayleigh');
[ray_target_h, ray_p, ray_target_stats] = chi2gof(target, 'CDF', ray_pd, 'emin',3);
%dataN = random('gamma', pd.a, pd.b, 10000000, 1);
%[gamma_p, gamma_t] = ksdensity(dataN); % pt = probability density, t = increment over range of data
%plot(gamma_t, gamma_p);
% gamma
gam_pd = fitdist(target, 'gamma');
[gamma_target_h, gam_p, gamma_target_stats] = chi2gof(target, 'CDF', gam_pd, 'emin', 3);
%dataN = random('gamma', pd.a, pd.b, 10000000, 1);
%[gamma_p, gamma_t] = ksdensity(dataN); % pt = probability density, t = increment over range of data
%plot(gamma_t, gamma_p);
% nakagami
nak_pd = fitdist(target, 'nakagami');
[nak_target_h, nak_p, nak_target_stats] = chi2gof(target, 'CDF', nak_pd, 'emin', 4);
%dataN = random('gamma', pd.a, pd.b, 10000000, 1);
%[gamma_p, gamma_t] = ksdensity(dataN); % pt = probability density, t = increment over range of data
%plot(gamma_t, gamma_p);
% rician
ric_pd = fitdist(target, 'rician');
[ric_target_h, ric_p, ric_target_stats] = chi2gof(target, 'CDF', ric_pd, 'emin', 3);
%dataN = random('gamma', pd.a, pd.b, 10000000, 1);
%[gamma_p, gamma_t] = ksdensity(dataN); % pt = probability density, t = increment over range of data
%plot(gamma_t, gamma_p);


% TARGET NO PRESENT DENSITIES:
% rayleigh
ray_pd = fitdist(no_target, 'rayleigh');
[ray_target_h, ray_p, ray_target_stats] = chi2gof(no_target, 'CDF', ray_pd);
% gamma
gam_pd = fitdist(no_target, 'gamma');
[gam_target_h, gam_p, gam_target_stats] = chi2gof(no_target, 'CDF', gam_pd);
% nakagami
nak_pd = fitdist(no_target, 'nakagami');
[nak_target_h, nak_p, nak_target_stats] = chi2gof(no_target, 'CDF', nak_pd);
% rician
ric_pd = fitdist(no_target, 'rician');
[ric_target_h, ric_p, ric_target_stats] = chi2gof(no_target, 'CDF', ric_pd);


%% Plotting Theoretical Densities 
figure;
hold on;
grid on;
% Target Found Best Fit Theoretical (NAKAGAMI)
dataN = random('nakagami', 2.00199, 15.0637, 10000000, 1);
[target_theor_p, target_theor_t] = ksdensity(dataN); % pt = probability density, t = increment over range of data
plot(target_theor_t, target_theor_p, 'LineWidth', 2);
% Data Target
[target_p, target_t] = ksdensity(target); % pt = probability density, t = increment over range of data
plot(target_t, target_p, '--', 'LineWidth', 2);
% No Target
dataN = random('gamma', 5.36663, 0.391895, 10000000, 1);
[notarget_theor_p, notarget_theor_t] = ksdensity(dataN); % pt = probability density, t = increment over range of data
plot(notarget_theor_t, notarget_theor_p, 'LineWidth', 2);
% Data No Target
[notarget_p, notarget_t] = ksdensity(no_target); % pt = probability density, t = increment over range of data
plot(notarget_t, notarget_p, '--', 'LineWidth', 2);

% Plotting Labels
title("Theoretical Fit Densities - Team 7")
xlabel("values");
ylabel("pdf");
legend('theoretical fit (Target Present):nakagami(2.00199, 15.0637)', 'data(Target Present)', 'theoretical fit (Target Absent):gamma(5.36663, 0.391895)', 'data(Target Absent)');

%% Finding The Optimum
figure;
hold on;
grid on;
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
plot(prob(:, 2), prob(:, 1), 'k');
plot(opt_PF, opt_PD, 'r*');
% ROC PLOT from estimated
theor_target = random('nakagami', 2.00199, 15.0637, 100000, 1);
theor_no_target = random('gamma', 5.36663, 0.391895, 100000, 1);

gs0 = [zeros(100000, 1) theor_no_target];
gs1 = [ones(100000, 1) theor_target];
disp('Created Stuff ... ');
gs = sortrows([gs0; gs1], 2, 'descend');
counts = [0 0];
for i = 1:length(gs)
    disp(i);
    counts(i+1, :) = [sum(gs(1:i) == 1) sum(gs(1:i) == 0)];
end
disp('Done Iterating');
prob = [counts(:, 1)/length(theor_target) counts(:, 2)/length(theor_no_target)];
dist = sqrt(prob(:, 2).^2 + (1-prob(:, 1)).^2);
[M, I] = min(dist);
opt_pt = prob(I, :);
opt_PD = opt_pt(1, 1);
opt_PF = opt_pt(1, 2);
thr = gs(I, :);
plot(prob(:, 2), prob(:, 1), 'k');
plot(opt_PF, opt_PD, '--');

title("ROC Plot - Team 7");
xlabel("Prob. False Alaram");
ylabel("Prob. Detection");
legend("input data Az = 0.8648", "theoretical fit Az = 0.8345");
% Finding Area of ROC for theor.
total_area = 0;
last = -1;
for i = 1:length(prob(:,1))
    if prob(i, 1) > last
        last = prob(i, 1);
        addit = 1.0e-5 * (1-prob(i, 2));
        total_area = total_area + addit;
    end
end

figure;
hold on;
grid on;
plot(n, pn, 'r');
plot(t, pt, 'k--');
title("Estimated Densities (Optimum) - Team 7")
xlabel("Input Data");
ylabel("Estimated PDF");
axis([0 9.5 0 0.45])

[M , I] = min(dist);
opt_pt = prob(I, :);
opt_PD = opt_pt(1, 1);
opt_PF = opt_pt(1, 2);

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
plot(thr(2), p1, "mo")
legend("Target Absent", "Target Present", "PM", "PF", "Thr = 2.5368");
text(5, 0.25, "Perfomance Index = 0.9221");

%% Creating Std Devation vs. Number of Samples (Target Not Present) Curve
az = [];
stds =  [];
for n=0:5:40
    prob = [counts(:, 1)/length(target) counts(:, 2)/(30+n)];
    % Finding Area of ROC
    total_area = 0;
    last = -1;
    for i = 1:length(prob(:,1))
        if prob(i, 1) > last
            last = prob(i, 1);
            addit = 0.0333 * (1-prob(i, 2));
            total_area = total_area + addit;
        end
    end
    % Store all areas
    az = [az total_area];
    % Store all std deviations of areas
    A1 = (total_area)/(2-total_area);
    A2 = (2*total_area^2)/(1+total_area);
    std_i = (total_area*(1-total_area)+(30-1)*(A1-total_area^2) + (30+n-1)*(A2-total_area^2)) / (30*(30+n));
    stds = [stds sqrt(std_i)];
end
n=30:5:70;
plot(n, stds);
xlabel('Number of Samples with (Target Absent)');
ylabel('Standard Deviation of ROC area');
title('Hanley and McNeill, Std Deviation with Respect to Target Absent Values');
text(45, 0.06, {"Sample Size (Target Absent)N1 = 70", "Sample Size (Target Present)N2 = 30", "Az = 0.8648"});

