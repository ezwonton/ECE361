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

%% Finding The Optimum Plots
% Calculating ROC curve
gs0 = [zeros(70, 1) no_target];
gs1 = [ones(30, 1) target];
gs = sortrows([gs0; gs1], 2, 'descend');
counts = [0 0];
for i = 1:length(gs)
    counts(i+1, :) = [sum(gs(1:i) == 1) sum(gs(1:i) == 0)];
end
prob = [counts(:, 1)/length(target) counts(:, 2)/length(no_target)];
dist = sqrt(prob(:, 2).^2 + (1-prob(:, 1)).^2);
[~, I] = min(dist);
opt_pt = prob(I, :);
opt_PD = opt_pt(1, 1);
opt_PF = opt_pt(1, 2);
thr = gs(I, :);
opt_dist = dist(I);
figure;
hold on;
grid on;
plot(prob(:, 2), prob(:, 1), 'k');
plot(opt_PF, opt_PD, 'b*');
title("Data - Team 7")
xlabel("PF");
ylabel("PD");

roc_int = prob([(prob(:, 2) == 16/70) (prob(:, 2) == 16/70)]);
dist_int = dist(prob(:, 2) == 16/70);
plot(roc_int(2, 1), roc_int(1, 1), 'ks');
legend({"ROC: Az (area under curve)= 0.8648 | STD Deviation = 0.0451", "Opt.Pt. [PF, PD]=[0.2714, 0.80]", "Intersection: [PF, PD]=[0.2286, 0.70]"}, 'Location', 'southeast');
text(0.4, 0.5, "PPV (Optical Oper. Point) = 0.5581");
text(opt_PF, opt_PD, "  <--  Dist. to top left corner: 0.3372", 'FontWeight', 'bold', 'Color', 'r')
text(roc_int(2, 1), roc_int(1, 1), "  <--  Dist. to top left corner: 0.3772", 'FontWeight', 'bold', 'Color', 'r')

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

figure;
hold on;
grid on;
plot(n, pn, 'r');
plot(t, pt, 'k--');
title("Estimated Densities (Optimum) - Team 7")
xlabel("Input Data");
ylabel("Estimated PDF");
axis([0 9.5 0 0.45])

% Plotting optimum pdf
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
figure;
plot(n, stds, 'r');
xlabel('Number of Samples with (Target Absent)');
ylabel('Standard Deviation of ROC area');
title('Hanley and McNeill, Std Deviation with Respect to Target Absent Values');
text(45, 0.06, {"Sample Size (Target Absent) N1 = 70", "Sample Size (Target Present) N2 = 30", "Az = 0.8648"});

