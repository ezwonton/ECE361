%% Eric Wan - ezw23@drexel.edu
clear, clc, warning('off'), close all

%% Question 10
%{
Data collected from a machine vision lab to see the efficieny of an object
recognition system is given (HW5_data_shankar). The IR data set consist of
a column of samples (one for each of you as indicated by your last name),
the first 40 values are the responses when there is no target in the field
of view of the robot with remaining 30 values are the responses when there
is target in field of view. As shown in the lecture, obtain the receiver
operating characteristic curve for the system, determine the area under the
ROC curve, optimal operational point, the positive predictive value
corresponding to the
%}
data = xlsread("HW5_data_shankar_Spring.xlsx", 1, "CF:CF");
no_target = data(1:40);
target = data(41:70);

gs0 = [zeros(40, 1) data(1:40)];
gs1 = [ones(30, 1) data(41:70)];
gs = sortrows([gs0 ; gs1], 2, 'descend');
counts = [0 0];
for i = 1:length(gs)
%    counts = [counts; sum(gs(1:i) == 1) sum(gs(1:i) == 0)];
    counts(i+1,:) = [sum(gs(1:i) == 1) sum(gs(1:i) == 0)];
end
prob = [counts(:,1)/length(target) counts(:,2)/length(no_target)];
dist = sqrt(prob(:,2).^2 + (1-prob(:,1)).^2);
[M , I] = min(dist);
opt_pt = prob(I,:);
opt_PD = opt_pt(1,1);
opt_PF = opt_pt(1,2);
thr = gs(I,:);
figure;
hold on;
grid on;
plot(prob(:,2),prob(:,1));
plot(opt_PF, opt_PD, 'r*');
title("Data - Wan")
xlabel("PF");
ylabel("PD");
legend({"ROC", "Opt.Pt. [PF, PD]=[0.15, 0.80]. Thr.=3.7254"}, 'Location', 'southeast');

[pn, n] = ksdensity(no_target); % pn = probability density, n = increment over range of data
[pt, t] = ksdensity(target); % pt = probability density, t = increment over range of data
figure;
hold on;
grid on;
plot(n, pn);
plot(t, pt);
title("Estimated Densities")
xlabel("input data");
ylabel("estimated pdf");
axis([0 20 0 0.4])

x = 0:0.001:thr(2);
p = polyfit(t,pt,20);
p1 = polyval(p,x);
area(x,p1);

x = thr(2):0.001:6.0;
p = polyfit(n,pn,20);
p1 = polyval(p,x);
area(x,p1);

p = polyfit(n,pn,20);
p1 = polyval(p,thr(2));
plot(thr(2), p1, "g*")
legend("target absent", "target present", "PM", "PF", "Thr=3.7254");
