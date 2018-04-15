% Eric Wan
% ezw23@drexel.edu
clear, clc

%% Homework Set 1
% Question 11
selections = 10^6;
randnum = zeros(selections, 1);
for i = 1:selections
    randnum(i) = rand(1);
end
step = 0.10:0.05:0.95;
prob = zeros(size(step, 2), 1);
c = 1;
for i = step
    H = sum(randnum < i);
    prob(c) = H/selections;
    c = c + 1;
end
figure;
hold on;
plot(step, prob, '--k');
plot(step, 1-prob, ':k');
legend("P(H)", "P(T)");
title("Question 11");
ylabel("Probabilities");


% Question 12
die = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 5, 5, 6];
%hold on;
figure;
subplot(1, 2, 1);
histogram(die);
title("Biased Die")
xlabel("Number")
ylabel("Appearances")

randnum = zeros(10^6,1);
for i = 1:10^6
    randnum(i) = randi([1,6],1);
end
subplot(1, 2, 2);
histogram(randnum);
title("Fair Die")
xlabel("Number")
ylabel("Appearances")
