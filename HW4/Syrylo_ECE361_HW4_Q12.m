%% Nicholas Syrylo Question 14
clear
clc
close
%{
12. You are given a data set consisting of 200 entries. Check if the data fits a Nakagami, gamma,
Weibull if the data set is completely positive. If data set contains negative values, test if it
normal or Laplacian. Laplacian is not a built in pdf in Matlab.
%}

% data = xlsread('HW4_data_shankar_Spring','BZ2:BZ201');
data = xlsread('HW4_data_shankar_Spring','CL2:CL201');

% a = mean(data);
% b = sqrt(var(data)/2);
% 
% [data_ks, data_points] = ksdensity(data, 'NumPoints', 200);
% 
% lap_pdf = zeros(length(data),1);
% lap_cdf = zeros(length(data),1);
% [lap_ks, lap_points] = ksdensity(lap_pdf, 'NumPoints', 200);
% 
% normal_fitdist = fitdist(data, 'normal');
% normal_pdf = pdf('normal',data_points, normal_fitdist.mu, normal_fitdist.sigma);
% [normal_ks, normal_points] = ksdensity(normal_pdf, 'NumPoints', 200);
% 
% % Data(i) = x
% 
% for i = 1:length(data)
%      lap_pdf(i) = 1/(2*b)*exp(-abs(data(i)-a)/b);
%      if data(i) > a
%          lap_cdf(i) = 1-0.5*exp(-(data(i)-a)/b);
%      elseif data(i) < a
%          lap_cdf(i) = 0.5*exp(-(data(i)-a)/b);
%      end
% end
% 
% 
% 
% lap_MSE = 0;
% normal_MSE = 0;
% 
% for i = 1:length(data)
%     lap_MSE = lap_MSE + (data_ks(i) - lap_pdf(i))^2;
%     normal_MSE = normal_MSE + (data_ks(i) - normal_pdf(i))^2;
% end
% lap_MSE = 1/length(data) * lap_MSE
% 
% normal_MSE = 1/length(data) * normal_MSE
% 
% plot(data_points, data_ks)
% hold on
% plot(data_points, normal_pdf)
% plot(data_points,lap_ks)
% 
% legend('Data', 'Normal', 'Lap')

%% Shankar's

a = mean(data);
b = sqrt(var(data)/2);

[data_ks, data_points] = ksdensity(data, 'NumPoints', 200);

normal_fitdist = fitdist(data, 'normal');
normal_pdf = pdf('normal',data_points, normal_fitdist.mu, normal_fitdist.sigma);
msen = (1/length(data)*sum((data_ks-normal_pdf).^2))

lap_pdf = (1/(2*b))*exp(-abs(data_points-a)/b);
%lap_pdf = (1/(2*b))*exp(-abs(sort(data)-a)/b);
msel  = (1/length(data)*sum((data_ks-lap_pdf).^2))

plot(data_points, data_ks)
hold on
plot(data_points, normal_pdf)
plot(data_points,lap_pdf)

legend('Data', 'Normal', 'Lap')

