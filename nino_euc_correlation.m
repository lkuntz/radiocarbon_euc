% Loads in temperature and uvelocity data from a control run of the GCM to 
% determine the correlation coefficent and linear fit.

load('EUCcont.mat', 'temp'); % nino 3.4 temp
load('EUCcont.mat', 'uvel'); % lon x lat (only on eq) x depth x time
load('EUCcont.mat', 'ulat'); 
load('EUCcont.mat', 'ulon');

[m, ind_220] = min(abs(ulon-220));
uvel220 = uvel(ind_220,:,:,:);
uvel220 = squeeze(uvel220);
uvel220(uvel220<0)= 0;
max_uvel220 = max(uvel220);

scatter(max(uvel220),temp);
corrcoef(squeeze(max(uvel220))', squeeze(temp))

p = polyfit(squeeze(max(uvel220))', squeeze(temp), 1);
v1 = linspace(20, 160);
t1 = polyval(p, v1);
hold on;
plot(v1, t1, 'Color', 'black', 'LineWidth', 2);

xlabel('Maximum Velocity (cm/s)', 'FontSize', 16);
ylabel('Temperature in Nino3.4 (C)', 'FontSize', 16);
