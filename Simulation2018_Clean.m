%This simulation of C14 is an idealized box model of the system that tries
%to recreate and parameterize the time evolution of C14. This is the clean
%copy.

clear;

%Determine which datasource to use for EUC flow rates
EUC_data = 'EUC_nino';
%EUC_data = 'EUC_soda';
%EUC_data = 'EUC_oras';

%Determine whether to use data, climatology, or constant value of mixing
Mixing_Input = 'data';
% Mixing_Input = 'Clim';
%Mixing_Input = 'Const';

%Determine whether to use data, climatology or constant value of wind
%stress for upwelling
WindStress_Input = 'data';
%WindStress_Input = 'Clim';
%WindStress_Input = 'Const';

%Whether to also plot climatology prebomb
plot_calibration = false;
% plot_calibration = true

[GalModeled, GalDICModeled, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'EUC_data', EUC_data, ...
    'Mixing_Input', Mixing_Input, 'WindStress_Input', WindStress_Input, ...
    'plot_calibration', plot_calibration);

% plot the model results
f = figure;
set(f,'Units','normalized');
set(f,'Position',[0 0 1 1]);

plot(time(1:length(GalModeled)-1), ConcentrationToDelta14(GalModeled(1:end-1),GalDICModeled(1:end-1)), 'LineWidth',2);

%determine yearly maximum and minimum values of C14
GalModeled_permil = ConcentrationToDelta14(GalModeled,GalDICModeled);
count = 1;
Max = []; Min = []; MaxTime = []; MinTime =[];
for i=25+6:12:length(time)
    lastind = i+11;
    if lastind>length(time)
        lastind=length(time);
    end
    Max(count) = max(GalModeled_permil(i:lastind));
    MaxTime(count) = time(i-1+find(Max(count)==GalModeled_permil(i:lastind),1,'first'));
    Min(count) = min(GalModeled_permil(i:lastind));
    MinTime(count) = time(i-1+find(Min(count)==GalModeled_permil(i:lastind),1,'first'));
    
    MaxObs(count) = max(GalC14(i:lastind));
    MaxTimeObs(count) = time(i-1+find(MaxObs(count)==GalC14(i:lastind),1,'first'));
    MinObs(count) = min(GalC14(i:lastind));
    MinTimeObs(count) = time(i-1+find(MinObs(count)==GalC14(i:lastind),1,'first'));
    
    count = count+1;
end

ind76 = (1976-1958);

hold on;

%calculate and plot fit lines for max and min values before and after 1976
[a] = polyfit(MaxTime(1:ind76-1),Max(1:ind76-1),1);
plot(MaxTime(1:ind76),a(2)+a(1)*MaxTime(1:ind76),'r--','LineWidth',2);
[a] = polyfit(MaxTime(ind76:end),Max(ind76:end),1);
plot(MaxTime(ind76:end),a(2)+a(1)*MaxTime(ind76:end),'r--','LineWidth',2);

[a] = polyfit(MinTime(1:ind76-1),(Min(1:ind76-1)),1);
plot(MaxTime(1:ind76),a(2)+a(1)*MaxTime(1:ind76),'r--','LineWidth',2);
[a] = polyfit(MinTime(ind76:end),(Min(ind76:end)),1);
plot(MaxTime(ind76:end),a(2)+a(1)*MaxTime(ind76:end),'r--','LineWidth',2);

xlabel('Year'); ylabel('\Delta^{14}C')
set(gca,'FontSize',16);
ylim([-90 100]);

% plot observed values
obs_color = [.1, .1, .1, .3];
plot(time,ConcentrationToDelta14(GalC14,DIC_surface),'Color',obs_color)
scatter(time,ConcentrationToDelta14(GalC14,DIC_surface), 12, 'Marker','o',...
    'MarkerEdgeColor',obs_color(1:3),'MarkerFaceColor', obs_color(1:3),...
    'MarkerFaceAlpha', obs_color(4), 'MarkerEdgeAlpha', obs_color(4));

%plot max and min best fit lines from observations
%{
[a] = polyfit(MaxTimeObs(1:ind76-1),ConcentrationToDelta14(MaxObs(1:ind76-1),2),1);
plot(MaxTimeObs(1:ind76-1),a(2)+a(1)*MaxTimeObs(1:ind76-1),'k--');
[a] = polyfit(MaxTimeObs(ind76:end),ConcentrationToDelta14(MaxObs(ind76:end),2),1);
plot(MaxTimeObs(ind76:end),a(2)+a(1)*MaxTimeObs(ind76:end),'k--');

[a] = polyfit(MinTimeObs(1:ind76-1),ConcentrationToDelta14(MinObs(1:ind76-1),2),1);
plot(MinTimeObs(1:ind76-1),a(2)+a(1)*MinTimeObs(1:ind76-1),'k--');
[a] = polyfit(MinTimeObs(ind76:end),ConcentrationToDelta14(MinObs(ind76:end),2),1);
plot(MinTimeObs(ind76:end),a(2)+a(1)*MinTimeObs(ind76:end),'k--');
%}


%%


% plot the model results
f = figure;
set(f,'Units','normalized');
set(f,'Position',[0 0 1 1]);


patch_color = [0    0.4470    0.7410  .4];
%--------------------Deep DIC---------------------------
subplot(5, 1, 1);
hold on;
[GalModeled_lower, GalDICModeled_lower, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'DIC_deep', 1.7);

[GalModeled_upper, GalDICModeled_upper, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'DIC_deep', 2.7);

patch([time(1:length(GalModeled_upper)-1) fliplr(time(1:length(GalModeled_lower)-1))],...
    [ConcentrationToDelta14(GalModeled_upper(1:end-1),GalDICModeled_upper(1:end-1)) ...
    fliplr(ConcentrationToDelta14(GalModeled_lower(1:end-1),GalDICModeled_lower(1:end-1)))], ...
    patch_color(1:3))
alpha(0.4) 

[GalModeled, GalDICModeled, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0);
plot(time(1:length(GalModeled)-1), ConcentrationToDelta14(...
    GalModeled(1:end-1),GalDICModeled(1:end-1)), 'LineWidth',2, 'Color', ...
    patch_color(1:3));

ylabel('\Delta^{14}C');
title('Sensitivity to Deep Ocean DIC');
set(gca,'FontSize',16);
ylim([-90 100]);
set(gca,'XTick',[])

% plot observed values
obs_color = [.1, .1, .1, .3];
plot(time,ConcentrationToDelta14(GalC14,DIC_surface),'Color',obs_color)
scatter(time,ConcentrationToDelta14(GalC14,DIC_surface), 12, 'Marker','o',...
    'MarkerEdgeColor',obs_color(1:3),'MarkerFaceColor', obs_color(1:3),...
    'MarkerFaceAlpha', obs_color(4), 'MarkerEdgeAlpha', obs_color(4));

%--------------------EUC DIC---------------------------
subplot(5, 1, 2);
hold on;
[GalModeled_lower, GalDICModeled_lower, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'DIC_vt', 2);

[GalModeled_upper, GalDICModeled_upper, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'DIC_vt', 2.1);

patch([time(1:length(GalModeled_upper)-1) fliplr(time(1:length(GalModeled_lower)-1))],...
    [ConcentrationToDelta14(GalModeled_upper(1:end-1),GalDICModeled_upper(1:end-1)) ...
    fliplr(ConcentrationToDelta14(GalModeled_lower(1:end-1),GalDICModeled_lower(1:end-1)))], ...
    patch_color(1:3))
alpha(0.4) 

[GalModeled, GalDICModeled, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0);
plot(time(1:length(GalModeled)-1), ConcentrationToDelta14(...
    GalModeled(1:end-1),GalDICModeled(1:end-1)), 'LineWidth',2, 'Color', ...
    patch_color(1:3));

ylabel('\Delta^{14}C');
title('Sensitivity to EUC DIC');
set(gca,'FontSize',16);
ylim([-90 100]);
set(gca,'XTick',[])

% plot observed values
obs_color = [.1, .1, .1, .3];
plot(time,ConcentrationToDelta14(GalC14,DIC_surface),'Color',obs_color)
scatter(time,ConcentrationToDelta14(GalC14,DIC_surface), 12, 'Marker','o',...
    'MarkerEdgeColor',obs_color(1:3),'MarkerFaceColor', obs_color(1:3),...
    'MarkerFaceAlpha', obs_color(4), 'MarkerEdgeAlpha', obs_color(4));

%--------------------dC14 Prebomb---------------------------
subplot(5, 1, 3);
hold on;
[GalModeled_lower, GalDICModeled_lower, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'dC14_prebomb', -80);

[GalModeled_upper, GalDICModeled_upper, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'dC14_prebomb', -110);

patch([time(1:length(GalModeled_upper)-1) fliplr(time(1:length(GalModeled_lower)-1))],...
    [ConcentrationToDelta14(GalModeled_upper(1:end-1),GalDICModeled_upper(1:end-1)) ...
    fliplr(ConcentrationToDelta14(GalModeled_lower(1:end-1),GalDICModeled_lower(1:end-1)))], ...
    patch_color(1:3))
alpha(0.4) 

[GalModeled, GalDICModeled, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0);
plot(time(1:length(GalModeled)-1), ConcentrationToDelta14(...
    GalModeled(1:end-1),GalDICModeled(1:end-1)), 'LineWidth',2, 'Color', ...
    patch_color(1:3));

ylabel('\Delta^{14}C');
title('Sensitivity to Prebomb \Delta^{14}C levels');
set(gca,'FontSize',16);
ylim([-90 100]);
set(gca,'XTick',[])

% plot observed values
obs_color = [.1, .1, .1, .3];
plot(time,ConcentrationToDelta14(GalC14,DIC_surface),'Color',obs_color)
scatter(time,ConcentrationToDelta14(GalC14,DIC_surface), 12, 'Marker','o',...
    'MarkerEdgeColor',obs_color(1:3),'MarkerFaceColor', obs_color(1:3),...
    'MarkerFaceAlpha', obs_color(4), 'MarkerEdgeAlpha', obs_color(4));

%--------------------dC14 Prebomb---------------------------
subplot(5, 1, 4);
hold on;
[GalModeled_lower, GalDICModeled_lower, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'ml_depth', 15);

[GalModeled_upper, GalDICModeled_upper, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'ml_depth', 50);

patch([time(1:length(GalModeled_upper)-1) fliplr(time(1:length(GalModeled_lower)-1))],...
    [ConcentrationToDelta14(GalModeled_upper(1:end-1),GalDICModeled_upper(1:end-1)) ...
    fliplr(ConcentrationToDelta14(GalModeled_lower(1:end-1),GalDICModeled_lower(1:end-1)))], ...
    patch_color(1:3))
alpha(0.4) 

[GalModeled, GalDICModeled, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0);
plot(time(1:length(GalModeled)-1), ConcentrationToDelta14(...
    GalModeled(1:end-1),GalDICModeled(1:end-1)), 'LineWidth',2, 'Color', ...
    patch_color(1:3));

ylabel('\Delta^{14}C');
title('Sensitivity to Mixed Layer Depth');
set(gca,'FontSize',16);
ylim([-90 100]);
set(gca,'XTick',[])

% plot observed values
obs_color = [.1, .1, .1, .3];
plot(time,ConcentrationToDelta14(GalC14,DIC_surface),'Color',obs_color)
scatter(time,ConcentrationToDelta14(GalC14,DIC_surface), 12, 'Marker','o',...
    'MarkerEdgeColor',obs_color(1:3),'MarkerFaceColor', obs_color(1:3),...
    'MarkerFaceAlpha', obs_color(4), 'MarkerEdgeAlpha', obs_color(4));

%--------------------Bioactivity---------------------------
subplot(5, 1, 5);
hold on;
[GalModeled_lower, GalDICModeled_lower, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'biouptakerate', 8*10^-3);

[GalModeled_upper, GalDICModeled_upper, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0, 'biouptakerate', 8*10^-5);

patch([time(1:length(GalModeled_upper)-1) fliplr(time(1:length(GalModeled_lower)-1))],...
    [ConcentrationToDelta14(GalModeled_upper(1:end-1),GalDICModeled_upper(1:end-1)) ...
    fliplr(ConcentrationToDelta14(GalModeled_lower(1:end-1),GalDICModeled_lower(1:end-1)))], ...
    patch_color(1:3))
alpha(0.4) 

[GalModeled, GalDICModeled, GalC14, DIC_surface, time] = ...
    radiocarbon_model_simulation(0);
plot(time(1:length(GalModeled)-1), ConcentrationToDelta14(...
    GalModeled(1:end-1),GalDICModeled(1:end-1)), 'LineWidth',2, 'Color', ...
    patch_color(1:3));

ylabel('\Delta^{14}C');
title('Sensitivity to Surface Biological Update');
set(gca,'FontSize',16);
xlabel('Year');
ylim([-90 100]);

% plot observed values
obs_color = [.1, .1, .1, .3];
plot(time,ConcentrationToDelta14(GalC14,DIC_surface),'Color',obs_color)
scatter(time,ConcentrationToDelta14(GalC14,DIC_surface), 12, 'Marker','o',...
    'MarkerEdgeColor',obs_color(1:3),'MarkerFaceColor', obs_color(1:3),...
    'MarkerFaceAlpha', obs_color(4), 'MarkerEdgeAlpha', obs_color(4));
%%

%----------------------INVERSE SIMULATION-----------------------------------%
syms x
%monthly timestep through the model
for i=1:length(GalC14)-1
    %parameter calculation for DIC solubility - air-gas exchange
    [Ko,k] = DIC_solubility(Nino3(i),34.78,Wind_eq(i));
    
    %invasion of CO2 and C14 from atmosphere
    InvasionCO2 = k*Ko*AtmCO2(i)/1e6; %[mol/m^2/month]
    InvasionC14 = InvasionCO2*(Delta14Todelta(AtmC14(i),AtmC13(i))/1000+1)*Rstd;
    
    %biological update of DIC and C14
    BioUptake = 20.09*10^-3*365/12/ml_depth; %From Chai et al Fig. 7
    BioUptakeC14 = BioUptake*GalModeled(i)/GalDICModeled(i);
    
    %outgasing of CO2 and C14 from ocean
    EvasionCO2 = k*Ko*carbonate_calculation(Nino3(i),GalDICModeled(i))/1e6;
    EvasionC14 = EvasionCO2*GalModeled(i)/GalDICModeled(i);
    
    eqn = GalC14(i+1) == (vReplace(i)+x)*(x*DeepC14+x*GalC14(i)+(1-2*x)*EUCc14(i))+(1-vReplace(i)-x)*GalC14(i)+(InvasionC14-EvasionC14)/ml_depth-BioUptakeC14;
    solutions = vpasolve(eqn, x);
    [val, idx]=min(abs(solutions-1));
    vel(i) = real(max(solutions(idx),0));
end
veln = sqrt(vel/.2)*mean(vEUC)/mean(sqrt(vel/.2));

f = figure;
set(f,'Units','normalized');
set(f,'Position',[0 0 1 1]);

plot(time(1:length(veln)), veln, 'LineWidth',2);
xlabel('Year'); ylabel('Maximum EUC Velocity (cm/s)')
set(gca,'FontSize',16);
title('Inferred EUC Peak Velocity')