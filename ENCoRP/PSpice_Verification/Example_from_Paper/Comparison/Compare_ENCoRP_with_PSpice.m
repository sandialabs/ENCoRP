%
% Housekeeping
%
clear
close all
clc

%
% Plot all VI comparisons (0) OR cycle through one at a time (1)
%
all_or_1by1 = 0;             % 0 = all, 1 (nonzero) = 1 by 1
freqdomain = [1,1e9];        % Frequency domain

%
% Load in SPICE data
%
SD = csvread('PSpice_Data.csv',2);  
fs = SD(:,1);
VIs = SD(:,2:end);

%
% Load in data from ENCoRP
%
load('ENCoRP_Data.mat'); 
f = Data.Freq;
VIo = abs(Data.NodeVI); 
LV = abs(Data.ListVI)+1;
LV2 = LV(1,:);

%
% Subplot size for plotting if all
%
if ~all_or_1by1
    if ceil(sqrt(size(LV,2)))*(ceil(sqrt(size(LV,2)))-1) > size(LV,2)
        sz = [ceil(sqrt(size(LV,2))),(ceil(sqrt(size(LV,2)))-1)];
    else
        sz = [ceil(sqrt(size(LV,2))),ceil(sqrt(size(LV,2)))];
    end
end

%
% Plot
%
figure
for i = 1:size(LV,2)

    % If all, create subplots appropriately
    if ~all_or_1by1
        subplot(sz(2),sz(1),i)
    end

    % Plot PSpice data vs ENCoRP data on log-frequency scale
    semilogx(f,VIo(:,i),'LineWidth',2)
    hold on
    plot(fs,VIs(:,i),'--','LineWidth',2)
    if exist('freqdomain','var') 
        xlim(freqdomain)
    end

    % Assign ylabel
    if find(LV2(i)==LV2,1)==i
        ylabel('Voltage (V)')
    else
        ylabel('Current (A)')
    end

    % Create x labels
    if all_or_1by1 || (i > size(LV,2)-sz(1))
        xlabel('Frequency (Hz)')
    end

    % If one by one, clean up figure
    if all_or_1by1
        legend('ENCoRP','PSpice')
        set(gca,'FontSize',14)
        box on
        drawnow;
        pause
        clf
    end
end
if ~all_or_1by1
    legend('ENCoRP','PSpice')
else
    close all
end
