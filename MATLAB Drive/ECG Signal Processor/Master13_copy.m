%% 1:Load Data
%General Idea is just loading all the data into two arrays and then
%clearing it for a new file.
% t- prefix denotes time instance
% ecg- prefix denotes ecg reading data
load("BoxBreathing.mat")
t_box = b1(:,1);
ecg_box = b1(:,2);

clear b1

load("Exercise.mat")
t_ex = b1(:,1);
ecg_ex = b1(:,2);

clear b1

load("Resting.mat")

t_rest = b1(:,1);
ecg_rest = b1(:,2);

clear b1

%% 2:Plot Signals
t0_rest = 10;    % start time for resting
t0_ex   = 20;    % start time for exercise
t0_box  = 15;    % start time for box breathing

% Masks for 60-second windows
mask_rest_60 = (t_rest >= t0_rest) & (t_rest <= t0_rest + 60);
mask_ex_60   = (t_ex   >= t0_ex)   & (t_ex   <= t0_ex   + 60);
mask_box_60  = (t_box  >= t0_box)  & (t_box  <= t0_box  + 60);
figure;
 
% 1) Resting
subplot(3,1,1);
plot(t_rest(mask_rest_60), ecg_rest(mask_rest_60), 'LineWidth', 1.5);
title('Resting (60 s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

% 2) Exercise
subplot(3,1,2);
plot(t_ex(mask_ex_60), ecg_ex(mask_ex_60), 'LineWidth', 1.5);
title('Exercise (60 s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

% 3) Box Breathing
subplot(3,1,3);
plot(t_box(mask_box_60), ecg_box(mask_box_60), 'LineWidth', 1.5);
title('Box Breathing (60 s)', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('ECG (mV)', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

%% 3:Install Signal Processing Toolbox




%% 3:Feature Extraction
%% Calculate sampling rate
Fs = 1/mean(diff(t_rest));

[pks_rest, locs_rest] = findpeaks(ecg_rest, ...
                                  'MinPeakHeight',    0.5, ...    % threshold in mV
                                  'MinPeakDistance',  0.3 * Fs);  % 300 ms * sampling rate
[pks_ex, locs_ex] = findpeaks(ecg_ex, ...
    'MinPeakHeight',0.5, ...
    'MinPeakDistance', 0.3 * Fs);

[pks_box, locs_box] = findpeaks(ecg_box, ...
    'MinPeakHeight', 0.5, ...
    'MinPeakDistance', 0.3 * Fs);

%% Step 4: Data Visual of Heartbeat peaks using findpeaks() [default configs]
%% RESTING HEART RATE
% Plot 2 subplots of the Resting heart rate with the second indicating
% the hearbeat peaks.

% Plot data: x = instance in time; y = wave-amplitude data (use only first 100 values)
% x-value Vector -> t_rest 
% y-value Vector -> ecg_rest

% Formatting Customization configs
windowMargin = 0.5; % border thickness between plot and window edge. 50% margin

% General plot configs
figure;
sgtitle("Resting Heart rate data from ECG");
xTimeVec = 1:100; % Chosen x-range for plots

%% HELPER LAMBDA functions
% Implement function for formatting plot window
% Use 20% borders (from the bottom of the plot to the edge of the window):

% Bottom border: -50% of minimum y-value
% Top border: +50% of maximum y-value

% formula: addBorder = [min(yValArg)*(-50%), max(yValArg)*(50%)]

% Args: yVal - vector of values
% Returns: scaled upper and lower limits for plot clarity
% NOTE: Window margin is a constant throughot all plots
addBorder = @(yVal) [(min(yVal) * (1 + windowMargin)), (max(yVal)  * (1 + windowMargin))]; 


% Implement lambda function for creating offset between plot and peak
% points marker
% Utilizes y-value range to scale offset between peak-marker and plot points.

% formula: clarityOffset = peakVal*(maxY - minY)(percentOffset)

% Args: peakVal - vector of all peaks on `yVals`
%       yVals -  vector of y-values for plot
% Offset between plot and marker
percentOffset = 0.05; % 5% offset on range
clarityOffset = @(peakVal, yVals) peakVal + (max(yVals) - min(yVals))*percentOffset;

%% Plot 1
subplot(2, 1, 1);

yValRest = ecg_rest;
% Use a solid blue line for plot
plot(t_rest, yValRest, "b-");

title("Resting");

% set x and y axis range limits
% Include window margin on plot for for data formatting in ylim
yLimRes = addBorder(yValRest);

ylim(yLimRes);
% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);

%% Plot 2
subplot(2, 1, 2);

% Use a solid blue line for plot
plot(t_rest, yValRest, "b-");
hold on;

% Heart Rate peaks
% Use findpeaks to get points 
[pks_rest_def, locs_rest_def] = findpeaks(yValRest);

% Get time-vec of peak points by indexing `t_rest`
pTimeVecRest = t_rest(locs_rest_def);

% Formatted peakValues
% - Add-in a margin right below the peak indicators for clarity
formatPeakValRes = clarityOffset(pks_rest_def, yLimRes);

% plot peak points
scatter(pTimeVecRest, formatPeakValRes, [], "v", "MarkerEdgeColor", "b", "MarkerFaceColor", "b");
hold off;

% Title subplot
title("Resting with peaks");

% set x and y axis range limits
% Include window margin on plot for for data formatting
ylim(yLimRes);

% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);
grid on;

%% EXERCISE HEART RATE
% Plot 2 subplots of the exercise heart rate with the second indicating
% the hearbeat peaks.

% Plot data: x = instance in time; y = wave-amplitude data
% x-value Vector -> t_ex
% y-value Vector -> ecg_ex

% General plot configs
figure;
sgtitle("Exercise Heart rate data from ECG");

%% Plot 1
subplot(2, 1, 1);

yValsEx = ecg_ex;
% Use a solid blue line for plot
plot(t_ex, yValsEx, "b-");

title("Exercise");
% set x and y axis range limits
% Include window margin on plot for for data formatting in ylim
yLimEx = addBorder(yValsEx);

ylim(yLimEx);
% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);

%% Plot 2
subplot(2, 1, 2);

% Use a solid blue line for plot
plot(t_ex, yValsEx, "b-");
hold on;

% Heart Rate peaks
% Use findpeaks to get points 
[pks_ex_def, locs_ex_def] = findpeaks(yValsEx);

% Get time-vec of peak points by indexing `t_ex`
pTimeVecEx = t_ex(locs_ex_def);

% Formatted peakValues
% - Add-in a margin right below the peak indicators for clarity
formatPeakValEx = clarityOffset(pks_ex_def, yLimEx);

% plot peak points
scatter(pTimeVecEx, formatPeakValEx, [], "v", "MarkerEdgeColor", "b", "MarkerFaceColor", "b");
hold off;

% Title subplot
title("Exercise with peaks");

% set x and y axis range limits
% Include window margin on plot for for data formatting
ylim(yLimEx);

% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);
grid on;

%% 7: Histogram and Boxplot Analysis

% Extract R peak times in 120-130s window for resting
 window_mask_rest = (t_rest(locs_rest) >= 120) & (t_rest(locs_rest) <= 130);
 peak_times_rest = t_rest(locs_rest(window_mask_rest));

 % Calculate intervals
 intervals_rest_manual = peak_times_rest(2: end) - peak_times_rest(1:end-1);

 % Calculate intervals (using diff)
 intervals_rest_diff = diff(peak_times_rest);

 % Convert to heart rate (BPM)
 heart_rate_rest = 60 ./ intervals_rest_diff;

 % Plot histogram
 figure;
 histogram(heart_rate_rest);
 xlabel('Heart Rate (BPM)', 'FontSize', 14);
 ylabel('Count', 'FontSize', 14);
 title('Resting heart rate distribution', 'FontSize', 16);

 %% Extract R peak times in 120-130s window for BOX BREATHING
window_mask_box = (t_box(locs_box) >= 120) & (t_box(locs_box) <= 130);
peak_times_box = t_box(locs_box(window_mask_box));

% Calculate intervals (manual method)
intervals_box_manual = peak_times_box(2:end) - peak_times_box(1:end-1);

% Calculate intervals (using diff)
intervals_box_diff = diff(peak_times_box);

% Convert to heart rate (BPM)
heart_rate_box = 60 ./ intervals_box_diff;

%% Plot histograms for all three conditions
figure;

% Resting histogram
subplot(3,1,1);
histogram(heart_rate_rest);
xlabel('Heart Rate (BPM)', 'FontSize', 14);
ylabel('Count', 'FontSize', 14);
title('Resting Heart Rate Distribution', 'FontSize', 16);
grid on;

% Exercise histogram
subplot(3,1,2);
histogram(heart_rate_ex);
xlabel('Heart Rate (BPM)', 'FontSize', 14);
ylabel('Count', 'FontSize', 14);
title('Exercise Heart Rate Distribution', 'FontSize', 16);
grid on;

% Box Breathing histogram
subplot(3,1,3);
histogram(heart_rate_box);
xlabel('Heart Rate (BPM)', 'FontSize', 14);
ylabel('Count', 'FontSize', 14);
title('Box Breathing Heart Rate Distribution', 'FontSize', 16);
grid on;

%% Create boxchart with combined data
% Combine all heart rate data
combinedData = [heart_rate_rest; heart_rate_ex; heart_rate_box];

% Create group labels
group = [ones(length(heart_rate_rest), 1); 
         2*ones(length(heart_rate_ex), 1); 
         3*ones(length(heart_rate_box), 1)];

% Plot boxchart
figure;
boxchart(group, combinedData);
xlabel('Condition', 'FontSize', 14);
ylabel('Heart Rate (BPM)', 'FontSize', 14);
title('Heart Rate Comparison Across Conditions', 'FontSize', 16);
xticklabels({'Resting', 'Exercise', 'Box Breathing'});
grid on;
set(gca, 'FontSize', 14);