%File Name: Master.m
%Dependencies: Signal processing toolbox
% Steps 1-3:Noah Nir
% Steps 4-6:Will Dakare
% Steps 7-9:Sam Yeboah-Asi

% William: A struggle we had were finding a efficient way to make the plot window
% clear. We solved this by making helper functions. 

% Noah: Another issue we had
% was to properly load the mat files in order to use them properly. This
% was solved by clearing b1 after each load.

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


%% 7 — Histogram and Boxplot (match project PDF style)

%  Extract R peak times within a 120–130 s window 
window_rest = (t_rest(locs_rest) >= 120) & (t_rest(locs_rest) <= 130);
window_ex   = (t_ex(locs_ex)   >= 120) & (t_ex(locs_ex)   <= 130);

peak_times_rest = t_rest(locs_rest(window_rest));
peak_times_ex   = t_ex(locs_ex(window_ex));

%  Compute intervals (R–R intervals) 
intervals_rest = diff(peak_times_rest);
intervals_ex   = diff(peak_times_ex);

%  Convert intervals to heart rate (BPM) 
heart_rate_rest = 60 ./ intervals_rest;
heart_rate_ex   = 60 ./ intervals_ex;

%%  HISTOGRAMS 
figure;

histogram(heart_rate_rest);
xlabel('Heart Rate (BPM)', 'FontSize', 14);
ylabel('Count', 'FontSize', 14);
title('Resting Heart Rate Distribution', 'FontSize', 16);
grid on;
hold on;
% Exercise histogram
histogram(heart_rate_ex);
xlabel('Heart Rate (BPM)', 'FontSize', 14);
ylabel('Count', 'FontSize', 14);
title('Exercise Heart Rate Distribution', 'FontSize', 16);
grid on;
legend("Resting Heart Rate", "Exercise Heart Rate Distribution");
hold off;
%%  BOXCHART (combined)
combined = [heart_rate_rest; heart_rate_ex];
groups = [ones(length(heart_rate_rest),1); ...
          2*ones(length(heart_rate_ex),1)];

figure;
boxchart(groups, combined);
xticklabels({'Resting','Exercise'});
ylabel('Heart Rate (bpm)');
title('Heart Rate Comparison Across Conditions');
set(gca,'FontSize',14);
grid on;

%% 8: Compute the derivative of different phases of the signal

% Find T peaks in resting data
[pks_T_rest, locs_T_rest] = findpeaks(ecg_rest, ...
    'MinPeakHeight', 0.15, ...       
    'MinPeakDistance', 0.15 * Fs);

% Remove R peaks from T peak list
is_T_peak = true(length(locs_T_rest), 1);
for i = 1:length(locs_T_rest)
    % Check if this peak is too close to any R peak
    time_diff = abs(t_rest(locs_T_rest(i)) - t_rest(locs_rest));
    if any(time_diff < 0.1) 
        is_T_peak(i) = false;
    end
end
% Keep only the true T peaks
locs_T_rest = locs_T_rest(is_T_peak);
pks_T_rest = pks_T_rest(is_T_peak);

% Initialize arrays to store speeds
R_speeds = zeros(length(locs_rest), 1);
T_speeds = zeros(length(locs_T_rest), 1);

% Compute speed for R pulses using a for loop
for i = 1:length(locs_rest)
    % Get indices 10 samples before and after R peak
    idx_before = locs_rest(i) - 10;
    idx_after = locs_rest(i) + 10;
    
    % Make sure we don't go out of bounds
    if idx_before < 1
        idx_before = 1;
    end
    if idx_after > length(ecg_rest)
        idx_after = length(ecg_rest);
    end
    
    % Compute absolute voltage change
    voltage_change = abs(ecg_rest(idx_after) - ecg_rest(idx_before));
    
    % Compute time difference
    time_change = t_rest(idx_after) - t_rest(idx_before);
    
    % Speed = voltage change / time change (mV/s)
    R_speeds(i) = voltage_change / time_change;
end

% Compute speed for T pulses using a for loop
for i = 1:length(locs_T_rest)
    % Get indices 10 samples before and after T peak
    idx_before = locs_T_rest(i) - 10;
    idx_after = locs_T_rest(i) + 10;
    
    % Make sure we don't go out of bounds
    if idx_before < 1
        idx_before = 1;
    end
    if idx_after > length(ecg_rest)
        idx_after = length(ecg_rest);
    end
    
    % Compute absolute voltage change
    voltage_change = abs(ecg_rest(idx_after) - ecg_rest(idx_before));
    
    % Compute time difference
    time_change = t_rest(idx_after) - t_rest(idx_before);
    
    % Speed = voltage change / time change (mV/s)
    T_speeds(i) = voltage_change / time_change;
end

% Plot histogram comparing R pulse speeds vs T pulse speeds
%  Histogram for R pulse speeds 
figure;
histogram(R_speeds, 'BinWidth', 0.5);   
title('R Pulse Speed Distribution');
xlabel('Pulse Speed (mV/s)');
ylabel('Count');
set(gca, 'FontSize', 14);
grid on;

%  Histogram for T pulse speeds 
figure;
histogram(T_speeds, 'BinWidth', 0.5);
title('T Pulse Speed Distribution');
xlabel('Pulse Speed (mV/s)');
ylabel('Count');
set(gca, 'FontSize', 14);
grid on;

%% 9: Run a two-sided t-test using ttest2

% Run two t-tests using a for loop with conditional statements
for test_num = 1:2
    if test_num == 1
        % Test 1: Compare heart rate during rest vs exercise
        [h, p] = ttest2(heart_rate_rest, heart_rate_ex);
        
        fprintf('\n Test 1: Resting vs Exercise Heart Rate ===\n');
        fprintf('Null Hypothesis: The two datasets are NOT different\n');
        
        if h == 1
            fprintf('Result: We REJECT the null hypothesis (h = %d)\n', h);
            fprintf('Conclusion: There IS statistically significant evidence that\n');
            fprintf('            resting and exercise heart rates are DIFFERENT.\n');
            fprintf('p-value: %.6f\n', p);
        else
            fprintf('Result: We CANNOT reject the null hypothesis (h = %d)\n', h);
            fprintf('Conclusion: There is NOT enough evidence to claim the heart rates are different.\n');
            fprintf('p-value: %.6f\n', p);
        end
        
    else
        % Test 2: Compare two different periods of the same resting data
        % Split resting heart rate data into two halves
        midpoint = floor(length(heart_rate_rest) / 2);
        rest_period1 = heart_rate_rest(1:midpoint);
        rest_period2 = heart_rate_rest(midpoint+1:end);
        
        [h, p] = ttest2(rest_period1, rest_period2);
        
        fprintf('\n Test 2: First Half vs Second Half of Resting Heart Rate ===\n');
        fprintf('Null Hypothesis: The two datasets are NOT different\n');
        
        if h == 1
            fprintf('Result: We REJECT the null hypothesis (h = %d)\n', h);
            fprintf('Conclusion: There IS statistically significant evidence that\n');
            fprintf('            the two resting periods are DIFFERENT.\n');
            fprintf('p-value: %.6f\n', p);
        else
            fprintf('Result: We CANNOT reject the null hypothesis (h = %d)\n', h);
            fprintf('Conclusion: There is NOT enough evidence to claim the two resting periods are different.\n');
            fprintf('            This makes sense - both periods are from the same resting condition!\n');
            fprintf('p-value: %.6f\n', p);
        end
    end
end

