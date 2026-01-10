%% HEADER DOCUMENTATION
%File Name: Master.m
% Dependencies: 
% - Signal processing toolbox
% - Statistical and Machine Learning toolbox

% Steps 1-3:Noah Nir
% Steps 4-6:Will Dakare
% Steps 7-9:Sam Yeboah-Asi

% William: A struggle we had were finding a efficient way to make the plot window
% clear. We solved this by making helper functions. There were also issues 
% with overlapping the 3-wave cycles on each other to measure trigger integrity. 

% Noah: Another issue we had
% was to properly load the mat files in order to use them properly. This
% was solved by clearing b1 after each load.

% Limitations:
% - Script sequence has limited modularity for plotting because some
% variables are arbitrarily chosen for the dataset and observed trends

% Improvements:
% - Creating an intelligent seqeunce that implements modulated data
% analyis and visualization choices based on the data signal type.

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

% Limit 1:
xlim([t0_rest t0_rest + 60]);
y = ecg_rest(mask_rest_60); 
yr = max(y) - min(y);
ylim([min(y) - 0.1*yr, max(y) + 0.1*yr]);

% 2) Exercise
subplot(3,1,2);
plot(t_ex(mask_ex_60), ecg_ex(mask_ex_60), 'LineWidth', 1.5);
title('Exercise (60 s)', 'FontSize', 16);
ylabel('ECG (mV)', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

% Limit 2:
xlim([t0_ex t0_ex + 60]);
y = ecg_ex(mask_ex_60);
yr = max(y) - min(y);
ylim([min(y) - 0.1*yr, max(y) + 0.1*yr]);


% 3) Box Breathing
subplot(3,1,3);
plot(t_box(mask_box_60), ecg_box(mask_box_60), 'LineWidth', 1.5);
title('Box Breathing (60 s)', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 14);
ylabel('ECG (mV)', 'FontSize', 14);
set(gca, 'FontSize', 14);
grid on;

% Limit 3:
xlim([t0_box t0_box + 60]);
y = ecg_box(mask_box_60);
yr = max(y) - min(y);
ylim([min(y) - 0.1*yr, max(y) + 0.1*yr]);

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
window_margin = 0.5; % border thickness between plot and window edge. 50% margin

% General plot configs
figure;
sgtitle("Resting Heart rate data from ECG");
x_time_vec = 1:100; % Chosen x-range for plots

%% HELPER functions
% Implement function for formatting plot window
% Use 20% borders (from the bottom of the plot to the edge of the window):

% Bottom border: -50% of minimum y-value
% Top border: +50% of maximum y-value

% formula: addBorder = [min(yValArg)*(-50%), max(yValArg)*(50%)]

% Args: yVal - vector of values
% Returns: scaled upper and lower limits for plot clarity
% NOTE: Window margin is a constant throughot all plots
add_border = @(y_val) [(min(y_val) * (1 + window_margin)), (max(y_val)  * (1 + window_margin))]; 


% Implement lambda function for creating offset between plot and peak
% points marker
% Utilizes y-value range to scale offset between peak-marker and plot points.

% formula: clarityOffset = peakVal*(maxY - minY)(percentOffset)

% Args: peakVal - vector of all peaks on `yVals`
%       yVals -  vector of y-values for plot
% Offset between plot and marker
percent_offset = 0.10; % 10% offset on range
clarity_offset = @(peak_val, y_vals) peak_val + (max(y_vals) - min(y_vals))*percent_offset;


function fit_window = return_window_arr(time_vec, y_vec)
% Formats peakValues
% Window size: x-vals (101 to 110)
% - Add-in a margin right below the peak indicators for clarity
% - Designate args for `clarityOffset` to be specific to the window range
% Args:
%   timeVec - vector of double-values for plot
%   yVec - vector of y-values for plot 
% Returns:
%   array of 2 values [min max] fit the 101 to 110 window

% compare error, because find only works with exact values: 0.3 != 0.30000
indx_start = find(abs(time_vec - 101) < 0.0001); 
indx_end = find(abs(time_vec - 110) < 0.0001);

range_of_vals = y_vec(indx_start: indx_end); % y-Values for a window of values on the x-axis.
fit_window = [min(range_of_vals)  max(range_of_vals)];
end

%% Plot 1
subplot(2, 1, 1);

y_val_rest = ecg_rest;
% Use a solid blue line for plot
plot(t_rest, y_val_rest, "b-");

title("Resting");

% set x and y axis range limits
% Include window margin on plot for for data formatting in ylim
y_lim_res = add_border(y_val_rest);

ylim(y_lim_res);
% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);

%% Plot 2
subplot(2, 1, 2);

% Use a solid blue line for plot
plot(t_rest, y_val_rest, "b-");
hold on;

% Heart Rate peaks
% Use findpeaks to get points 
[pks_rest_def, locs_rest_def] = findpeaks(y_val_rest);

% Get time-vec of peak points by indexing `t_rest`
p_time_vec_rest = t_rest(locs_rest_def);

fit_window_rest = return_window_arr(t_rest, y_val_rest);
format_peak_val_res = clarity_offset(pks_rest_def, fit_window_rest);

% plot peak points
scatter(p_time_vec_rest, format_peak_val_res, [], "v", "MarkerEdgeColor", "b", "MarkerFaceColor", "b");
hold off;

% Title subplot
title("Resting with peaks");

% set x and y axis range limits
% Include window margin on plot for for data formatting
ylim(y_lim_res);

% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);
grid on;

% save plot
saveas(gcf, "resting_with_default_findpeaks.png");

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

y_vals_ex = ecg_ex;
% Use a solid blue line for plot
plot(t_ex, y_vals_ex, "b-");

title("Exercise");
% set x and y axis range limits
% Include window margin on plot for for data formatting in ylim

ylim([-0.5 0.8]);
% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);

%% Plot 2
subplot(2, 1, 2);

% Use a solid blue line for plot
plot(t_ex, y_vals_ex, "b-");
hold on;

% Heart Rate peaks
% Use findpeaks to get points 
[pks_ex_def, locs_ex_def] = findpeaks(y_vals_ex);

% Get time-vec of peak points by indexing `t_ex`
p_time_vec_ex = t_ex(locs_ex_def);

% Formatted peakValues
% - Add-in a margin right below the peak indicators for clarity
% - Designate range on y-vector to make borders specific to the window
fit_window_ex = return_window_arr(t_ex, y_vals_ex);
format_peak_val_ex = clarity_offset(pks_ex_def, fit_window_ex); 

% plot peak points
scatter(p_time_vec_ex, format_peak_val_ex, [], "v", "MarkerEdgeColor", "b", "MarkerFaceColor", "b");
hold off;

% Title subplot
title("Exercise with peaks");

% set x and y axis range limits
% Include window margin on plot for for data formatting
ylim([-0.5 0.8]);

% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);
grid on;

% save plot
saveas(gcf, "exercise_with_default_findpeaks.png");

%% STEP 5: Data Visual of Heartbeat peaks using findpeaks() [MinPeakHeight = 0.3]
%% RESTING HEART RATE
% Plot 2 subplots of the Resting heart rate with the second indicating
% the hearbeat peaks.

% Plot data: x = instance in time; y = wave-amplitude data (use only first 100 values)
% x-value Vector -> t_rest 
% y-value Vector -> ecg_rest

% Formatting Customization configs
window_margin = 0.5; % border thickness between plot and window edge. 50% margin

% General plot configs
figure;
sgtitle("Resting heart rate data from ECG");
%% Plot 1
subplot(2, 1, 1);

y_val_rest = ecg_rest;
% Use a solid blue line for plot
plot(t_rest, y_val_rest, "b-");

title("Resting");

% set x and y axis range limits
% Include window margin on plot for for data formatting in ylim
y_lim_res = add_border(y_val_rest);

ylim(y_lim_res);
% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);

%% Plot 2
subplot(2, 1, 2);

% Use a solid blue line for plot
plot(t_rest, y_val_rest, "b-");
hold on;

% Heart Rate peaks
% Use findpeaks to get points 
[pks_rest_def, locs_rest_def] = findpeaks(y_val_rest, "MinPeakHeight", 0.3);

% Get time-vec of peak points by indexing `t_rest`
p_time_vec_rest = t_rest(locs_rest_def);

fit_window_rest = return_window_arr(t_rest, y_val_rest);
format_peak_val_res = clarity_offset(pks_rest_def, fit_window_rest);

% plot peak points
scatter(p_time_vec_rest, format_peak_val_res, [], "v", "MarkerEdgeColor", "b", "MarkerFaceColor", "b");
hold off;

% Title subplot
title("Resting with peaks - Marked R-wave");

% set x and y axis range limits
% Include window margin on plot for for data formatting
ylim(y_lim_res);

% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);
grid on;

% save plot
saveas(gcf, "resting_with_formatted_findpeaks.png");

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

y_vals_ex = ecg_ex;
% Use a solid blue line for plot
plot(t_ex, y_vals_ex, "b-");

title("Exercise");
% set x and y axis range limits
% Include window margin on plot for for data formatting in ylim

ylim([-0.5 0.8]);
% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);

%% Plot 2
subplot(2, 1, 2);

% Use a solid blue line for plot
plot(t_ex, y_vals_ex, "b-");
hold on;

% Heart Rate peaks
% Use findpeaks to get points 
[pks_ex_def, locs_ex_def] = findpeaks(y_vals_ex, "MinPeakHeight", 0.3);

% Get time-vec of peak points by indexing `t_ex`
p_time_vec_ex = t_ex(locs_ex_def);

% Formatted peakValues
% - Add-in a margin right below the peak indicators for clarity
% - Designate range on y-vector to make borders specific to the window
fit_window_ex = return_window_arr(t_ex, y_vals_ex);
format_peak_val_ex = clarity_offset(pks_ex_def, fit_window_ex); 

% plot peak points
scatter(p_time_vec_ex, format_peak_val_ex, [], "v", "MarkerEdgeColor", "b", "MarkerFaceColor", "b");
hold off;

% Title subplot
title("Exercise with peaks - Marked R-wave");

% set x and y axis range limits
% Include window margin on plot for for data formatting
ylim([-0.5 0.8]);

% x-value window from 101 to 100 for data visual clarity
xlim([101 110]);
grid on;

% save plot
saveas(gcf, "exercise_with_formatted_findpeaks.png");

%% Step 6: Visualize trigger redundancy in each wave-instance cycle
% Plot the same wave instances over intervals where all R, T and P-waves
% occur. This will show the integrity of the signal trigger for each wave
% occurence:

%% Step 6.1: Find the interval of 3-wave cycle using diff()
% Use diff to find the interval of occurence then take the average:
% For resting:
sum_rest_cycle = diff(p_time_vec_rest);
avg_rest_cycle = mean(sum_rest_cycle);

% For exercise
sum_ex_cycle = diff(p_time_vec_ex);
avg_ex_cycle = mean(sum_ex_cycle);

% For box breathing
p_time_vec_box = t_box(locs_box);
sum_box_cycle = diff(p_time_vec_box);
avg_box_cycle = mean(sum_box_cycle);

%% Step 6.2: Plot overlayed plot with formatted transparency
% For resting:
figure;% parameters (tweak if you want)
preSec  = 0.20;   % seconds before R peak
postSec = 0.60;   % seconds after R peak

Fs    = 1/mean(diff(t_rest));   % sampling frequency
preS  = round(preSec * Fs);
postS = round(postSec * Fs);

figure; hold on; box on;
set(gcf,'Name','Aligned beats — Resting','NumberTitle','off');

%% Step 6.2: Plot overlayed plot with formatted transparency
% For resting:
for k = 1:numel(locs_rest)
    idx = locs_rest(k); % sample index of R peak

    % define window edges
    idxStart = max(1, idx - preS);      % truncate if before start
    idxEnd   = min(numel(ecg_rest), idx + postS); % truncate if past end

    % Build y-vector from window edges
    y_vector = ecg_rest(idxStart:idxEnd);

    % baseline-correct using pre-trigger samples
    numPreSamples = min(preS, length(y_vector));
    y_vector = y_vector - mean(y_vector(1:numPreSamples));

    % Build x-vector to match y-vector length
    t_rel = (-numPreSamples:(length(y_vector)-numPreSamples-1)) / Fs; % relative to R-peak

    % plot overlapping graph with 50% formatted transparency
    set(gcf,'Renderer','opengl');       % helps MATLAB respect alpha
    plot(t_rel, y_vector, 'Color', [1 0 0 0.25], 'LineWidth', 0.8);
    hold on;
end

xlabel('Time from R-peak (s)');
ylabel('Amplitude (a.u.)');
xlim([-preSec postSec]);
title('Overlapped beats centered on R trigger — Resting');
hold off;

% for box breathing
figure;
for k = 1:numel(locs_box)
    idx = locs_box(k); % sample index of R peak

    % define window edges
    idxStart = max(1, idx - preS);      
    idxEnd   = min(numel(ecg_box), idx + postS); 

    % Build y-vector from window edges
    y_vector = ecg_box(idxStart:idxEnd);

    % baseline-correct using pre-trigger samples
    numPreSamples = min(preS, length(y_vector));
    y_vector = y_vector - mean(y_vector(1:numPreSamples));

    % Build x-vector to match y-vector length
    t_rel = (-numPreSamples:(length(y_vector)-numPreSamples-1)) / Fs; % relative to R-peak

    % plot overlapping graph with 50% formatted transparency
    set(gcf,'Renderer','opengl');       
    plot(t_rel, y_vector, 'Color', [1 0 0 0.25], 'LineWidth', 0.8);
    hold on;
end

xlabel('Time from R-peak (s)');
ylabel('Amplitude (a.u.)');
xlim([-preSec postSec]);
title('Overlapped beats centered on R trigger — Box');
hold off;

% for exercise
figure;
for k = 1:numel(locs_ex)
    idx = locs_ex(k); % sample index of R peak

    % define window edges
    idxStart = max(1, idx - preS);      
    idxEnd   = min(numel(ecg_ex), idx + postS); 

    % Build y-vector from window edges
    y_vector = ecg_ex(idxStart:idxEnd);

    % baseline-correct using pre-trigger samples
    numPreSamples = min(preS, length(y_vector));
    y_vector = y_vector - mean(y_vector(1:numPreSamples));

    % Build x-vector to match y-vector length
    t_rel = (-numPreSamples:(length(y_vector)-numPreSamples-1)) / Fs; % relative to R-peak

    % plot overlapping graph with 50% formatted transparency
    set(gcf,'Renderer','opengl');       
    plot(t_rel, y_vector, 'Color', [0 0 1 0.25], 'LineWidth', 0.8); % blue for exercise
    hold on;
end

xlabel('Time from R-peak (s)');
ylabel('Amplitude (a.u.)');
xlim([-preSec postSec]);
title('Overlapped beats centered on R trigger — Exercise');
hold off;

%% 7: Histogram and Boxplot Analysis
% Extract R peak times in 120-130s window for RESTING
window_mask_rest = (t_rest(locs_rest) >= 120) & (t_rest(locs_rest) <= 130);
peak_times_rest = t_rest(locs_rest(window_mask_rest));

% Calculate intervals (manual method)
intervals_rest_manual = peak_times_rest(2:end) - peak_times_rest(1:end-1);

% Calculate intervals (using diff)
intervals_rest_diff = diff(peak_times_rest);

% Convert to heart rate (BPM)
heart_rate_rest = 60 ./ intervals_rest_diff;

%% Extract R peak times in 120-130s window for EXERCISE
window_mask_ex = (t_ex(locs_ex) >= 120) & (t_ex(locs_ex) <= 130);
peak_times_ex = t_ex(locs_ex(window_mask_ex));

% Calculate intervals (manual method)
intervals_ex_manual = peak_times_ex(2:end) - peak_times_ex(1:end-1);

% Calculate intervals (using diff)
intervals_ex_diff = diff(peak_times_ex);

% Convert to heart rate (BPM)
heart_rate_ex = 60 ./ intervals_ex_diff;

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
figure;
hold on;
histogram(R_speeds, 'FaceColor', 'b', 'FaceAlpha', 0.6, 'DisplayName', 'R Pulse Speed');
histogram(T_speeds, 'FaceColor', 'r', 'FaceAlpha', 0.6, 'DisplayName', 'T Pulse Speed');
hold off;

xlabel('Pulse Speed (mV/s)', 'FontSize', 14);
ylabel('Count', 'FontSize', 14);
title('Speed Comparison: R Pulse vs T Pulse (Resting Data)', 'FontSize', 16);
legend('FontSize', 12, 'Location', 'best');
grid on;
set(gca, 'FontSize', 14);

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

