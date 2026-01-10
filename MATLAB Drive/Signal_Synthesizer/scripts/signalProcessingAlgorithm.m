% READ ME:
% signalProcessingAlgorithm.m - Processes data, makes comparisons inferences
% and conclusions on signal quality, robustness and classification of
% signal's quality disruption resistance to noise (measure of noise's effect on signal quality).
%
% Team Members: [Noah Nir(Data-Structure), Samuel Yeboah-Asi(Data Visualization), William Dakare(Algorithm Developer)]
% Date: November 17-23 2025
%
% Description:
% This MATLAB script loads signal CSVs, computes quality metrics (RMS, peak-to-peak, zero crossings), 
% measures percent degradation caused by noise, runs simple filtering experiments, produces figures/plots, 
% and saves summary outputs. It also includes lightweight tests and helper methods for processing 
% columns as structured fields (e.g., allSignals.Clean.Sine1Hz).

% Inputs: Expected files at these paths
%   - These parts were technically created in this scrips but the code for its
%   creation is in helperScript
% ../data/signals/Clean/Clean.csv
% ../data/signals/LowNoise/lowNoise.csv
% ../data/signals/HighNoise/largeNoise.csv

% Outputs: Printed Analysis and Comparison graphs:
    % results/figures/scenario_comparison.png
    % results/figures/frequency_comparison.png
    % results/figures/quality_metrics.png
    % results/figures/filtering_effect.png
    % results/figures/signal_dashboard.png
    % signalAnalysisResults.mat (workspace variables)
    % signalAnalysisSummary.csv

%% Creats Files:
helperScript.projectCsvMaker;

%% Algorithm Development: William Dakare
%% EXTRACT ALL SIGNAL DATA from CSV files

% SET-UP
% Target Path: "C:\Users\dakar\MATLAB Drive\Signal_Synthesizer\data\signals\(Clean|LowNoise|HighNoise)"
% This script is in the "Signal_Synthesizer" directory. 
% To navigate from there to the target path create a relative base path
% Current path is "C:\Users\dakar\MATLAB Drive\Signal_Synthesizer\scripts

% I'll be using a structure to save each signal-types's corresponding
% category as attributes (For instance: clean.Sine10Hz).
% This approach is chosen because column header names can be formatted into
% attribute names which will increase vectorization.

%% CREATE MASTER STRUCTURE for ALL SIGNALS and Categories
% Purpose: Simplified and cohesive structure for indexing and data
% retrieval

% Make a dictionary of Dir : File pairs i.e. "Clean" : "Clean.csv"
% Create a dictionary: To be used in Structure object creation
allSignalsDirFilePair = containers.Map(["Clean", "LowNoise", "HighNoise"], ["Clean.csv", "lowNoise.csv", "largeNoise.csv"]);

allKeys = keys(allSignalsDirFilePair);

% Store each scenario as a table inside allSignals
allSignals = struct();

for index = 1:length(allKeys)
    dirName = allKeys{index};                    
    fileName = allSignalsDirFilePair(dirName);
    
    filePath = fullfile(relativePathBase, dirName, fileName);

    if ~isfile(filePath) % check if file exists
       error("File NOT FOUND. Go to ../README.md");
    end
    
     % Read CSV into table
    signalTable = readtable(filePath);
    
    % Make column names valid for table variable names
    signalTable.Properties.VariableNames = matlab.lang.makeValidName(signalTable.Properties.VariableNames);
    
    % Store the entire table for this scenario 
    allSignals.(dirName) = signalTable;
end



%% TEST for Data Validity & Table Length
% Get Clean -> Sine1Hz (first Column), and Triangle10Hz (last column) vector

% TEST:
% To test set arg for `shouldTest` to (logical 1 = True) | or (0 for false)
shouldTest(allSignals, false); % set to true to test
% Test results: PASS


% At this point all signal columns are in the structure allSignals in
% form-> allSignals.("Clean" | "LowNoise" | "HighNoise").("sine1Hz"..."tri10Hz")

%% CALCULATE SIGNAL QUALITY for ALL SIGNALS and SUB_CATEGORIES
% processSignals(processType) in helperScript.m processes all tables, by column, in the
% data folder and returns a structure-object that can be indexed by
% directory name(e.g, Clean, LowNoise) and column(e.g, Sine1Hz, Triangle5Hz)

%% TEST processSignals(processType) in helperScript.m with testProcessor(TrueOrFalse)
% testProcessor(TrueOrFalse) tests processSignals 
helperScript.testProcessor(false); % pass true to test

%% Root Mean Square Comparison Across Signals 

% Structure = allSignals(3-letter abreviation in uppercase: RMS, PTP, ZCS)->
% .(processType: "rms"|"peaktopeak"|"zerocross")
% .(directory: "Clean"...)
% .(column: "Sine1Hz"...)

% Root Mean Square for all signal types across all frequencies
allSignalsRMS = helperScript.processSignals("rms");

% Peak to Peak Amplitude for all signal types across all frequencies
allSignalsPTP = helperScript.processSignals("peaktopeak");

% Zero Crossings (sign change) for all signal types across all frequencies
allSignalsZCS = helperScript.processSignals("zerocross");

% Store RMS for each signal in a matrix
% Matrix has column format:
% Sine1Hz | Sine5Hz | Sine10Hz | Square1Hz | Square5Hz | Square10Hz | Triangle1Hz | Triangle5Hz | Triangle10Hz
fieldNames = ["Sine1Hz", "Sine5Hz", "Sine10Hz", ...
            "Square1Hz", "Square5Hz", "Square10Hz", ...
            "Triangle1Hz", "Triangle5Hz", "Triangle10Hz"];

allRmsMat = zeros(9, 3); % Memory preallocation
row = 1;
col = 1;
for signal = ["Clean", "LowNoise", "HighNoise"]
    for field = fieldNames
       noisySignalRms = allSignalsRMS.rms.(signal).(field);
       allRmsMat(row, col) = noisySignalRms;
       row = row + 1; % go down the vector column
    end
    row = 1; % Reset row value after each column 
    col = col + 1;
end

% Finished building allRmsMat with cols: Clean | LowNoise | HighNoise; and
% rows: Sine1Hz | ... | Triangle10Hz

% Calculate Percentage Degradation on Noisy values
rpercentDegMat = zeros(9, 2); % Preallocation. Also where frequencies will be stored
% Columns: LowNoise | HighNoise

cleanRms = allRmsMat(:,1);
lowNoiseRms = allRmsMat(:,2);
highNoiseRms = allRmsMat(:,3);

%% CALCULATING AND COMPARING % of DEGRADATION using Weighted Sums
% Percent Degradation Formula: ((clean_value - noisy_value) / clean_value) × 100
lowNoisePercentDeg = percentDegradation(cleanRms, lowNoiseRms);
highNoisePercentDeg = percentDegradation(cleanRms, highNoiseRms);

% Reshape row-vectors (1x9) low/highNoisePercentDeg into (3 x 3) matrix for
% the structure: 
% [Sine1Hz... Sine10Hz;
% Square1z ...;
% ... Triangle10Hz]
% use this for row-wise comparison and analysis 
     
lowDegMatrix  = reshape(lowNoisePercentDeg, 3, 3);  
highDegMatrix = reshape(highNoisePercentDeg,3,3);

fprintf("--- Noise Resistance Classification and Most Robust Signal Identification---\n");

%% PERCENTAGE DEGRADATION COMPARISION (Column-Wise)

% Find most affected frequency (average across waveforms)
freqLabels = ["1Hz","5Hz","10Hz"];

meanLowPerFreq  = mean(lowDegMatrix); %Get means 1 x 3
meanHighPerFreq = mean(highDegMatrix);

medianLowPerFreq  = median(lowDegMatrix, 1);
medianHighPerFreq = median(highDegMatrix, 1);

% Combined score (weighted sum: 60% mean, 40% median) - higher = worse |
% lowest = best/most robust 
combinedLowScore  = 0.6 * meanLowPerFreq  + 0.4 * medianLowPerFreq;
combinedHighScore = 0.6 * meanHighPerFreq + 0.4 * medianHighPerFreq;

% Find best & worst index 
maxLow  = max(combinedLowScore);
maxHigh = max(combinedHighScore);

% Find index postion of the maxLow/High in combinedLowScore
% This index maps on directly to freqLabels, to locate freq. most affected
% by noise
lowMaxPosition = find(maxLow == combinedLowScore);
highMaxPosition = find(maxLow == combinedHighScore);

if ~(highMaxPosition == lowMaxPosition)
    fprintf('\nMost affected frequency (avg across waveforms):\n');
    fprintf('Low Noise  : %s', freqLabels(lowMaxPosition));
    fprintf('High Noise : %s', freqLabels(highMaxPosition));
else
    fprintf("Frequency affected most by noise in both scenarios..." + ...
        "(high and low noise): %s", freqLabels(lowMaxPosition));
end
disp(" ");

%% SIGNAL ROBUSTNESS COMPARISION (Row-Wise)
% Very similar to % Degradation comparison sequence but row-wise for signal
% type and we check for minimum values for max robustness
signalLabels = ["Sine", "Square", "Triangle"]; % Adjust according to your row order

% Row-wise mean & median
meanLowPerSignal  = mean(lowDegMatrix, 2); %Get means per signal (row)
meanHighPerSignal = mean(highDegMatrix, 2);

medianLowPerSignal  = median(lowDegMatrix, 2);
medianHighPerSignal = median(highDegMatrix, 2);

% Same formula as above
combinedLowScoreSig  = 0.6 * meanLowPerSignal  + 0.4 * medianLowPerSignal;
combinedHighScoreSig = 0.6 * meanHighPerSignal + 0.4 * medianHighPerSignal;

% Find best index 
minLowSig  = min(combinedLowScoreSig);
minHighSig = min(combinedHighScoreSig);

% Find index postion of the minLowSig/HighSig in combinedLowScoreSig
% This index maps on directly to signalLabels, to locate most robust
% signals under noise constriants
lowMaxPositionSig = find(minLowSig == combinedLowScoreSig);
highMaxPositionSig = find(minHighSig == combinedHighScoreSig);

% Display results
if ~(highMaxPositionSig == lowMaxPositionSig)
    fprintf('\nMost robust signal type (across low and large noise scenarios):\n');
    fprintf('Low Noise  : %s\n', signalLabels(lowMaxPositionSig));
    fprintf('High Noise : %s\n', signalLabels(highMaxPositionSig));
else
    fprintf("Most robust signal in Either scenario is: %s Wave\n", ...
        signalLabels(lowMaxPositionSig));
end
disp(" ");


%% SIGNAL STATISTICS
scenarioNames = ["Low Noise", "High Noise"];
signalLabels  = ["Sine", "Square", "Triangle"]; % adjust per row order

% Pre-allocate stats arrays (rows = signals, cols = scenarios)
numSignals   = size(lowDegMatrix, 1);
numScenarios = 2; % Low vs High Noise

meanVals = zeros(numSignals, numScenarios);
stdVals  = zeros(numSignals, numScenarios);
varVals  = zeros(numSignals, numScenarios);

% Loop over each signal (row-wise)
for signal = 1:numSignals
    % Low Noise
    signalLow = lowDegMatrix(signal, :);
    meanVals(signal, 1) = mean(signalLow);       % Mean value
    stdVals(signal, 1)  = std(signalLow);        % Standard deviation
    varVals(signal, 1)  = var(signalLow);        % Variance

    % High Noise
    signalHigh = highDegMatrix(signal, :);
    meanVals(signal, 2) = mean(signalHigh);
    stdVals(signal, 2)  = std(signalHigh);
    varVals(signal, 2)  = var(signalHigh);
end

fprintf('--- Signal Statistics per Scenario ---\n');
for signal = 1:numSignals
    fprintf('%s:\n', signalLabels(signal));
    for scenario = 1:numScenarios
        fprintf('  %s - Mean: %.4f, Std: %.4f, Var: %.4f\n', ...
                scenarioNames(scenario), meanVals(signal, scenario), stdVals(signal, scenario), varVals(signal, scenario));
    end
end
disp(" ");

%% Moving Average Filtering
windowSize = 5; % smoothing window

% Example: Apply to first high-noise signal (e.g., Sine)
sineSigIndex = 1; % Sine row
triSigIndex = 3; % Triangle wave row

sineSignalHigh = highDegMatrix(sineSigIndex, :);
triSignalHigh = highDegMatrix(triSigIndex, :);

% RMS before filtering (sine and triangle)
rmsSineBefore = rms(sineSignalHigh);
rmsTriBefore = rms(triSignalHigh);

% Apply moving average filter
signalSineFiltered = movmean(sineSignalHigh, windowSize);
signalTriFiltered = movmean(sineSignalHigh, windowSize);

% RMS after filtering (sine and triangle)
rmsSineAfter = rms(signalSineFiltered);
rmsTriAfter = rms(signalTriFiltered);

% Improvement percentage
sineImprovementPercent = ((rmsSineBefore - rmsSineAfter) / rmsSineBefore) * 100;
triImprovementPercent = ((rmsTriBefore - rmsTriAfter) / rmsTriBefore) * 100;

fprintf('--- Filtering Results for %s (High Noise) ---\n', signalLabels(sineSigIndex));
fprintf("\t\t  Sine Wave | Triangle Wave")
disp(" ")
fprintf('RMS before filtering: %.4f %.4f\n', rmsSineBefore, rmsTriBefore);
fprintf('RMS after filtering : %.4f %.4f\n', rmsSineAfter, rmsTriAfter);
fprintf('Improvement:\t\t %.2f%% %.2f%%\n', sineImprovementPercent, triImprovementPercent);
disp(" ");

%% Save Analysis Outputs
% Outputs all go to the `results` directory inside the parent directory
parentDir = fileparts(pwd); % Signal_sythesizer

% Define target results directory (../Signal_synthesizer/results)
resultsDir = fullfile(parentDir, "results");

% Guard clause: check if `results` exists.
if ~exist(resultsDir, "dir")
    mkdir(resultsDir);
end

% Save .mat file in results directory
save(fullfile(resultsDir, 'signalAnalysisResults.mat'));

% Create summary table
summaryTable = table();
summaryTable.Signal  = signalLabels';
summaryTable.RMS_LowNoise  = rms(lowDegMatrix, 2);
summaryTable.RMS_HighNoise  = rms(highDegMatrix, 2);
summaryTable.PeakToPeak_LowNoise = max(lowDegMatrix, [], 2) - min(lowDegMatrix, [], 2);
summaryTable.PeakToPeak_HighNoise = max(highDegMatrix, [], 2) - min(highDegMatrix, [], 2);
% Write table to CSV
writetable(summaryTable, fullfile(resultsDir, 'signalAnalysisSummary.csv'));
% Summary Table
disp(" ");
fprintf('--- Summary Table ---\n');
disp(summaryTable);
disp("**NOTE: RMS is shorthand for Root Mean Squared**");
fprintf('Analysis saved to signalAnalysisResults.mat (workspace variables) and signalAnalysisSummary.csv\n (summary)');
disp(' ');

%% Comparative Analysis Questions
disp(" ");
disp("---------Comparative Analysis Questions-------");
%1. How does noise affect signal quality (measured by RMS)?
fprintf(" - Noise misrepresents signal points, making it seem stronger on at random intervals," + ...
    " thus reducing its quality.\n\n");

%2. Which signal type (sine vs. square vs. triangle) is most affected by noise?
fprintf(" - From previous analysis, the square wave is less robust than the sine and triangles waves, so it is more affected by noise.\n\n");

%3. How effective is the moving average filter at reducing noise?

fprintf(" - The moving average filter is very effective at reducing noise.\n" + ...
    "\tWith signal quality improvement rates of about 8%% for robust wave-types like sine-wave\n" + ...
    "\tAnd up to 16%% for less robust waves like the Triangle wave.\n\n");
%4. Do higher frequencies show more degradation than lower frequencies?
fprintf(" - Yes. Higher frequencies show more degradation than lower frequencies\n" + ...
    "\tThis is shown by 10Hz having the highest degradation rate compared to 1 and 5 Hz.\n" + ...
    "\tIn another sense, the spontaneity of signal highs and lows in higher\n " + ...
    "\tfrequency signals makes its measured signals more prone to higher data misinterpretation\n\n");

%% LOCAL FUNCTIONS

function shouldTest(allSignals, TrueOrFalse)
% Determines if test case should run or not
% TrueOrFalse accepts bool: true | false

  if TrueOrFalse
    for i = ["Clean", "LowNoise", "HighNoise"]
        disp(i);
        cols = ["Sine1Hz", "Triangle10Hz"];
        for col = cols
            disp(col);  % MATLAB treats string arrays as concatenated arrays, so specify the string content in the cell
            disp(allSignals.(i).(col)(1:5));% first 5 samples
        end
        % Check for length of each Table
        fprintf("-------TABLE LENGTH CHECK for %s------", i);
        disp(" ");

        numRows = height(allSignals.(i));
        numCols = width(allSignals.(i));
        if numRows == 1000 && numCols == 10
            disp("PASS: Row length is 1000 and Number of Columns is 10");
        else
            fprintf("FAIL: Row length, %s AND Number of Columns is, %s", numRows, numCols);
        end
    end
  end
end

function degradation = percentDegradation(cleanVal, noisyVal)
% Percent Degradation Formula: ((clean_value - noisy_value) / clean_value) × 100
% Calculated the percentage degradation, and returns a value, array or
% matrix, depending on input type.
   degradation = abs(((cleanVal - noisyVal) ./ cleanVal) * 100);

end


%% VISUALIZATION (Sam)

%% FIGURE PATH SETUP

% Get script directory 
thisFile = mfilename('fullpath');
if isempty(thisFile)
    scriptDir = pwd;
else
    scriptDir = fileparts(thisFile);
end

% Navigate to project root (one level up from /scripts)
projectRoot = fileparts(scriptDir);

% Figures directory
figDir = fullfile(projectRoot, 'results', 'figures');

% Create if needed
if ~exist(figDir, 'dir')
    mkdir(figDir);
end

 
%%  FIGURE 1: Multi-Scenario Comparison (1 Hz Sine)

% Extract actual time-series signals from master table
clean1Hz     = allSignals.Clean.Sine1Hz;
lowNoise1Hz  = allSignals.LowNoise.Sine1Hz;
highNoise1Hz = allSignals.HighNoise.Sine1Hz;

figure;
subplot(3,1,1);
plot(clean1Hz, 'LineWidth', 1.4); 
title('Clean Signal (1 Hz Sine Wave)');
xlabel('Samples'); ylabel('Amplitude');
ylim([-2 2]); grid on;

subplot(3,1,2);
plot(lowNoise1Hz, 'LineWidth', 1.4); 
title('Low Noise Signal (1 Hz Sine Wave)');
xlabel('Samples'); ylabel('Amplitude');
ylim([-2 2]); grid on;

subplot(3,1,3);
plot(highNoise1Hz, 'LineWidth', 1.4); 
title('High Noise Signal (1 Hz Sine Wave)');
xlabel('Samples'); ylabel('Amplitude');
ylim([-2 2]); grid on;

saveas(gcf, fullfile(figDir, 'scenario_comparison.png'));

%%  FIGURE 2: Frequency Comparison (Clean Sine Waves)

% Build clean sine waves for 0–2 seconds
t = 0:0.001:2;   % 1000 Hz sampling for visualization
sig1  = sin(2*pi*1*t);
sig5  = sin(2*pi*5*t);
sig10 = sin(2*pi*10*t);

figure;
hold on;
plot(t, sig1, 'LineWidth', 1.5, 'DisplayName', '1 Hz');
plot(t, sig5, '--', 'LineWidth', 1.5, 'DisplayName', '5 Hz');
plot(t, sig10, ':', 'LineWidth', 1.5, 'DisplayName', '10 Hz');
hold off;

title('Frequency Comparison of Clean Sine Waves');
xlabel('Time (s)'); ylabel('Amplitude');
xlim([0 2]);
grid on; legend show;

saveas(gcf, fullfile(figDir, 'frequency_comparison.png'));

%%  FIGURE 3: Signal Quality Metrics (RMS Only)

% RMS degradation matrices already computed:
% lowDegMatrix  (3×3) rows = Sine/Square/Triangle, cols = 1/5/10 Hz
% highDegMatrix "

rmsLow  = rms(lowDegMatrix,  2);
rmsHigh = rms(highDegMatrix, 2);

figure;
bar([rmsLow rmsHigh]);
set(gca, 'XTickLabel', {'Sine', 'Square', 'Triangle'});
ylabel('RMS Degradation (%)');
legend({'Low Noise', 'High Noise'}, 'Location', 'northwest');
title('RMS-Based Signal Quality Comparison');
grid on;

saveas(gcf, fullfile(figDir, 'quality_metrics.png'));

%%  FIGURE 4: Before/After Filtering (High Noise Sine)

% Use full noisy time-series (not degradation matrix)
highNoiseSine = allSignals.HighNoise.Sine1Hz;
filteredSine  = movmean(highNoiseSine, 5);

figure;
subplot(2,1,1);
plot(highNoiseSine, 'LineWidth', 1.3);
title('High Noise 1 Hz Sine (Before Filtering)');
xlabel('Samples'); ylabel('Amplitude'); grid on;

subplot(2,1,2);
plot(filteredSine, 'LineWidth', 1.4); hold on;
plot(highNoiseSine, '--', 'LineWidth', 1);
hold off;
title('After Moving Average Filtering');
xlabel('Samples'); ylabel('Amplitude');
legend('Filtered', 'Original'); grid on;

saveas(gcf, fullfile(figDir, 'filtering_effect.png'));

%%  FIGURE 5: Signal Dashboard (4 Panels)

figure('Position',[100 100 1200 800]);

% Panel 1: Degradation by Frequency
subplot(2,2,1);
bar([mean(lowDegMatrix); mean(highDegMatrix)]');
title('Avg % Degradation per Frequency');
xlabel('Frequency (Hz)'); ylabel('Degradation (%)');
set(gca,'XTickLabel',{'1','5','10'});
legend({'Low Noise','High Noise'}); grid on;

% Panel 2: Zero Crossings (1 Hz Sine)
zcLow  = sum(diff(sign(clean1Hz))~=0);
zcHigh = sum(diff(sign(highNoise1Hz))~=0);
subplot(2,2,2);
bar([zcLow zcHigh]);
title('Zero Crossings — 1 Hz Sine');
set(gca, 'XTickLabel', {'Clean','High Noise'}); ylabel('Count');
grid on;

% Panel 3: Peak-to-Peak
ptpLow  = peak2peak(clean1Hz);
ptpHigh = peak2peak(highNoise1Hz);
subplot(2,2,3);
bar([ptpLow ptpHigh]);
title('Peak-to-Peak Amplitude — 1 Hz Sine');
set(gca,'XTickLabel',{'Clean','High Noise'});
ylabel('Amplitude'); grid on;

% Panel 4: RMS Comparison (Raw Time Series)
subplot(2,2,4);
rmsComp = [rms(highNoiseSine), rms(filteredSine), rms(clean1Hz)];
bar(rmsComp);
set(gca,'XTickLabel',{'Noisy','Filtered','Original'});
ylabel('RMS'); title('Filtering Effectiveness'); grid on;

saveas(gcf, fullfile(figDir, 'signal_dashboard.png'));


%% EXTENSION
%% Signal Degradation Analysis: Sine, Square, Triangle Waves
% Author: William Dakare
% Purpose: Compare robustness of different waveform types under noise
% Date: 2025-11-23

%% SET UP
% Define time vector and frequency
fs = 1000;             % Sampling frequency in Hz
t = 0:1/fs:1;          % 1 second duration
freq = 5;              % Signal frequency in Hz

% Generate clean signals
cleanSine     = sin(2*pi*freq*t);          % Sine wave
cleanSquare   = square(2*pi*freq*t);      % Square wave
cleanTriangle = sawtooth(2*pi*freq*t,0.5);% Triangle wave (symmetric)

% Generate noisy signals (add Gaussian noise)
noiseLevel = 0.3; % Standard deviation of noise
noisySine     = cleanSine + noiseLevel*randn(size(t));
noisySquare   = cleanSquare + noiseLevel*randn(size(t));
noisyTriangle = cleanTriangle + noiseLevel*randn(size(t));

%% STEP 1: Calculate Degradation Percentage
% Formula: ((std_noisy - std_clean) / std_clean) * 100

% Degradation for Sine
degradSine = (std(noisySine) - std(cleanSine)) / std(cleanSine) * 100;

% Degradation for Square
degradSquare = (std(noisySquare) - std(cleanSquare)) / std(cleanSquare) * 100;

% Degradation for Triangle
degradTriangle = (std(noisyTriangle) - std(cleanTriangle)) / std(cleanTriangle) * 100;

%% STEP 2: Compare Degradation Across Waveforms
% Combine degradation values into array for comparison
combinedDegrad = [degradSine, degradSquare, degradTriangle];
waveLabels = ["Sine","Square","Triangle"];

% Find the waveform with minimum degradation (most robust)
% This index maps directly to waveLabels
minDegrad = min(combinedDegrad);
robustWaveIdx = find(minDegrad == combinedDegrad); 
robustWave = waveLabels(robustWaveIdx);

fprintf('Most robust waveform under noise: %s (%.2f%% degradation)\n', robustWave, minDegrad);

%% STEP 3: Plot Side-by-Side Signals
figure;
subplot(3,1,1);
plot(t,noisySine,'r'); hold on;
plot(t,cleanSine,'k--'); % Overlay clean signal for reference
title('Sine Wave (Noisy vs Clean)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Noisy','Clean');

subplot(3,1,2);
plot(t,noisySquare,'b'); hold on;
plot(t,cleanSquare,'k--');
title('Square Wave (Noisy vs Clean)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Noisy','Clean');

subplot(3,1,3);
plot(t,noisyTriangle,'g'); hold on;
plot(t,cleanTriangle,'k--');
title('Triangle Wave (Noisy vs Clean)');
xlabel('Time (s)'); ylabel('Amplitude');
legend('Noisy','Clean');

%% STEP 4: Overlaid Plot for Direct Comparison
figure;
plot(t,noisySine,'r','DisplayName','Sine'); hold on;
plot(t,noisySquare,'b','DisplayName','Square');
plot(t,noisyTriangle,'g','DisplayName','Triangle');
xlabel('Time (s)'); ylabel('Amplitude');
title('Noisy Signals Comparison (Overlaid)');
legend;

%% STEP 5: Summary Table of Degradation
degradTable = table(waveLabels', combinedDegrad', 'VariableNames', {'Waveform','DegradationPercent'});
disp('Degradation Summary Table:');
disp(degradTable);

%% STEP 6: Interpretation
% Sine waves tend to maintain smooth curves and are less affected by noise
% Square waves have sharp transitions and are most sensitive to noise spikes
% Triangle waves have linear slopes, hence they are moderate robust, but
% still suceptible to signal disruption
