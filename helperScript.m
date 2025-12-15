% This will be the helperScript with all our data creation and processing functions. To make sure all
% the functions are accesible in the current directory, make all the
% functions as methods in a class 'helperScipt' matching the file name
% Function-methods contained: (Call all functions in `signalProcessingAlgorithm.m`)
%  - processSignals(processType)
%  - testProcessor(TrueOrFalse) Tests `processSignals()`
% to run refer to class `helperScript` first
% e.g, helperScript.processSignals(processType)

classdef helperScript
   methods(Static)
%% Create Data 
        
     function projectCsvMaker
        % Scenario File Creator
        %
        % Author: Noah Nir
        % Date: November 2025
        %
        % Description:
        %   This script will create CSV's for different noise distributions
        %   
        % Outputs: CSV File
        
        %% Scenario A: Setup
        duration = 10;
        samplingRate = 100;
        Time = linspace(0, duration, duration * samplingRate);
        A = 1;
        
        parentPath = fileparts(pwd); % Get parent path, this is where we want the directory for our outputs to be
        targetPath = fullfile(parentPath, 'data', 'signals'); % get into the target path "signals"
        
        if ~exist(targetPath, 'dir') % don't make new ones on reruns
            mkdir(targetPath);
        end
        
        
        % We will generate 3 wave-types(sine, square and triangle) at 3 different
        % frequencies: 1, 5 & 10Hz
        %% 1 Hz Waves
        f1 = 1;
        Sine1Hz      = A * sin(2*pi*f1*Time);
        Square1Hz    = sign(sin(2*pi*f1*Time));
        Triangle1Hz  = (2/pi) * asin(sin(2*pi*f1*Time));
        
        %% 5 Hz Waves
        f5 = 5;
        Sine5Hz      = A * sin(2*pi*f5*Time);
        Square5Hz    = sign(sin(2*pi*f5*Time));
        Triangle5Hz  = (2/pi) * asin(sin(2*pi*f5*Time));
        
        %% 10 Hz Waves
        f10 = 10;
        Sine10Hz     = A * sin(2*pi*f10*Time);
        Square10Hz   = sign(sin(2*pi*f10*Time));
        Triangle10Hz = (2/pi) * asin(sin(2*pi*f10*Time));
        
        %% Create Table and Write CSV (Scenario A Clean Signals)
        T = table( ...
            Time.', ...
            Sine1Hz.', Sine5Hz.', Sine10Hz.', ...
            Square1Hz.', Square5Hz.', Square10Hz.', ...
            Triangle1Hz.', Triangle5Hz.', Triangle10Hz.', ...
            'VariableNames', { ...
            'Time', ...
            'Sine1Hz', 'Sine5Hz', 'Sine10Hz', ...
            'Square1Hz', 'Square5Hz', 'Square10Hz', ...
            'Triangle1Hz', 'Triangle5Hz', 'Triangle10Hz' ...
            });
        
        writetable(T, fullfile(targetPath, "Clean/",'Clean.csv'));
        disp("CSV file 'Clean.csv' created successfully!");
        
        %% Scenario B: Setup (Low Noise)
        noiseLevel = 0.1;
        
        % Add 10% of a random value, to each signal and frequency-variant as noise to simulate low noise level
        %% Low Noise 1 Hz
        Sine1Hz     = Sine1Hz     + noiseLevel * randn(size(Time));
        Square1Hz   = Square1Hz   + noiseLevel * randn(size(Time));
        Triangle1Hz = Triangle1Hz + noiseLevel * randn(size(Time));
        
        %% Low Noise 5 Hz
        Sine5Hz     = Sine5Hz     + noiseLevel * randn(size(Time));
        Square5Hz   = Square5Hz   + noiseLevel * randn(size(Time));
        Triangle5Hz = Triangle5Hz + noiseLevel * randn(size(Time));
        
        %% Low Noise 10 Hz
        Sine10Hz     = Sine10Hz     + noiseLevel * randn(size(Time));
        Square10Hz   = Square10Hz   + noiseLevel * randn(size(Time));
        Triangle10Hz = Triangle10Hz + noiseLevel * randn(size(Time));
        
        %% Create Table (Low Noise) and Write CSV
        T = table( ...
            Time.', ...
            Sine1Hz.', Sine5Hz.', Sine10Hz.', ...
            Square1Hz.', Square5Hz.', Square10Hz.', ...
            Triangle1Hz.', Triangle5Hz.', Triangle10Hz.', ...
            'VariableNames', { ...
            'Time', ...
            'Sine1Hz', 'Sine5Hz', 'Sine10Hz', ...
            'Square1Hz', 'Square5Hz', 'Square10Hz', ...
            'Triangle1Hz', 'Triangle5Hz', 'Triangle10Hz' ...
            });
        
        writetable(T, fullfile(targetPath, "LowNoise/",'lowNoise.csv'));
        disp("CSV file 'lowNoise.csv' created successfully!");
        
        %% Scenario C: Setup (Large Noise)
        noiseLevel = 0.5;
        
        % Add 50% of a random value, to each signal and frequency-variant as noise to simulate low noise level
        
        %% High Noise 1 Hz
        Sine1Hz     = Sine1Hz     + noiseLevel * randn(size(Time));
        Square1Hz   = Square1Hz   + noiseLevel * randn(size(Time));
        Triangle1Hz = Triangle1Hz + noiseLevel * randn(size(Time));
        
        %% High Noise 5 Hz
        Sine5Hz     = Sine5Hz     + noiseLevel * randn(size(Time));
        Square5Hz   = Square5Hz   + noiseLevel * randn(size(Time));
        Triangle5Hz = Triangle5Hz + noiseLevel * randn(size(Time));
        
        %% High Noise 10 Hz
        Sine10Hz     = Sine10Hz     + noiseLevel * randn(size(Time));
        Square10Hz   = Square10Hz   + noiseLevel * randn(size(Time));
        Triangle10Hz = Triangle10Hz + noiseLevel * randn(size(Time));
        
        %% Create Table (High Noise) and Write CSV
        T = table( ...
            Time.', ...
            Sine1Hz.', Sine5Hz.', Sine10Hz.', ...
            Square1Hz.', Square5Hz.', Square10Hz.', ...
            Triangle1Hz.', Triangle5Hz.', Triangle10Hz.', ...
            'VariableNames', { ...
            'Time', ...
            'Sine1Hz', 'Sine5Hz', 'Sine10Hz', ...
            'Square1Hz', 'Square5Hz', 'Square10Hz', ...
            'Triangle1Hz', 'Triangle5Hz', 'Triangle10Hz' ...
            });
        
        writetable(T, fullfile(targetPath,"HighNoise/" ,'largeNoise.csv'));
        disp("CSV file 'largeNoise.csv' created successfully!");
     end

%% Process Data
    function processedData = processSignals(processType)
    % processSignals Compute signal metrics modularly for multiple signal types
    % Inputs:
    %   processType - 'RMS', 'peakToPeak', or 'zeroCross'
    %   relativePathBase - relative path to your signal folders (e.g., '.\Signals')
    %
    % Output:
    %   processedData - struct containing processed metrics
    %
    % Requirements:
    % A. RMS (Root Mean Square) Value:
    %   • Measures the "average power" of the signal
    %   • Formula: sqrt(mean(signal.^2))
    %   • Calculate for all signals (3 frequencies × 3 wave types × 3 scenarios = 27 combinations)
    %
    % B. Peak-to-Peak Amplitude:
    %   • Difference between maximum and minimum values: max(signal) - min(signal)
    %   • For clean sine waves, expect approximately 2.0
    %   • Noise will increase this value
    %
    % C. Zero Crossings:
    %   • Count how many times signal crosses zero (changes from positive to negative or vice versa)
    %   • Use diff(sign(signal)) to detect sign changes
    %   • Count non-zero elements in result
    %   • For 1 Hz signal over 10 seconds, expect ~20 crossings
        
        % Guard clause against invalid input 
        validInputs = ["rms", "peaktopeak", "zerocross"];
        
        if ~ismember(processType, validInputs)
            error("Invalid input. Valid inputs: ''rms'', ''peaktopeak'', ''zerocross''");
        end

        % Make a dictionary of Dir : File pairs i.e. "Clean" : "Clean.csv"
        allSignalsDirFilePair = dictionary(["Clean", "LowNoise", "HighNoise"], ...
                                           ["Clean.csv", "lowNoise.csv", "largeNoise.csv"]);
    
        allKeys = keys(allSignalsDirFilePair);
    
        % Initialize master output struct
        processedData = struct();
        % Create a subfield for the type of processing
        processedData.(processType) = struct();
    
        % Loop over each signal scenario/folder
        for key = 1:length(allKeys) 
            dirName = allKeys{key};         % e.g., "Clean"
            fileName = allSignalsDirFilePair(dirName); % corresponding CSV file name, e.g., "Clean.csv"
            
            relativePathBase = "..\data\signals"; % Add these at the end for appropriate dir: (Clean|LowNoise|HighNoise)
            filePath = fullfile(relativePathBase, dirName, fileName);
    
            % Read CSV into table
            signalTable = readtable(filePath);
    
            % Make column names valid for struct fields (CORE DEPENDENCY for this structure approach)
            signalTable.Properties.VariableNames = matlab.lang.makeValidName(signalTable.Properties.VariableNames);
    
            signalStruct = struct();  % temporary struct for this folder/scenario
            
            % Loop through columnns and process based on `processType`
            for col = 2:width(signalTable) % Ignore time column
                colName = signalTable.Properties.VariableNames{col};  % column header
                signalVec = signalTable.(colName);                    % column vector
    
                % Perform the requested metric calculation
                switch lower(processType) % Make case uniform (processType is case insensitive)
                    case 'rms'
                        % RMS (Root Mean Square) value
                        signalStruct.(colName) = sqrt(mean(signalVec.^2));
                    case 'peaktopeak'
                        % Peak-to-Peak amplitude
                        signalStruct.(colName) = max(signalVec) - min(signalVec);
                    case 'zerocross'
                        % Zero-crossings
                        % abs() because sign-function returns positive and negative values when signs change
                        % Account for both
                        signalStruct.(colName) = sum(abs(diff(sign(signalVec))) > 0); 
                    otherwise
                        % error handling: unexpected `processType`
                        error('Unknown processType: %s', processType);
                end
    
            end
    
            % Store the processed struct under this scenario
            processedData.(processType).(dirName) = signalStruct;
    
        end
    end

    %% ------------------ TEST for Data Validity ------------------
    % Call in signalProcessingAlgorithm.m to test
    
    function testProcessor(TrueOrFalse)
       if TrueOrFalse
        for i = ["Clean", "LowNoise", "HighNoise"]
          disp(i);
          processTypeArr = ["rms", "peaktopeak", "zerocross"];
          
          % Test all cases that are processed in `processSignals`
          for processType = processTypeArr
        
              processedData = helperScript.processSignals(processType);
              data = processedData;
    
              if isfield(data, processType)
                  field = processedData.(processType).(i);
                  
                  % Check if attribute contains 'Sine10Hz'
                  if isfield(field, 'Sine10Hz')
                      data = field.Sine10Hz;
                      disp(data(1:min(5,end))); % show first 5 samples or fewer
                  else
                      fprintf('  Sine10Hz not present for process "%s". Available fields:\n', processType);
                      disp(fieldnames(field));
                  end
              else
                  fprintf('  No data for process "%s".\n', processType);
              end
          end
        end
       end
    end
  end
end
