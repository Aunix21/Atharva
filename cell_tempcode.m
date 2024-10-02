scalingfactor = 2/3;
% Load the data from the Excel file
data = readmatrix('time_current_onelap.xlsx');  % Excel filename with current vs time data

% Extract and scale current profile
currentProfile_singleLap = data(:, 2) * scalingfactor;

% Number of laps and total data points
numLaps = 18; %total number of laps
dataPointsPerLap = length(currentProfile_singleLap);
numDataPoints = numLaps * dataPointsPerLap;

% Preallocate temperature array and initialize
temperatureProfile = zeros(numDataPoints, 1);
temperatureProfile(1) = 25;  % Initial temperature
heat_transfer_coeff=50;      % in W/m2
Area_of_ht_cell=0.0092;      % in m2
Temp_air_initial=25;
mass_of_cell=0.07;           %in kg
specific_heat_cell=840;      %Cp in J/Kkg
k1=1/(mass_of_cell*specific_heat_cell);
R_cell=0.007;                %in ohm

% Precompute constants
timeStep = 0.0001;
currentFactor = R_cell * k1 * timeStep;
temperatureFactor = -k1 * heat_transfer_coeff * Area_of_ht_cell * timeStep;
constantFactor =  heat_transfer_coeff * Area_of_ht_cell * Temp_air_initial * k1 * timeStep;

% Expand the current profile for 18 laps
currentProfile = repmat(currentProfile_singleLap, numLaps, 1);

% Calculate temperature profile 
for index = 2:numDataPoints
    previousTemperature = temperatureProfile(index - 1);
    current = currentProfile(index);
    
    % Calculate new temperature
    newTemperature = currentFactor * current^2 ...
                     + temperatureFactor * previousTemperature ...
                     + previousTemperature + constantFactor;
    
    % Assign the calculated temperature to the output
    temperatureProfile(index) = newTemperature;
    
    % Check if the current index is the end of a lap
    if mod(index, dataPointsPerLap) == 0
        lapNumber = index / dataPointsPerLap;
        fprintf('Temperature at the end of Lap %d: %.2f°C\n', lapNumber, temperatureProfile(index));
    end
end

%% SAVING DATA INTO EXCEL FILES

%Save the current profile to a single Excel file
fileNameCurrent = 'time_vs_currentProfile.xlsx';
writematrix(currentProfile_singleLap, fileNameCurrent);

% Save data for each lap's temperature profile to separate Excel files
for lap = 1:numLaps
    startIndex = (lap - 1) * dataPointsPerLap + 1;
    endIndex = lap * dataPointsPerLap;
    
    % Extract the temperature data for the current lap
    temperatureProfileLap = temperatureProfile(startIndex:endIndex);
    
    % Adjust timeProfileLap to start from 0
    timeProfileLap = ((startIndex:endIndex) - startIndex)' * timeStep;
    
    % Save the lap data to Excel
    fileNameTemperature = sprintf('temperatureProfile_Lap%d.xlsx', lap);
    
    % Write the temperature data to Excel file
    writematrix([timeProfileLap, temperatureProfileLap], fileNameTemperature);
end

%% PLOTTING GRAPHS

% Plot the current profile over time
timeProfile = (0:numDataPoints-1)' * timeStep;
tiledlayout(2,1);
nexttile
plot(timeProfileLap, currentProfile_singleLap);
xlabel('Time');
ylabel('Current');
title('Current vs Time (1 Lap)');

% Plot the temperature over time
nexttile
plot(timeProfile, temperatureProfile);
xlabel('Time');
ylabel('Temperature');
title('Temperature vs Time (18 Laps)');