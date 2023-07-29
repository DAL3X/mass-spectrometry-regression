%% Default constants
% Constants needed to find and read in the data file 
firstDataLine = 379;                        % First line in the file that contains measurements
inputFiletype = "*.dat";                    % File type of the data file
timeAxisIdentifier = "Time Relative (sec)"; % The name of the variable that should be used as the time axis (for regression and section division)
valueAxisIdentifier = "44_amu";             % Having this word in the variable name identifies the value variables on which the regression is performed
temperaturAxisIdentifier = "Temperature";   % The name of the variable which contains the temperature information

% Constants needed to divide the data into regression sections.
gradientDistance = 350;        % Distance (in number of data points) in which gradients get calculated
roomTemperatureChange = 1.3;   % Threshold (in the unit of temperatureAxis) for detection of changes from, or to room temperature
heatTemperatureChange = 1.8;   % Threshold (in the unit of temperatureAxis) for detection of changes from, or to maximum heating temperature


%% Select and import the data file
[file,path] = uigetfile(inputFiletype);
try
    opts = detectImportOptions(strcat(path, file),'ReadVariableNames', true, 'NumHeaderLines', firstDataLine-2, 'DecimalSeparator', ',', 'VariableNamingRule', 'preserve');
    data = readtable(strcat(path, file), opts);
    clear opts;
catch
    disp('Could not load data. Shutting down!');
    return % Ends the script if the data could not be loaded
end


%% Identify all columns of interest in the data
timeAxis = data(:,timeAxisIdentifier);
temperatureAxis = data(:,temperaturAxisIdentifier);
valueAxesNames = []; 
% Find all column names in the data which contain the valueAxisIdentifier
for name = data.Properties.VariableNames
    if contains(name, valueAxisIdentifier)
        valueAxesNames = cat(2, valueAxesNames, name);
    end
end
if isempty(valueAxesNames)
    disp('No value axes found. Please check if the value axis identifier is set correctly. Shutting down!');
    return % Ends the script if no values were found.
end

%% Find heating intervals and divide the data accordingly
% Get the max temperature (max heating) and min temperature (room temperature)
maxTemp = max(table2array(temperatureAxis));
minTemp = min(table2array(temperatureAxis));

% Find the heating intervals
borders = getSectionsBorders(table2array(timeAxis), table2array(temperatureAxis), minTemp, maxTemp, roomTemperatureChange, heatTemperatureChange, gradientDistance);

% Divide the whole data table in section tables and save them in a cell array
tables = cell(length(borders)+1,1);
tables{1,:} = data(1:borders(1)-1,:); % First section is determined outside of loop
for i = 1:(length(borders))
    if i ~= length(borders)
        tables{i+1,:} = data(borders(i):borders(i+1)-1,:);
    else
        tables{i+1,:} = data(borders(i):end,:);
    end
end


%% Show all graphs
figure();
dat = data(:,valueAxesNames(1));
plot(table2array(timeAxis),table2array(dat),'.');
set(gca,'FontSize',30);
xlabel("Zeit in Sekunden (x10^5)");
ylabel("Ausschlag im Massenspektrum");

figure();
plot(table2array(timeAxis), table2array(temperatureAxis));
set(gca,'FontSize',30)
xlabel("Zeit in Sekunden (x10^5)");
ylabel("Temperatur in °C");

figure();
plot(table2array(timeAxis), table2array(temperatureAxis));
set(gca,'FontSize',30);
xlabel("Zeit in Sekunden (x10^5)");
ylabel("Temperatur in °C");
lines = table2array(timeAxis);
xline(lines(borders),'r');

%% Finds the heating sections and returns the border indicies between them
% time and temp are the time and temperature vectors.
% gradDist determines the gradient interval width.
% state determines the starting state (0, 1, 2, or 3).
% The thresholds determine the thresholds for changes in the gradients.

function secBorders = getSectionsBorders(time, temp, roomTemp, heatTemp, roomTemperatureChange, heatTemperatureChange, gradDist)
    derivative = deriv(time, temp, gradDist);
    index = 1;
    secBorders = [];
    searchState = 0;
    % This while loop searches for borders between the states 
    % (0) room temperature,
    % (1) actively increasing temperature,
    % (2) holding increased temperature,
    % (3) decreasing the temperature
    while index < (length(temp) - gradDist)
        if (searchState == 0 && temp(index) >= (roomTemp + roomTemperatureChange) && derivative(index) > 0) ... 
        || (searchState == 1 && temp(index) >= (heatTemp - heatTemperatureChange)) ...
        || (searchState == 2 && temp(index) < (heatTemp - heatTemperatureChange) && derivative(index) < 0) ...
        || (searchState == 3 && temp(index) < (roomTemp + roomTemperatureChange))
        % Found change (rising heat -> constant heat) in state 1 or
        % (constant heat -> decreasing heat) in state 2 or
        % (room temperature -> rising heat) in state 0 or 
        % (decreasing heat -> room temperature) in state 3
            secBorders = cat(1, secBorders, index);
            searchState = mod(searchState + 1, 4);
            index = index + 20;
        end
        index = index + 1;
    end
end


%% Calculates the discrete-time derivative averaged over a certain amount of data points in order to deal with errors in measurements.
% The vectors are both required to have a dimension of 1xN with the same N.
% The derivative is stepwise calculated.
% In each step the discrete-time derivatives for the next width data points get calculated and averaged. 
% Returns a 1xN vector containing all the (N - width) gradients followed by width times the last calculated gradient.

function deriv = deriv(time,value,width)
    % Check if the time and value vectors are column vectors with the same size
    if size(time,1) ~= size(value,1) || size(time,2) ~= 1 || size(value,2) ~= 1
        error('Error while calculating the discrete-time derivative: Incompatible dimensions.')
    end
    length = size(time,1);
    deriv = zeros(length, 1);
    for index = 1:length
        if index > (length-width)
            % Fill the last interval entries with the last calculated gradient
            deriv(index,1) = deriv(index-1,1);
        else
            % Calculate the gradients and average them
            gradientSum = 0;
            for subIndex = 1:width    
                gradientSum = gradientSum + ((value(index + subIndex,1) - value(index,1)) / (time(index + subIndex,1) - time(index,1)));
            end
            deriv(index,1)= gradientSum / width;
        end
    end
end