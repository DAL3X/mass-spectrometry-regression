%% Default constants
% Constants needed to find and read in the data file 
firstDataLine = 379;                        % First line in the file that contains measurements
inputFiletype = "*.dat";                    % File type of the data file
timeAxisIdentifier = "Time Relative (sec)"; % The name of the variable that should be used as the time axis (for regression and section division)
valueAxisIdentifier = "amu";             % Having this word in the variable name identifies the value variables on which the regression is performed
temperaturAxisIdentifier = "Temperature";   % The name of the variable which contains the temperature information

% Constants needed to divide the data into regression sections.
gradientDistance = 350;        % Distance (in number of data points) in which gradients get calculated
roomTemperatureChange = 1.3;   % Threshold (in the unit of temperatureAxis) for detection of changes from, or to room temperature
heatTemperatureChange = 1.8;   % Threshold (in the unit of temperatureAxis) for detection of changes from, or to maximum heating temperature

% Constants for optional visualizations and formats
showRegressionGraphs = true;    % Determines if the script outputs all the regression data as graphs
showTemperatureGraph = false;     % Determines if the script outputs the temperature data as a graph (set on true to verify heating sections)
consoleOutput = false;           % Determines if the regression results get displayed in the command window
fileOutput = true;              % Determines if the regression results get saved in files
outputFormat = ".txt";           % File type of the output regression data
selectionWindowSize = [350,250]; % Determines the function specify GUI window size

% Constants needed for the regression
weighted = false;            % Determines if the regression algorithm uses weights
weightFuction = 'bisquare';  % If weighted = true, this determines the used weight function
                             % Possible functions are : 'andrews', 'bisquare', 'cauchy', 'fair', 'huber','logistic', 'talwar', and 'welsch'
regressionWidth = 0.8;  % Determines how much of a section is used for the regression step
                        % The values get cut on both ends equally (e.g. using 0.8 cuts 10% of data from each end) 

% All supported functions
constant = {(@(b,x) b(1)+0.*x), 1};  % a
linear = {(@(b,x) b(1)+b(2).*x), 2};    % a+bx
quadratic = {(@(b,x) b(1)+b(2).*x+b(3).*(x.^2)), 3};    % a+bx+cx^2
cubic = {(@(b,x) b(1)+b(2).*x+b(3).*(x.^2)+b(4).*(x.^3)), 4};    % a+bx+cx^2+dx^3 
exponential = {(@(b,x) b(1)+b(2).*exp(x.*b(3))), 3};   % a+b*exp(c*x)
logarithm = {(@(b,x) b(1)+b(2).*reallog(x+abs(b(3)))), 3};    % a+b*log(x+c) with c >= 0
root = {(@(b,x) b(1)+b(2).*realsqrt(x+abs(b(3)))), 3};    % a+b*root(x+c) with c >= 0
logistic = {(@(b,x) (b(1)./(1+b(2).*(exp(-b(3).*(x)))))), 3}; % a/(1+b*e^(-c*x))
functions = {constant, linear, quadratic, cubic, exponential, logarithm, root, logistic};


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


%% Get the user selection for specific functions on specific intervals
[valIndex, valTf] = listdlg('PromptString','Choose (multiple) value axes:','ListString',valueAxesNames,'ListSize',selectionWindowSize);
selectionComplete = valTf;
funcSpecify = [];
for valNum = valIndex
    valName = valueAxesNames(valNum);
    secStringList = cell(length(tables),1);
    % Find the section borders
    for i = 1:length(tables)
        tabldata = table2array(tables{i}(:,timeAxisIdentifier));
        secstr = "Section " + i + " (" + tabldata(1) + " - " + tabldata(end) + ")";
        secStringList{i} = secstr;
    end
    secString = "Choose sections for " + valName;
    [secIndex, secTf] = listdlg('PromptString',secString,'ListString',secStringList,'ListSize',selectionWindowSize);
    selectionComplete = selectionComplete & secTf;
    for secNum = secIndex
        funcStringList = cell(length(functions),1);
        funcString = "Choose function for " + valName + " section " + secNum;
        for funcI = 1:length(functions)
            funcStringList{funcI} = func2str(functions{funcI}{1});
        end
        [funcIndex, funcTf] = listdlg('PromptString',funcString,'ListString',funcStringList,'ListSize',selectionWindowSize);
        selectionComplete = selectionComplete & funcTf;
        funcSpecify = cat(2, funcSpecify,  [valNum,secNum,funcIndex]);
    end
end
if selectionComplete == 0
    % No selections were made or something went wrong during the selection
    funcSpecify = [];
    disp("Insufficient selection. The script will try to determine all functions automatically!");
end


%% Use regression on all of the section for all columns of interest
warning('off'); % Disabled to prevent console spamming
outputData = cell(length(valueAxesNames),1);
funcSpecifyCounter = 1;
for valueNr = 1:length(valueAxesNames)
    valueName = valueAxesNames(valueNr);
    outputSectionData = cell(length(tables),1);
    for section = 1:length(tables)
        tbl = tables{section};
        X = table2array(tbl(:,timeAxisIdentifier));
        Y = table2array(tbl(:,valueName));
        % Cut the values according to the regressionWidth
        cutValue = (length(X)*(1-regressionWidth))/2;
        cutX = X(cutValue+1:length(X)-cutValue);
        cutY = Y(cutValue+1:length(Y)-cutValue);
        % Set the beginning of the time axis to 0 for every interval
        cutXAdjusted = arrayfun(@(x) (x-cutX(1))+realmin, cutX);
        maxRValue = intmin;
        bestModel = NonLinearModel;
        % Check only for specific functions, if specified
        funcIndices = 1:length(functions);
        if funcSpecifyCounter < length(funcSpecify) && funcSpecify(funcSpecifyCounter) == valueNr && funcSpecify(funcSpecifyCounter + 1) == section
            funcIndices = funcSpecify(funcSpecifyCounter + 2);
            funcSpecifyCounter = funcSpecifyCounter + 3;
        end
        for funcIndex = funcIndices
            % At this point we have a section of a data and time column and a function with parameters
            % This is where the actual regression happens
            func = functions{funcIndex};
            modelfun = func{1};
            beta0 = zeros(func{2},1); % DO NOT CHANGE to anything other than zero to avoid overflows
            if weighted
                opts = statset('fitnlm');
                opts.Robust = 'on';
                opts.RobustWgtFun = weightFuction;
                mdl = fitnlm(cutXAdjusted,cutY,modelfun,beta0,'Options',opts);
            else
                mdl = fitnlm(cutXAdjusted,cutY,modelfun,beta0);
            end
            if (mdl.Rsquared.Ordinary > maxRValue) || (maxRValue == intmin)
                bestModel = mdl;
                maxRValue = mdl.Rsquared.Adjusted;
            end
        end
        outputSectionData{section} = {X,Y,cutX,bestModel};
    end
    outputData{valueNr} = {valueName,outputSectionData};
end
warning('on'); % No more console spamming. Re-enable warnings
output(outputData, timeAxisIdentifier, timeAxis, temperatureAxis, borders, showRegressionGraphs, ...
    showTemperatureGraph, consoleOutput, fileOutput, outputFormat, weighted, weightFuction, regressionWidth);


%% Output all graphs and regression data according to settings
% Tons of messy output stuff. Should probably stay untouched

function output(regData, timeAxisIdentifier, timeAxis, tempAxis, borders, showRegressionGraphs, ...
showTemperatureGraph, consoleOutput, fileOutput, outputFormat, weighted, weightFuction, regressionWidth)
    if showTemperatureGraph
        % Plot the temperature graph to visually verify the section borders
        figure();
        plot(table2array(timeAxis), table2array(tempAxis));
        lines = table2array(timeAxis);
        xline(lines(borders));
    end
    if fileOutput
        dirPath = uigetdir;
    end
    for valueNr = 1:length(regData)
        valueName = regData{valueNr}{1};
        if fileOutput
            filepath = string(dirPath) + "/" + valueName + outputFormat;
            if isfile(filepath)
                % Delete if file present
                delete(filepath);
            end
            fid = fopen(filepath, 'a');
            fprintf(fid, "Weighted: " + weighted + "\n");
            if weighted
                fprintf(fid, "Weightfunction: " + weightFuction + "\n");
            end
            fprintf(fid, "Cut time interval: " + regressionWidth + "\n");
            fprintf(fid, "\n");
        end
        if showRegressionGraphs
            % Create a new regression graph for this data column if enabled
            figure();
        end
        if consoleOutput
            disp("Regression data for: " + valueName);
        end
        for section = 1:length(regData{valueNr}{2})
            X = regData{valueNr}{2}{section}{1};
            Y = regData{valueNr}{2}{section}{2};
            cutX = regData{valueNr}{2}{section}{3};
            bestModel = regData{valueNr}{2}{section}{4};
            timeDescription = X(1) + " - " + X(end) + " (adjusted to " + cutX(1) + " - " + cutX(end) + ")";
            coefficients = string(table2cell(bestModel.Coefficients(:,'Estimate')));
            for i = 1:length(coefficients)
                if contains(bestModel.Formula.Expression, "x + abs(b" + i)
                    % If the coefficient was abs, change to positiv number
                    coefficients(i) = abs(double(coefficients(i)));
                end
                desc = "b" + i;
                coefficients(i) = desc + ":" + coefficients(i);
            end
            coefficients = join(coefficients,', ');
            if showRegressionGraphs
                % Add the section to the regression graph
                plot(X,Y,'.');
                hold on;
                plot(cutX,bestModel.Fitted,'k','LineWidth',4);
            end
            if consoleOutput
                disp("Section: " + section);
                disp("Time interval: " + timeDescription);
                disp("Formula: " + bestModel.Formula.Expression);
                disp("Coefficients: " + coefficients);
                missing = ismissing(bestModel.MSE);
                if missing(1)
                    disp("MSE: 0");
                else
                    disp("MSE: " + bestModel.MSE);
                end
                disp("R^2 value ordinary/adjusted: " + bestModel.Rsquared.Ordinary + "/" + bestModel.Rsquared.Adjusted);
                disp(" ");
            end
            if fileOutput
            fprintf(fid, "Section: " + section + "\n");
            fprintf(fid, "Time interval: " + timeDescription + "\n");
            fprintf(fid, "Formula: " + bestModel.Formula.Expression + "\n");
            fprintf(fid, "Coefficients: " + coefficients + "\n");
            missing = ismissing(bestModel.MSE);
            if missing(1)
                fprintf(fid, "MSE: 0\n");
            else
                fprintf(fid, "MSE: " + bestModel.MSE + "\n");
            end
            fprintf(fid, "R^2 value ordinary/adjusted: " + bestModel.Rsquared.Ordinary + "/" + bestModel.Rsquared.Adjusted + "\n");
            fprintf(fid, "\n");
            end
        end
        if showRegressionGraphs
            % Finalize the plot attributes
            readableValueName = strrep(valueName,'_',' ');
            readableTimeName = strrep(timeAxisIdentifier,'_',' ');
            title("Regression for: " + readableValueName);
            xlabel(readableTimeName);
            ylabel(readableValueName);
            set(gca,'FontSize',30);
            lines = table2array(timeAxis);
            xline(lines(borders),'r');
            hold off;
        end
        if fileOutput
            fclose(fid);
        end
    end
end


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