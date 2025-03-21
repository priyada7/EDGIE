function [temperature,shortwave] = importW(fileName,t) %,shortwave,windspeed,humidity
% importWeather imports and processes weather data from a CSV file
% generated by OikoLab.
%
% Input:
%   fileName, the name of the OikoLab CSV weather file.
%   t, the datetime span.
%
% Output:
%   temperature, the outdoor temperature in C
%   shortwave, the global solar shortwave irradiance on a horizontal surface in kW/m^2


% import raw data
opts = detectImportOptions(fileName);
opts.PreserveVariableNames = 1;
weatherData = readtable(fileName,opts);

% extract data and convert units
timeStamp = weatherData{:,1}; % time stamp
if timeStamp(1).Year < 2000, timeStamp.Year = timeStamp.Year + 2000; end
timeZoneAdjust = weatherData{:,5};
timeStamp = timeStamp + hours(timeZoneAdjust); % convert UTC to eastern
temperature = weatherData{:,6}; % outdoor air temperature, C
shortwave = weatherData{:,7}/1000; % global horizontal shortwave irradiance, kW/m^2
% humidity = weatherData{:,9}; % relative humidity
% windspeed = weatherData{:,10}; % wind speed, m/s

% fill any missing data
temperature = fillmissing(temperature,'linear');
shortwave = fillmissing(shortwave,'linear');
% humidity = fillmissing(humidity,'linear');
% windspeed = fillmissing(windspeed,'linear');

% pack the data into a timetable object
TT = timetable(timeStamp,temperature,shortwave);%,windspeed,humidity);

% retime to the desired time steps
TT = retime(TT,t);
temperature = TT{:,1};
shortwave = TT{:,2};
% windspeed = TT{:,3};
% humidity = TT{:,4};

% fill any missing data again
temperature = fillmissing(temperature,'previous');
shortwave = fillmissing(shortwave,'previous');
% humidity = fillmissing(humidity,'linear');
% windspeed = fillmissing(windspeed,'linear');


end

