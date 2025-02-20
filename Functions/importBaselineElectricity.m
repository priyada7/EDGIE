function [P,summerPeak,fullP,powerData] = importBaselineElectricity(fileName,tStart,tEnd,weatherTime,n1,desiredState)
% this function performs the baseline or no electrical load calculation 
%
% Input:
%  filename, name of the load file 
%  tStart, start date for simulation
%  tEnd, end date for simulation
%  weatherTime, timespan of weather file
%  n1, number of homes
%
% Output:
%  P, Kx1 vector for neighnorhood base electrical power, kW
%  summerPeak, summer peak aggregate power, kW

fileName = 'cleanedMFRED.xlsx';  

opts = detectImportOptions(fileName);
opts.PreserveVariableNames = 1;
rawData = readtable(fileName,opts);

% % extract and retime aggregate power data
 powerTime = rawData{:,1};                                   % time stamp
if powerTime.Year(1) < 2000
    powerTime.Year = 2000+powerTime.Year;                   % fix year from e.g. 20 to 2020
end

% %to match the timeline of 2021 weather data
 if powerTime.Year(1) < 2022
     powerTime.Year = powerTime.Year+2;                   % fix year from e.g. 20 to 2020
 end

powerTime = powerTime - hours(5);                           % convert UTC to eastern
iPower = 2:3:size(rawData,2);                               % indices of kW columns
individualPower = rawData{:,iPower};                        % individual power profiles
aggregatePower = sum(individualPower,2);
powerData = timetable(powerTime,aggregatePower);
powerData = fillmissing(powerData,'linear');
powerData = retime(powerData,weatherTime);
powerData = fillmissing(powerData,'linear');

% rescale aggregate power data


West = {'Montana', 'Idaho', 'Wyoming' ,'Nevada', 'Utah', 'Colorado', 'Arizona', 'New Mexico','Washington', 'Oregon', 'California', 'Alaska', 'Hawaii'};
MidWest = {'Ohio', 'Indiana', 'Illinois', 'Michigan', 'Wisconsin','Minnesota', 'Iowa', 'Missouri', 'North Dakota', 'South Dakota', 'Nebraska', 'Kansas'};
NorthEast ={'Maine', 'New Hampshire', 'Vermont', 'Massachusetts', 'Rhode Island', 'Connecticut','New York', 'Pennsylvania', 'New Jersey'};

if any(strcmpi(desiredState, West))
    annualElectricityPerHome = (175e12-15e12-14e12-33e12)/16.97e6;% RECS 2020 Table CE4.5 annual kWh per single-family detached home in west pg 4
elseif any(strcmpi(desiredState, MidWest))
    annualElectricityPerHome = (206e12-28e12-20e12-27e12)/22.09e6;
elseif any(strcmpi(desiredState, NorthEast))
    annualElectricityPerHome = (114e12-11e12-10e12-17e12)/11.23e6;
else
    annualElectricityPerHome = (459e12-58e12-53e12-112e12)/30.29e6;
end


neighborhoodElectricity = n1*annualElectricityPerHome;      % neighborhood annual kWh
powerData.aggregatePower = neighborhoodElectricity*powerData.aggregatePower/sum(powerData.aggregatePower);


% extract a representative load period
P = powerData.aggregatePower(powerData.powerTime>=tStart & powerData.powerTime<tEnd);


% rescale aggregate power data to match E. Bitar's Ithaca pilot, https://www.eesi.org/files/Eilyan_Bitar_Slides_20210625.pdf
normP = (P-min(P))/(max(P)-min(P));
minP = n1*80/35;                                            % minimum power on selected winter day
maxP = n1*120/35;                                           % maximum power on selected winter day
fullP = powerData.aggregatePower;                           % extract full year of unscaled power
fullP = (fullP-min(P))/(max(P)-min(P));                     % apply same transforms to fullP as selected winter power, P
fullP = minP + (maxP-minP)*fullP;
summerPeak = quantile(fullP,0.999);                         % summer peak aggregate power, kW
P = minP + (maxP-minP)*normP;

end