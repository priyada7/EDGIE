function [P,summerPeak] = importBaselineElectricity(fileName,tStart,tEnd,weatherTime,n1)
opts = detectImportOptions(fileName);
opts.PreserveVariableNames = 1;
rawData = readtable(fileName,opts);

% extract and retime aggregate power data
powerTime = rawData{:,1};                                   % time stamp
if powerTime.Year(1) < 2000
    powerTime.Year = 2000+powerTime.Year;                   % fix year from e.g. 20 to 2020
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
annualElectricityPerHome = 115e9/(10.8e6);                  % RECS 2015 Table CE4.2 annual kWh per single-family detached home in Northeast
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