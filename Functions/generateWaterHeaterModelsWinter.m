function [R,a,w,theta,Tset,eta,pMaxHP,pMaxR] ...
    = generateWaterHeaterModelsWinter(waterFile,L,dt,tStart,tEnd)
% generateWaterHeaterModels generates thermal models of heat-pump water
% heaters, based primarily on water withdrawal data from Hendron's NREL
% tool.

% Input:
%   waterFile, the DHW Event Generator CSV file name
%   L, the number of buildings to model
%   timeSpan, the datetime span
%
% Output:
%   timeSpan, a Kx1 time span, datetime object
%   R, a 1xL vector of thermal resistances, C/kW
%   a, a 1xL vector of discrete-time dynamics parameters
%   w, a KxL matrix of exogenous thermal powers, kW
%   theta, a Kx1 vector of outdoor temperatures, C
%   Tset, a KxL matrix of indoor temperature setpoints, C
%   eta, a KxL matrix of heat pump coefficients of performance
%   pMisc, a KxL matrix of electric powers from miscellaneous loads, kW
%   p, a 1xL vector of maximum electric powers, kW
%   H, a Kx1 vector of indicators of the heating season

startTime = datetime(2021,1,1,0,0,0);           % start time
endTime = datetime(2022,1,1,0,0,0) - hours(dt); % end time
timeSpan = (startTime:hours(dt):endTime)';      % time span as datetime

%% set high-level parameters
K = length(timeSpan);
dt = hours(timeSpan(2) - timeSpan(1));
t = (0:dt:365*24)'; % time span as float, h
Tset = f2c(trirnd(120,130,1,L)); % tank water temperature, C

%% import thermal power withdrawal (from Hendron's NREL tool)
w0 = importWater(waterFile,timeSpan);

%% shuffle days to randomize withdrawal patterns over loads
% data storage
w = zeros(K,L); % thermal power withdrawals, one column per load, kW

% shuffle
w0 = reshape(w0,24/dt,365); % reshape one-load data to one column per day
for i=1:L
    w(:,i) = vecX(w0(:,randperm(365))); % shuffle columns and store the result as one stacked column
end

%% define tank parameters
% thermal capacitance
gallonsPerCubicMeter = 264.17; % gallons per cubic meter, gal/m^3
waterDensity = 997; % water density, kg/m^3
waterSpecificHeat = 0.001163056; % specific heat of water, kWh/kg/C
tankVolume = (50 + randi([0,1],1,L)*30)/gallonsPerCubicMeter; % tank volume, either 50 or 80 gallons, converted to m^3
C = waterDensity*waterSpecificHeat*tankVolume; % tank thermal capacitance, kWh/C

% thermal resistance
tankRadius = 0.25; % tank radius, m
tankArea = 1.1*tankVolume/tankRadius; % tank vertical surface area, m^2, from V = pi*r^2*h and A = 2*pi*r*h = 2*V/r
englishTankRValue = 6; % tank R-value in English units, F*ft^2/BTU/h
celsiusPerFahrenheit = 5/9;
squareMetersPerSquareFoot = 1/10.8;
kWPerBTUh = 1/3412;
metricTankRValue = englishTankRValue*... % tank R-Value in metric units, C*m^2/kW (1060-1410)
    celsiusPerFahrenheit*squareMetersPerSquareFoot/kWPerBTUh;
R = metricTankRValue./tankArea; % tank thermal resistance, C/kW
a = exp(-dt./(R.*C)); % discrete-time dynamics parameter
eta = repmat(3,K,L);


theta = repmat(f2c(60),K,1); % ambient temperature surrounding tank, C theta =f2c(60); % modified by pd on 27th march,2023
pMaxHP = repmat(0.5,1,L); % pMaxHP = 0.5; % heat pump electrical capacity, kW
pMaxR = repmat(9,1,L); % 

for i=1:length(timeSpan)
    if tStart == timeSpan(i)
        startIndex = i;
    end
    if tEnd == timeSpan(i)
        endIndex = i;
    end
end

%match the time step with weather data
eta = eta(startIndex:endIndex-1,:);
w = w(startIndex:endIndex-1,:);
theta = theta(startIndex:endIndex-1,:);



end

