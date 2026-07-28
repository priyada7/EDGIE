function [pWorkBase] = importBaseElectricityfromWorkPlace(n1,K,t)
% aggregate baseline power of work neighborhood, kW
workAnnualElectricityIntensity = 17;                        % kWh/sqft/yr, from CBECS 2012 Table PBA4
workSqftPerPerson = 250;
workMeanPower = 2*n1*workSqftPerPerson*workAnnualElectricityIntensity/8760; % kW
workLoadFactor = 0.33;                                      % ratio of peak power to mean power
workPeakPower = workMeanPower/workLoadFactor; % kW
rho = 0.7;                                                  % ratio of daily mean power to daily peak power
pWorkBase = rho*workPeakPower + (1-rho)*workPeakPower*sin(2*pi*(t(1:K)-5)/24);
end