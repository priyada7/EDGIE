function [P,fullP] = importBaselineElectricity(individualPower,tStart,tEnd,weatherTime,n1,desiredState,percentageAttached,percentageDetached)
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


West = {'Montana', 'Idaho', 'Wyoming' ,'Nevada', 'Utah', 'Colorado', 'Arizona', 'New Mexico','Washington', 'Oregon', 'California', 'Alaska', 'Hawaii'};
MidWest = {'Ohio', 'Indiana', 'Illinois', 'Michigan', 'Wisconsin','Minnesota', 'Iowa', 'Missouri', 'North Dakota', 'South Dakota', 'Nebraska', 'Kansas'};
NorthEast ={'Maine', 'New Hampshire', 'Vermont', 'Massachusetts', 'Rhode Island', 'Connecticut','New York', 'Pennsylvania', 'New Jersey'};

 if any(strcmpi(desiredState, West)) %https://www.eia.gov/consumption/residential/data/2020/c&e/pdf/ce4.2.pdf
    scaling_detached = (17e9+96e9)/(16.97e6*8760);
    scaling_attached = ((1e9+6e9)/(8760*1.69e6) + (1e9+5e9)/(8760*1.89e6) + (3e9+14e9)/(8760*5.70e6) )/3;
elseif any(strcmpi(desiredState, MidWest))
    scaling_detached = (18e9+115e9)/(8760*18.58e6);
    scaling_attached = ((1e9+6e9)/(8760*1.33e6) + (1e9+5e9)/(8760*1.95e6) + (2e9+10e9)/(8760*4.20e6) )/3;
elseif any(strcmpi(desiredState, NorthEast))
    scaling_detached = (10e9+66e9)/(11.23e6*8760);
    scaling_attached = ((2e9+8e9)/(8760*1.95e6) + (2e9+8e9)/(8760*3.15e6) + (3e9+11e9)/(8760*5.10e6) )/3;
else
    scaling_detached = (31e9+205e9)/(30.29e6*8760);
    scaling_attached = ((2e9+11e9)/(8760*2.48e6) + (1e9+8e9)/(8760*2.36e6) + (5e9+26e9)/(8760*7.83e6) )/3;
end

detachedIndividualPower = scaling_detached*individualPower.Power./mean(individualPower.Power);
attachedIndividualPower = scaling_attached*individualPower.Power./mean(individualPower.Power);

nDetached = round(n1*percentageDetached);

nMFRED= size(individualPower,2);

fullP = zeros(size(individualPower,1),n1);
for i = 1 : nDetached
fullP(:,i) = detachedIndividualPower(:,randi(nMFRED));
end
for i = nDetached + 1 : n1
fullP(:,i) = attachedIndividualPower(:,randi(nMFRED));
end
fullP = sum(fullP,2);
P = fullP(individualPower.powerTime>=tStart & individualPower.powerTime<tEnd); % extract a representative load period



end