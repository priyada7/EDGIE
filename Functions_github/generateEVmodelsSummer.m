function [ec,eMin,eMax,e0,alpha,groceryEnergy,errandEnergy] = generateEVmodelsSummer...
    (n2,nDays,K,theta,batteryDeg,dt,t,percentageSedans,commuteDistanceProfile,sedan_idx,pickup_idx)
% this function is used to generate input parameters for electric vehicle
% simulation
%
% Input:
%  n2, number of homes
%  nDays, number of days
%  K, number of time steps
%  theta, Kx1 vector of outdoor temperature, C
%  batteryDeg, indicator for battery degraation effects 
%  dt, time sep, h
%  t, (K+1)x1 vector of time span, h
%
% Output:
%  a2, 1xn1 matrix of discrete-time dynamics parameters
%  etac, 1xn2 vector of charge efficiency
%  etad, 1xn2 vector of discharge efficiency 
%  eMin, n2xK matrix of min-energy capacity, kWh 
%  eMax, n2xK matrix of max-energy capacity, kWh 
%  pcMax, n2xK matrix of charge capacity, kW
%  pdMax, n2xK matrix of discharge capacity, kW
%  tau, KxL matrix of dissipation rate, 1/h
%  e0, initial energy of ev battery, kWh
%  ec, 1xn2 vector of commute energy, kWh
%  atHome,indicator that vehicle's at home
%  atWork, indicator that vehicle's at work
%  onRoad, indicator that vehicle's on the road
%  tw, 1xn2 vector for commute to work time of day, h
%  th, 1xn2 vector forcommute to home time of day, h
sedans = round(percentageSedans*n2);
pickupTruck = n2-sedans;
etac = trirnd(0.9,0.95,1,n2);       % charge efficiency

% electric vehicle commutes
tw = round(trirnd(5,11,1,n2));      % commute to work time of day, h
th = round(trirnd(15,23,1,n2));     % commute to home time of day, h
if batteryDeg ==1
   commuteDistance = commuteDistanceProfile;% miles
    % source for drivng efficiency modeling https://doi.org/10.1016/j.trd.2019.07.025
  %  drivingEfficiency = zeros(1,nDays);
    for i=1:length(theta)
        if theta(i,1) < 22
            drivingEfficiency(i) = 0.3392 - 0.005238*theta(i,1) - 0.0001078*theta(i,1)^2 + 1.04710e-5*theta(i,1)^3 + 3.955e-7*theta(i,1)^4-1.362e-8*theta(i,1)^5 - 3.109e-10*theta(i,1)^6;
        else
            drivingEfficiency(i)  = 0.4211-0.01627*theta(i,1) + 0.0004229*theta(i,1)^2;
        end
    end
    sedanEff = drivingEfficiency-0.0318;
    pickupTruckEff = drivingEfficiency+0.22;
    
     sedans = round(percentageSedans*n2);
     pickupTruck = n2-sedans;
    
    alpha_sedan = trirnd(min(sedanEff),max(sedanEff),1,n2);
    alpha_pickup = trirnd(min(pickupTruckEff),max(pickupTruckEff),1,n2);
    
    ec1=commuteDistance.*alpha_sedan;
    ec2=commuteDistance.*alpha_pickup;
    % Generate random indices for selecting columns from ec1
    

    % Create ec by selecting columns from ec1 and ec2 based on the generated indices
    ec = [ ec2(:, pickup_idx),ec1(:, sedan_idx)];

    %% commute distance calc
     % Preallocate
     d_commute = zeros(1, n2);

     

     % Assign commute distance to sedans and pickups
     d_commute(sedan_idx)  = commuteDistance(sedan_idx);
     d_commute(pickup_idx) = commuteDistance(pickup_idx);
      alpha = [alpha_pickup(pickup_idx),alpha_sedan(sedan_idx)];

else
    ec = trirnd(3,4,1,n2);          % commute energy, kWh (12 miles, 0.3 miles/kWh)
end
eMax_sedan = 50 + trirnd(0,40,1,sedans);              % energy capacity, kWh 
eMax_pickup = 120 + trirnd(0,50,1,pickupTruck);

eMax =[eMax_pickup,eMax_sedan];

groceryEnergySedan = trirnd(2*2, 2*10, 1, n2).*alpha_sedan;  % 
groceryEnergyPickup = trirnd(2*2, 2*10, 1, n2).*alpha_pickup;  % 


errandEnergySedan = trirnd(2*11, 2*15, 1, n2).*alpha_sedan;  % 
errandEnergyPickup = trirnd(2*11, 2*15, 1, n2).*alpha_pickup;  % 

errandEnergy = [errandEnergyPickup(pickup_idx),errandEnergySedan(sedan_idx)];
groceryEnergy = [groceryEnergyPickup(pickup_idx),groceryEnergySedan(sedan_idx)];

eMin = trirnd(0.40,0.50,1,n2).*eMax;        % user-specified minimum energy, kWh


% compute total required trip energy
requiredEnergy = (ec + groceryEnergy + errandEnergy) ./ etac;

% raise eMin if too low
eMin = max(eMin, requiredEnergy);

% cap it at battery capacity
eMin = min(eMin, eMax);

% safety check
idxBad = find(eMin > eMax);
if ~isempty(idxBad)
    error('eMin > eMax at indices: %s', mat2str(idxBad));
end
infeasibleEV=[];
e0 = eMin + rand(1,n2).*(eMax-eMin);  % initial energy, kWhend
end
