function [ec] = generateEVmodelsSummer...
    (n2,nDays,K,theta,batteryDeg,dt,t,commuteDistance,percentageSedans)
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



% electric vehicle commutes
tw = round(trirnd(5,11,1,n2));      % commute to work time of day, h
th = round(trirnd(15,23,1,n2));     % commute to home time of day, h
if batteryDeg ==1
    commuteDistance = trirnd(0.9*commuteDistance,1.1*commuteDistance,1,n2);% miles
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
    
    
    ec1=commuteDistance.*trirnd(min(sedanEff),max(sedanEff),1,n2);
    ec2=commuteDistance.*trirnd(min(pickupTruckEff),max(pickupTruckEff),1,n2);

   sedans = round(percentageSedans*n2);
   pickupTruck = n2-sedans;

    % Generate random indices for selecting columns from ec1
    indices_ec1 = randperm(size(ec1, 2),sedans);

    % Create ec by selecting columns from ec1 and ec2 based on the generated indices
    ec = [ec1(:, indices_ec1), ec2(:, setdiff(1:size(ec2, 2), indices_ec1))];
else
    ec = trirnd(3,4,1,n2);          % commute energy, kWh (12 miles, 0.3 miles/kWh)
end
atHome = zeros(K,n2);               % indicator that vehicle's at home
atWork = zeros(K,n2);               % indicator that car's at work
for i=1:n2
    atHome(:,i) = mod(t(1:K),24)<tw(i) | mod(t(1:K),24)>th(i);
    atWork(:,i) = mod(t(1:K),24)>tw(i) & mod(t(1:K),24)<th(i);
end




end