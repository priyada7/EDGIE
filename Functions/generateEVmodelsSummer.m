function [ec,eMin,e0,alpha] = generateEVmodelsSummer...
    (n2,nDays,K,theta,batteryDeg,dt,t,percentageSedans,eMax,commuteDistanceProfile,sedan_idx,pickup_idx)
% this function is used to generate input parameters for electric vehicle
% simulation
%


sedans = round(percentageSedans*n2);
pickupTruck = n2-sedans;
etac = trirnd(0.9,0.95,1,n2);       % charge efficiency

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



eMin = trirnd(0.20,0.30,1,n2).*eMax;        % user-specified minimum energy, kWh
mask = eMin < (ec ./ etac);          % logical mask where eMin is too small
eMin(mask) = ec(mask) ./ etac(mask); % replace with minimum allowed
e0 = eMin + rand(1,n2).*(eMax-eMin);  % initial energy, kWh
end