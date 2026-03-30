function [a2,etac,etad,eMin,eMax,pcMax,pdMax,tau,e0,p0,ec,atHome,atWork,onRoad,tw,th,commuteDistanceProfile,sedan_idx,pickup_idx,alpha] = generateEVmodelsWinter...
    (n2,nDays,K,theta,batteryDeg,dt,t,commuteDistance,percentageSedans)
% this function is used to generate input parameters for electric vehicle
% simulation
%

%--------------------------------------------------------------------------
% INPUT ARGUMENTS
%--------------------------------------------------------------------------
% n2 : scalar (integer)
%     Number of EVs in the system.
%
% nDays : scalar (integer)
%     Number of simulated days.
%
% K : scalar (integer)
%     Total number of time steps (K = 24 × nDays for hourly models).
%
% theta : K × 1 vector
%     Ambient temperature profile (°C), used for winter driving efficiency.
%
% batteryDeg : scalar {0,1}
%     Battery degradation / temperature‑dependent efficiency flag:
%       0 → constant energy consumption
%       1 → temperature‑dependent driving efficiency
%
% dt : scalar
%     Length of each time step (hours).
%
% t : K × 1 vector
%     Time index (in hours), used to determine hour‑of‑day behavior.
%
% commuteDistance : scalar
%     Nominal one‑way commute distance (miles).
%
% percentageSedans : scalar in [0,1]
%     Fraction of EVs that are sedans (rest are pickup trucks).
%
%--------------------------------------------------------------------------
% OUTPUT ARGUMENTS
%--------------------------------------------------------------------------
% a2 : n2 × 1 vector
%     Discrete‑time battery self‑discharge coefficient.
%
% tau : n2 × 1 vector
%     Battery dissipation time constant (hours).
%
% etac, etad : n2 × 1 vectors
%     Charging and discharging efficiencies.
%
% pcMax : n2 × 1 vector
%     Maximum charging power capacity (kW).
%
% pdMax : n2 × 1 vector
%     Maximum discharging power capacity (kW).
%
% eMax : 1 × n2 vector
%     Maximum battery energy capacity (kWh).
%
% eMin : 1 × n2 vector
%     Minimum required battery energy (kWh), ensuring commute feasibility.
%
% e0 : 1 × n2 vector
%     Initial battery state of charge (kWh).
%
% p0 : 1 × n2 vector
%     Initial charging power (kW); only a subset of vehicles charge initially.
%
% ec : 1 × n2 vector
%     Energy consumption per commute (kWh).
%
% alpha : 1 × n2 vector
%     Vehicle‑specific driving energy intensity (kWh/mile).
%
% commuteDistanceProfile : 1 × n2 vector
%     Realized commute distances for all EVs (miles).
%
%--------------------------------------------------------------------------
% LOCATION AVAILABILITY MATRICES
%--------------------------------------------------------------------------
% atHome : K × n2 logical matrix
%     atHome(k,i) = 1 if EV i is at home during time step k.
%
% atWork : K × n2 logical matrix
%     atWork(k,i) = 1 if EV i is at work during time step k.
%
% onRoad : K × n2 logical matrix
%     onRoad(k,i) = 1 if EV i is driving during time step k.
%
%--------------------------------------------------------------------------
% COMMUTE TIMING
%--------------------------------------------------------------------------
% tw : 1 × n2 vector
%     Departure time to work (hour of day).
%
% th : 1 × n2 vector
%     Departure time to home (hour of day).
%
%--------------------------------------------------------------------------
% VEHICLE TYPE INDICES
%--------------------------------------------------------------------------
% sedan_idx : vector
%     Indices of EVs classified as sedans.
%
% pickup_idx : vector
%     Indices of EVs classified as pickup trucks.
%


% electric vehicle batteries
rho = trirnd(0.97,0.99,1,n2);       % fraction of charge remaining after 24 h
tau = -24./log(rho);                % dissipation rate, 1/h
a2 = exp(-dt./tau);                 % discrete-time dynamics parameter
etac = trirnd(0.9,0.95,1,n2);       % charge efficiency
etad = etac;                        % discharge efficiency
pcMax = 240*trirnd(24,48,1,n2)/1e3; % charge capacity, kW (level 2: 240 V, 24-32 A)

sedans = round(percentageSedans*n2);
pickupTruck = n2-sedans;
pdMax = trirnd(30,40,1,n2);                 % discharge capacity, kW (100 mile/h at 0.35 kWh/mile: 35 kW)


%e0 = eMin + trirnd(0,1,1,n2).*(eMax-eMin);  % initial energy, kWh
daysPerCharge = 5.5;                        % average number of days between full charges
 
%p0 = (rand(1,n2)<=1/daysPerCharge).*pcMax;  % initial charging power, kW

numToCharge = round(n2 / daysPerCharge);  % Calculate the number of cars to charge today 
p0 = zeros(1, n2);
selectIndicesp0 = randperm(n2, numToCharge); % Randomly select cars to be charged
p0(selectIndicesp0) = pcMax(selectIndicesp0);% Assign pcMax to the selected cars


% electric vehicle commutes
tw = round(trirnd(5,11,1,n2));      % commute to work time of day, h
th = round(trirnd(15,23,1,n2));     % commute to home time of day, h





atHome = zeros(K,n2);
atWork = zeros(K,n2);


hoursPerDay = 24;
for d = 1:nDays
    dayStart = (d-1)*hoursPerDay + 1;
    dayEnd = d*hoursPerDay;
    dayIndices = dayStart:dayEnd;
    
    % 30% of vehicles stay home all day
    stayHomeToday = randperm(n2, round(0.3*n2));
    
    % Mark them at home all 24h
    atHome(dayIndices, stayHomeToday) = 1;
    
    % The rest follow normal commute patterns
    commutesToday = setdiff(1:n2, stayHomeToday);
   % EVsAlwaysHomeEachDay{d} = stayHomeToday; % <--- for plotting
    for i = commutesToday
        for h = dayIndices
            hourOfDay = mod(t(h),24);
            if hourOfDay < tw(i) || hourOfDay > th(i)
                atHome(h,i) = 1;
            elseif hourOfDay > tw(i) && hourOfDay < th(i)
                atWork(h,i) = 1;
            end
        end
    end
    
end


onRoad = ~atHome & ~atWork; % indicator that vehicle's on the road

if batteryDeg ==1
    commuteDistance = trirnd(0.9*commuteDistance,1.1*commuteDistance,1,n2);% miles
    % source for drivng efficiency modeling https://doi.org/10.1016/j.trd.2019.07.025
    drivingEfficiency = zeros(1,length(theta));
    for i=1:length(theta)
        if theta(i,1) < 22
            drivingEfficiency(i) = 0.3392 - 0.005238*theta(i,1) - 0.0001078*theta(i,1)^2 + 1.04710e-5*theta(i,1)^3 + 3.955e-7*theta(i,1)^4-1.362e-8*theta(i,1)^5 - 3.109e-10*theta(i,1)^6;
        else
            drivingEfficiency(i)  = 0.4211-0.01627*theta(i,1) + 0.0004229*theta(i,1)^2;
        end
    end
    sedanEff = drivingEfficiency-0.0318;
    pickupTruckEff = drivingEfficiency+0.22;
    
    alpha_sedan = trirnd(min(sedanEff),max(sedanEff),1,n2);
    alpha_pickup = trirnd(min(pickupTruckEff),max(pickupTruckEff),1,n2);
    
    ec1=commuteDistance.*alpha_sedan;
    ec2=commuteDistance.*alpha_pickup;
 commuteDistanceProfile=commuteDistance;

   % Generate random indices for selecting columns from ec1
     indices_ec1 = randperm(size(ec1, 2),sedans);
 
     % Create ec by selecting columns from ec1 and ec2 based on the generated indices
     ec = [ec2(:, setdiff(1:size(ec2, 2), indices_ec1)),ec1(:, indices_ec1)];

     %% commute distance calc
     % Preallocate
     d_commute = zeros(1, n2);

     sedan_idx  = indices_ec1;
     pickup_idx = setdiff(1:n2, indices_ec1);

     % Assign commute distance to sedans and pickups
     d_commute(sedan_idx)  = commuteDistance(sedan_idx);
     d_commute(pickup_idx) = commuteDistance(pickup_idx);
     alpha = [alpha_pickup(pickup_idx),alpha_sedan(sedan_idx)];

else
    ec = trirnd(3,4,1,n2);          % commute energy, kWh (12 miles, 0.3 kWh/miles)
end

eMax_sedan = 50 + trirnd(0,40,1,sedans);              % energy capacity, kWh 
eMax_pickup = 120 + trirnd(0,50,1,pickupTruck);






eMax =[eMax_pickup,eMax_sedan];
eMin = trirnd(0.20,0.30,1,n2).*eMax;        % user-specified minimum energy, kWh
% Enforce minimum energy based on efficiency
mask = eMin < (ec ./ etad);          % logical mask where eMin is too small
eMin(mask) = ec(mask) ./ etad(mask); % replace with minimum allowed
e0 = eMin + rand(1,n2).*(eMax-eMin);  % initial energy, kWhend
end
