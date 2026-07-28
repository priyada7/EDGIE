function [a2,etac,etad,eMin,eMax,pcMax,pdMax,tau,e0,p0,ec,atHome,atWork,onRoad,tw,th,commuteDistanceProfile,sedan_idx,pickup_idx,alpha,atGrocery,groceryEnergy,atErrands,errandEnergy] = generateEVmodelsWinter...
    (n2,nDays,K,theta,batteryDeg,dt,t,commuteDistance,percentageSedans,design_temperatureWinter,design_temperatureSummer,warmupDays)
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
te = round(trirnd(5,22,1,n2));     % commute to home time of day, h

atHome = zeros(K,n2);
atWork = zeros(K,n2);
atErrands = zeros(K,n2);


maxTrips = 1;   % max errands per vehicle

errandsSetPerDay = cell(nDays,1);

hoursPerDay = 24;
for d = 1:nDays
    dayStart = (d-1)*hoursPerDay + 1;
    dayEnd = d*hoursPerDay;
    dayIndices = dayStart:dayEnd;

    
    numCommute = round(0.25*n2); %32
    numErrands = round(0.45*n2); %38
    perm = randperm(n2);
    commutesToday = perm(1:numCommute);
    errandsToday  = perm(numCommute+1:numCommute+numErrands);
    stayHomeToday = perm(numCommute+numErrands+1:end);
    errandsSetPerDay{d} = errandsToday;
    % Mark them at home all 24h
    atHome(dayIndices, stayHomeToday) = 1;
    
  
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

    
  for i = errandsToday
    atHome(dayIndices,i) = 1;

    numTrips = 1 ;%randi([1 3]);

    tStarts = randi([5,21],1,numTrips);
    tEnds   = tStarts+1;

    for j = 1:numTrips
        if tEnds(j) > tStarts(j)
            hIdx = dayStart + (tStarts(j):(tEnds(j)-1));
            atHome(hIdx,i) = 0;
            atErrands(hIdx,i) = 1;
        end
    end
end

    
end



%% Grocery runs


atGrocery = zeros(K, n2);
tGroceryOG = round(trirnd(5,22,1,n2)); 



for d = 1:nDays
    dayStart   = (d-1)*hoursPerDay + 1;
    dayIndices = dayStart : d*hoursPerDay;

    errandsToday = errandsSetPerDay{d};  

    for i = errandsToday

      tGrocery = tGroceryOG(i);

        for h = dayIndices
            hourOfDay = mod(t(h), 24);

          if hourOfDay == tGrocery
               
                hPrev = h - 1;
                hNext = h + 1;

                if hPrev >= 1 && hNext <= K && ...
                   atHome(hPrev, i) == 1 && ...
                   atHome(hNext, i) == 1 && ...
                   atHome(h, i) == 1

                    atHome(h, i)    = 0;
                    atGrocery(h, i) = 1;
                end
            end
        end

    end
end

onRoad = ~atHome & ~atWork & ~atErrands & ~atGrocery;


%%
if batteryDeg ==1
   commuteDistance = trirnd(0.9*commuteDistance,1.1*commuteDistance,1,n2);% miles
   %commuteDistance = trirnd(20,40,1,n2);

  
% 
% u = rand(1,n2);
% 
% idx1 = u <= 0.556;
% idx2 = u > 0.556 & u <= 0.878;
% idx3 = u > 0.878 & u <= 0.979;
% idx4 = u > 0.979;
% 
% % <10 miles
% commuteDistance(idx1) = trirnd(1,10,1,sum(idx1));
% 
% % 10-24 miles
% commuteDistance(idx2) = trirnd(10,24,1,sum(idx2));
% 
% % 25-50 miles
% commuteDistance(idx3) = trirnd(25,50,1,sum(idx3));
% 
% % >50 miles
% commuteDistance(idx4) = trirnd(50,55,1,sum(idx4));
   
    % source for drivng efficiency modeling https://doi.org/10.1016/j.trd.2019.07.025
    drivingEfficiency = zeros(1,length(theta));
    for i=1:length(theta)
        if theta(i,1) < 22
            drivingEfficiency(i) = 0.3392 - 0.005238*theta(i,1) - 0.0001078*theta(i,1)^2 + 1.04710e-5*theta(i,1)^3 + 3.955e-7*theta(i,1)^4-1.362e-8*theta(i,1)^5 - 3.109e-10*theta(i,1)^6;
        else
            drivingEfficiency(i)  = 0.4211-0.01627*theta(i,1) + 0.0004229*theta(i,1)^2;
        end
    end
    sedanEff = drivingEfficiency-0.0318; %0.09; %0.0318;
    pickupTruckEff = drivingEfficiency+0.22;% 0.10;%0.22;
    
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

groceryEnergySedan = trirnd(2*2, 2*10, 1, n2).*alpha_sedan;  % 
groceryEnergyPickup = trirnd(2*2, 2*10, 1, n2).*alpha_pickup;  % 


errandEnergySedan = trirnd(2*11, 2*15, 1, n2).*alpha_sedan;  % 
errandEnergyPickup = trirnd(2*11, 2*15, 1, n2).*alpha_pickup;  % 


errandEnergy = [errandEnergyPickup(pickup_idx),errandEnergySedan(sedan_idx)];
groceryEnergy = [groceryEnergyPickup(pickup_idx),groceryEnergySedan(sedan_idx)];

eMax =[eMax_pickup,eMax_sedan];
eMin = trirnd(0.40,0.45,1,n2).*eMax;        % user-specified minimum energy, kWh


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
