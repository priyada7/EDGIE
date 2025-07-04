clear
tic
batteryDeg =   1; % battery degradation with outside temperature is implemented
matpower   =   0; % set 1 if voltage regulations needs to be solved
transformer=   0; % set 1 if thermal model of transformer needs to be solved
optimization = 0; % set 1 for v2h
waterheater = 3;  % set 1 for resitance 2 for heat pump only 3 for hybrid
sizing = 3;       % set 1 for cooling 2 for heating 3 for max of heating or cooling 
waterfile    = 'DHWEventGeneratorOutput.csv';                     % load water scheduler file

stateFolders = {  
'MN', 'LOCATION OF WEATHER FILE';
};

for idx = 1:size(stateFolders, 1)
    stateAbbr = stateFolders{idx, 1};
    outputFolder = stateFolders{idx, 2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% input data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% timing (hour 0 is noon to make EV simulation easier)
ti = 0;                      % initial time, h
nDays = 7;                   % no. of days
tf = nDays*24;               % final time, h
dt = 1;                      % time step, h
K = tf/dt;                   % number of time steps
t = (0:dt:tf)';              % time span, h
n1 = 1000;                   % number of homes (= number of HPs)
L = n1;                      % number of water heater
ft2m2 = 0.092903;

% Read the CSV file
data = readtable('version13.xlsx');

% Filter cities for state_id = 'AZ'
stateAZ = data(strcmpi(data.state_id, stateAbbr), :);

% Extract city names for state_id = 'AZ'
arizonacities = stateAZ.city_ascii;
USAstateName = stateAZ.state_name;
USAcountyNAme = stateAZ.county_name; 
USAstatelat = stateAZ.lat;
USAstatelng = stateAZ.lng;
coolingtemp = stateAZ.x1__CoolingTemp___F_;
heatingtemp = stateAZ.x99__HeatingTemp___F_;
twoWayCommute = stateAZ.commuteDistance_miles_;
UvalueWall = stateAZ.Uwall;
UvalueWindow = stateAZ.Uwindow;
AttachedHome = stateAZ.AttachedHome_;
DetachedHome = stateAZ.DetachedHome_;
zeroCarinHome = stateAZ.x0Vehicle/100;
oneCarinHome = stateAZ.x1Vehicle/100;
twoCarinHome = stateAZ.x2Vehicle/100;
threeCarinHome = stateAZ.x3_Vehicle/100;
DetachedFloorArea = stateAZ.DetachedFloorArea;
AttachedFloorArea = stateAZ.AttachedFloorArea;

DetachedFloorArea(isnan(DetachedFloorArea))=0;
AttachedFloorArea(isnan(AttachedFloorArea))=0;

DetachedFloorArea(DetachedFloorArea>3000)=3000;
DetachedFloorArea(DetachedFloorArea<1000)=1000;

AttachedFloorArea(AttachedFloorArea>1700)=1700;
AttachedFloorArea(AttachedFloorArea<450)=450;


TodaysHeadroom = trirnd(1.15,1.36,length(USAcountyNAme),1);
FutureHeadroom = 1.2;

houseElecWH = stateAZ.ElectricWH_;

sedans = stateAZ.Cars_;
desiredState = USAstateName(1,1);


data = cell(length(arizonacities), 2);

fileName = 'cleanedMFREDdata.xlsx';  

opts = detectImportOptions(fileName);
opts.PreserveVariableNames = 1;
rawData = readtable(fileName,opts);

% % extract and retime aggregate power data
powerTime = rawData{:,1};                                   % time stamp
weatherTime = (datetime(2021,1,1,0,0,0):hours(1):datetime(2022,1,1,0,0,0))';
   
%to match the timeline of 2021 weather data
if powerTime.Year(1) < 2022
    powerTime.Year = powerTime.Year+1+1;                   % fix year from e.g. 20 to 2020
end

powerTime = powerTime - hours(5);                           % convert UTC to eastern
Power = rawData{:,2:end};                        % individual power profiles
Power = fillmissing(Power, 'linear'); % Or 'constant', 0
individualPower = timetable(powerTime,Power);
individualPower = retime(individualPower,weatherTime);
individualPower = fillmissing(individualPower,'linear');

for stateIdx = 27 %1:length(arizonacities) 27 is for Minneapolis, MN chech indices in variable stateAZ
    cityName = arizonacities{stateIdx};
    stateName = USAstateName{stateIdx};
    countyName = USAcountyNAme{stateIdx};
    commuteDistance = twoWayCommute(stateIdx);
    lat = USAstatelat(stateIdx);
    lng = USAstatelng(stateIdx);
    Uwall = UvalueWall(stateIdx);
    Uwindow = UvalueWindow(stateIdx);
   
    AreaAttached = ft2m2.*AttachedFloorArea(stateIdx);
    AreaDetached = ft2m2.*DetachedFloorArea(stateIdx);

   [RvalueDetached,RvalueAttached,floorAreaDetached,floorAreaAttached] = Rcalc(Uwall,Uwindow,AreaDetached,AreaAttached,n1);


   homeWithElectricWH = houseElecWH(stateIdx);
   percentageSedans = sedans(stateIdx);
   percentageAttached = AttachedHome(stateIdx);
   percentageDetached = DetachedHome(stateIdx);

   n2_0vehicle = zeroCarinHome(stateIdx);
   n2_1vehicle = oneCarinHome(stateIdx);
   n2_2vehicle = twoCarinHome(stateIdx);
   n2_3vehicle = threeCarinHome(stateIdx);


selectedHeadroom = TodaysHeadroom(stateIdx);

houseHPdetached = stateAZ.HP_Detached(stateIdx);
houseAuxdetached = stateAZ.Aux_Detached(stateIdx);

houseHPattached = stateAZ.HP_Attached(stateIdx);
houseAuxattached = stateAZ.Aux_Attached(stateIdx);
housingUnits = stateAZ.HousingUnits(stateIdx);
zone = stateAZ.Zone(stateIdx);
   
   % Replace spaces with underscores in the city name
    %cityName = strrep(cityName, ' ');
    weatherfileName = fullfile(outputFolder, sprintf('%s_HistoricalWeather.csv', cityName));
    weatherTime = datetime(2021, 1, 1, 0, 0, 0):hours(dt):datetime(2021, 12, 31, 23, 0, 0);
   [thetaFull, ~] = importW(weatherfileName, weatherTime);

    % Check if the file exists before attempting to import weather data
    if exist(weatherfileName, 'file') == 2         
        data{stateIdx,1} =  stateName;
        data{stateIdx, 2} = countyName;
        data{stateIdx, 3} = cityName;
        data{stateIdx, 4} = lat;
        data{stateIdx,5} = lng;
        data{stateIdx, 6} = mean(RvalueDetached);
        data{stateIdx, 7} = mean(RvalueAttached);

    else
        % File doesn't exist, skip processing and move to the next city
        fprintf('File not found for %s\n', cityName);
        continue; % Skip to the next iteration
    end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% download weather & basline load data %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Import weather data
weatherTime = (datetime(2021,1,1,0,0,0):hours(1):datetime(2022,1,1,0,0,0))';
[thetaFull, solarFull] = importW(weatherfileName, weatherTime);

% Reshape the data into a matrix with 365 rows (days) and 24 columns (hours)
daily_temperature_matrix = reshape(thetaFull(1:8760), 24, [])';


% Calculate the number of full weeks in the data
num_full_weeks = floor(size(daily_temperature_matrix, 1) / 7);

% Trim the data to keep only full weeks
weekly_temperature_matrix = daily_temperature_matrix(1:num_full_weeks * 7, :);

% Reshape the data into a matrix with 168 columns (24 hours/day * 7 days/week)
weekly_temperature_matrix = reshape(weekly_temperature_matrix', 168, [])';


% Calculate the minimum and maximum temperature for each week
weekely_min_temperature = min(weekly_temperature_matrix, [], 2);
weekely_max_temperature = max(weekly_temperature_matrix, [], 2);




data{stateIdx,8}=heatingtemp(stateIdx);
data{stateIdx,9}=coolingtemp(stateIdx);


% Design temperature 
design_temperatureWinter = f2c(heatingtemp(stateIdx)); 
design_temperatureSummer = f2c(coolingtemp(stateIdx)); 
% Find five days with daily minimum temperatures closest to the design temperature
[sortedTemp, sorted_indices] = sort(abs(weekely_min_temperature - design_temperatureWinter));
selected_days_indices = sorted_indices(1:5);
selected_days_temperatures = weekely_min_temperature(selected_days_indices,:);

[sortedTempCooling, sorted_indicesCooling] = sort(abs(weekely_max_temperature - design_temperatureSummer));
selected_days_indicesCooling = sorted_indicesCooling(1:5);
selected_days_temperaturesCooling = weekely_max_temperature(selected_days_indicesCooling,:);



tStartWinter = datetime(2021,1,1) + days((selected_days_indices(1)-1)*7);
tStartSummer = datetime(2021,1,1) + days((selected_days_indicesCooling(1)-1)*7);


tEndWinter = tStartWinter + hours(tf); 
tEndSummer = tStartSummer + hours(tf); 

ttWinter = (tStartWinter:hours(dt):tEndWinter)';
thetaWinter = thetaFull(weatherTime>=tStartWinter & weatherTime<tEndWinter);
solarWinter = solarFull(weatherTime>=tStartWinter & weatherTime<tEndWinter);


ttSummer= (tStartSummer:hours(dt):tEndSummer)';
thetaSummer = thetaFull(weatherTime>=tStartSummer & weatherTime<tEndSummer);
solarSummer = solarFull(weatherTime>=tStartSummer & weatherTime<tEndSummer);

%% import 'no electrification' load profile
[Pwinter,fullPwinter] = importBaselineElectricity(individualPower,tStartWinter,tEndWinter,weatherTime,n1,desiredState,percentageAttached,percentageDetached);
[Psummer,fullPsummer] = importBaselineElectricity(individualPower,tStartSummer,tEndSummer,weatherTime,n1,desiredState,percentageAttached,percentageDetached);

startTime = datetime(2021,1,1,0,0,0);           % start time
endTime = datetime(2022,1,1,0,0,0) - hours(dt); % end time
timeSpan = (startTime:hours(dt):endTime)';      % time span as datetime
tStart = startTime;
tEnd = endTime;

[pWorkBase] = importBaseElectricityfromWorkPlace(n1,K,t);

%% Heat pump sizing
       for i=1:length(coolingtemp)
           designTempCool = coolingtemp(i,1);
           check = 0; %identifier for soalr adjustments
           [Hp_sizeCooling_detached]=Hp_coolingLoad(thetaFull,solarFull,fullPsummer,n1,mean(RvalueDetached),designTempCool,mean(floorAreaDetached),check);
           check =1; %identifier for soalr adjustments
           [Hp_sizeCooling_attached]=Hp_coolingLoad(thetaFull,solarFull,fullPsummer,n1,mean(RvalueAttached),designTempCool,mean(floorAreaAttached),check);
           data{stateIdx,10}=Hp_sizeCooling_detached;
           data{stateIdx,11}=Hp_sizeCooling_attached;
       end
 
     
       for i=1:length(heatingtemp)
           designTempHeat = heatingtemp(i,1);
            check = 0; %identifier for soalr adjustments
           [Hp_sizeHeating_detached]=Hp_heatingLoad(thetaFull,solarFull,fullPwinter,n1,mean(RvalueDetached),designTempHeat,mean(floorAreaDetached),zone,check);
            check = 1;%identifier for soalr adjustments
            [Hp_sizeHeating_attached]=Hp_heatingLoad(thetaFull,solarFull,fullPwinter,n1,mean(RvalueAttached),designTempHeat,mean(floorAreaAttached),zone,check);
           data{stateIdx,12}=Hp_sizeHeating_detached;
           data{stateIdx,13}=Hp_sizeHeating_attached;
       end

       if sizing == 1
           hpSizeAttached = max(Hp_sizeCooling_attached);
           hpSizeDetached = max(Hp_sizeCooling_detached);
           [capLow_Detached_cooling,capHigh_Detached_cooling,capLow_Detached_heating,capHigh_Detached_heating]= select_best_hvac_Cooling(design_temperatureSummer,hpSizeDetached);
           [capLow_Attached_cooling,capHigh_Attached_cooling,capLow_Attached_heating,capHigh_Attached_heating]=select_best_hvac_Cooling(design_temperatureSummer,hpSizeAttached);
       elseif sizing == 2
           hpSizeAttached = max(Hp_sizeHeating_attached);
           hpSizeDetached = max(Hp_sizeHeating_detached);
           [capLow_Detached_cooling,capHigh_Detached_cooling,capLow_Detached_heating,capHigh_Detached_heating] = select_best_hvac_Heating(design_temperatureWinter,hpSizeDetached);
           [capLow_Attached_cooling,capHigh_Attached_cooling,capLow_Attached_heating,capHigh_Attached_heating] = select_best_hvac_Heating(design_temperatureWinter,hpSizeAttached);
       elseif sizing == 3
           if (Hp_sizeHeating_attached>Hp_sizeCooling_attached)
               hpSizeAttached =Hp_sizeHeating_attached;
               [capLow_Attached_cooling,capHigh_Attached_cooling,capLow_Attached_heating,capHigh_Attached_heating] = select_best_hvac_Heating(design_temperatureWinter,hpSizeAttached);
           else
                hpSizeAttached =Hp_sizeCooling_attached;
                [capLow_Attached_cooling,capHigh_Attached_cooling,capLow_Attached_heating,capHigh_Attached_heating]=select_best_hvac_Cooling(design_temperatureSummer,hpSizeAttached);
           end

           if Hp_sizeHeating_detached>Hp_sizeCooling_detached
               hpSizeDetached = Hp_sizeHeating_detached;
               [capLow_Detached_cooling,capHigh_Detached_cooling,capLow_Detached_heating,capHigh_Detached_heating] = select_best_hvac_Heating(design_temperatureWinter,hpSizeDetached);
           else
               hpSizeDetached = Hp_sizeCooling_detached;
               [capLow_Detached_cooling,capHigh_Detached_cooling,capLow_Detached_heating,capHigh_Detached_heating]= select_best_hvac_Cooling(design_temperatureSummer,hpSizeDetached);
           end
       end





%% generate building and hp parameters
    [a1_detached,C_detached,R_detached,eta1_detached,Tset_detached,qeSummer_detached,pMax1_detached,pMaxAux_detached] = generateBuildingModelsSummer(n1,dt,K,thetaSummer,t,Psummer,solarSummer,RvalueDetached,hpSizeDetached,floorAreaDetached);
    [eta2_detached,qeWinter_detached] = generateBuildingModelsWinter(zone,n1,K,thetaWinter,Pwinter,solarWinter,floorAreaDetached);


    [a1_attached,C_attached,R_attached,eta1_attached,Tset_attached,qeSummer_attached,pMax1_attached,pMaxAux_attached] = generateBuildingModelsSummer(n1,dt,K,thetaSummer,t,Psummer,solarSummer,RvalueAttached,hpSizeAttached,floorAreaAttached);
    [eta2_attached,qeWinter_attached] = generateBuildingModelsWinter(zone,n1,K,thetaWinter,Pwinter,solarWinter,floorAreaAttached);


%% generate water heater parameters
[Rw,aw,wwWinter,thetaw,Tsetw,etaw,pMaxHPWH,pMaxR] = generateWaterHeaterModelsWinter(waterfile,L,dt,tStartWinter,tEndWinter);
[wwSummer] = generateWaterHeaterModelsSummer(waterfile,L,dt,tStartSummer,tEndSummer);
[wwFullyear,thetawFullyear] = generateWaterHeaterModels(waterfile,L,dt);




%% electric vehicle parameters


 n2 =  1*round(n1*n2_1vehicle) + 2*round(n1*n2_2vehicle) + 3*(n1-round(n1*n2_0vehicle)-round(n1*n2_1vehicle)-round(n1*n2_2vehicle));

    [a2,etac,etad,eMin,eMax,pcMax,pdMax,tau,e0,p0,ecWinter,atHome,atWork,onRoad,tw,th] = generateEVmodelsWinter(n2,nDays,K,thetaWinter,batteryDeg,dt,t,commuteDistance,percentageSedans);
    [ecSummer] = generateEVmodelsSummer(n2,nDays,K,thetaSummer,batteryDeg,dt,t,commuteDistance,percentageSedans);

    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% baseline simulation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n1a=round(n1*percentageAttached);
n1d=n1-round(n1*percentageAttached);
ntot=n1a+n1d;
eta1(:,1:n1d)=eta1_detached(:,1:n1d);
eta1(:, n1d+1:n1d+n1a)=eta1_attached(:,1:n1a);
eta2(:,1:n1d)=eta2_detached(:,1:n1d);
eta2(:, n1d+1:n1d+n1a)=eta2_attached(:,1:n1a);
pMax1(:,1:n1d)=pMax1_detached(:,1:n1d);
pMax1(:, n1d+1:n1d+n1a)=pMax1_attached(:,1:n1a);
Tset(:,1:n1d)=Tset_detached(:,1:n1d);
Tset(:, n1d+1:n1d+n1a)=Tset_attached(:,1:n1a);
a1(:,1:n1d)=a1_detached(:,1:n1d);
a1(:, n1d+1:n1d+n1a)=a1_attached(:,1:n1a);
R(:,1:n1d)=R_detached(:,1:n1d);
R(:, n1d+1:n1d+n1a)=R_attached(:,1:n1a);
qeWinter(:,1:n1d)=qeWinter_detached(:,1:n1d);
qeWinter(:, n1d+1:n1d+n1a)=qeWinter_attached(:,1:n1a)./4;
qeSummer(:,1:n1d)=qeSummer_detached(:,1:n1d);
qeSummer(:, n1d+1:n1d+n1a)=qeSummer_attached(:,1:n1a)./4;
pMaxAux(:,1:n1d)=pMaxAux_detached(:,1:n1d);
pMaxAux(:, n1d+1:n1d+n1a)=pMaxAux_attached(:,1:n1a);

thetaHigh = f2c(47);
thetaLow = f2c(5);
capHigh = capHigh_Detached_heating;
capLow = capLow_Detached_heating; 
cap = @(x) capLow + (capHigh - capLow)/(thetaHigh - thetaLow)*(x - thetaLow);
cap2_detached_winter = repmat(cap(thetaWinter),1,n1d);

capHigh = capHigh_Attached_heating;
capLow = capLow_Attached_heating; 
cap = @(x) capLow + (capHigh - capLow)/(thetaHigh - thetaLow)*(x - thetaLow);
cap2_attached_winter = repmat(cap(thetaWinter),1,n1a);

cap2Winter(:,1:n1d)=cap2_detached_winter(:,1:n1d);
cap2Winter(:, n1d+1:n1d+n1a)=cap2_attached_winter(:,1:n1a);


capHigh = capHigh_Detached_cooling ;
capLow = capLow_Detached_cooling ;
thetaHigh = f2c(95);
thetaLow = f2c(82);
cap = @(x) capHigh - (capHigh - capLow)/(thetaHigh - thetaLow)*(x - thetaLow);
cap2=cap(thetaSummer);
cap2_detached_summer = repmat(cap(thetaSummer),1,n1d);
cap2_detached_summer_all = repmat(cap(thetaSummer),1,n1);

capHigh =capHigh_Attached_cooling ;
capLow = capLow_Attached_cooling;
cap = @(x) capHigh - (capHigh - capLow)/(thetaHigh - thetaLow)*(x - thetaLow);
cap2_attached_summer = repmat(cap(thetaSummer),1,n1a);
cap2_attached_summer_all = repmat(cap(thetaSummer),1,n1);

cap2Summer(:,1:n1d)=cap2_detached_summer(:,1:n1d);
cap2Summer(:, n1d+1:n1d+n1a)=cap2_attached_summer(:,1:n1a);

%% heat pumps

 [p1baseSummer,HPloadSummer,TbaseSummer] = heatPumpSimulationSummer(K,n1,cap2Summer,eta1,pMaxAux,Tset,a1,R,qeSummer,thetaSummer);
 [p1baseWinter,HPloadWinter,TbaseWinter,p_hp_Cap,auxCap] = heatPumpSimulationWinter(K,n1,cap2Winter,eta2,pMaxAux,Tset,a1,R,qeWinter,thetaWinter);
 

 %% Perform EV simulation for summer and winter
 [phBaseSummer,pwBaseSummer,eBaseSummer,pdBaseSummer] = evSimulationSummer(K,n2,e0,p0,eMin,eMax,pcMax,a2,tau,etac,etad,dt,ecSummer,onRoad,atHome,atWork);
 [phBaseWinter,pwBaseWinter,eBaseWinter,pdBaseWinter] = evSimulationWinter(K,n2,e0,p0,eMin,eMax,pcMax,a2,tau,etac,etad,dt,ecWinter,onRoad,atHome,atWork);

%% heat-pump water heaters
[p1basehpwhSummer,HPWHloadSummer,TbasewSummer,TsetwSummer] = heatPumpWaterHeaterSimulationSummer(K,L,pMaxHPWH,etaw,pMaxR,aw,Rw,wwSummer,thetaw,waterheater);
[p1basehpwhWinter,HPWHloadWinter,TbasewWinter,TsetwWinter] = heatPumpWaterHeaterSimulationWinter(K,L,pMaxHPWH,etaw,pMaxR,aw,Rw,wwWinter,thetaw,waterheater);


%% winter peak
check =0;
[electricPower_detachedWinter,waterHeatingLoadWinter]=winterPeakcalc(zone,thetaWinter,solarWinter,dt,Pwinter,n1,weatherTime,R_detached,cap2_detached_winter,HPWHloadWinter,houseHPdetached,p1basehpwhWinter,floorAreaDetached,pMaxAux_detached,houseAuxdetached,homeWithElectricWH,percentageDetached,check);
check =1;
[electricPower_attachedWinter,~]=winterPeakcalc(zone,thetaWinter,solarWinter,dt,Pwinter,n1,weatherTime,R_attached,cap2_attached_winter,HPWHloadWinter,houseHPattached,p1basehpwhWinter,floorAreaAttached,pMaxAux_attached,houseAuxattached,homeWithElectricWH,percentageAttached,check);

winterPeak =  quantile(  sum(waterHeatingLoadWinter,2) + Pwinter + electricPower_detachedWinter + electricPower_attachedWinter,0.999);


%% summer peak
    check=0;
    [electricPower_detached,~,waterHeatingLoad_detached]=summerPeakcalc(thetaSummer,solarSummer,dt,Psummer,n1,weatherTime,R_detached,cap2_detached_summer_all,HPWHloadSummer,homeWithElectricWH,p1basehpwhSummer,floorAreaDetached,check);
    check=1;
    [electricPower_attached,~,~]=summerPeakcalc(thetaSummer,solarSummer,dt,Psummer,n1,weatherTime,R_attached,cap2_attached_summer_all,HPWHloadSummer,homeWithElectricWH,p1basehpwhSummer,floorAreaAttached,check);
    
 summerPeak = quantile(  sum(waterHeatingLoad_detached,2) + Psummer + sum(electricPower_attached(:,1:round(n1*percentageAttached)),2) + sum(electricPower_detached(:,1:n1-round(n1*percentageAttached)),2),0.999);


if max(summerPeak,winterPeak)<1e3
    s = 1; % unit scalar, divide by s to convert from kW to plot units
    else
    s = 1e3; % for MW
end



%% data compile

data{stateIdx, 14}=max(summerPeak/s,winterPeak/s);
data{stateIdx, 15}=quantile(Pwinter/s + HPWHloadWinter/s + sum(phBaseWinter,2)/s + HPloadWinter/s,0.99);
data{stateIdx, 16}=quantile(Psummer/s + HPWHloadSummer/s + sum(phBaseSummer,2)/s + HPloadSummer/s,0.99);

data{stateIdx, 17}=-data{stateIdx, 16}+data{stateIdx, 15};
data{stateIdx,18} = selectedHeadroom*data{stateIdx, 14} - FutureHeadroom*max(data{stateIdx, 15},data{stateIdx, 16});


upgradeReqMW = ( FutureHeadroom*max(data{stateIdx, 15},data{stateIdx, 16}) - selectedHeadroom*data{stateIdx, 14});

     

if (upgradeReqMW <= 0)
    data{stateIdx,19} =  0;
else
    data{stateIdx,19} = (max(0,960*upgradeReqMW*s/n1));
end
data{stateIdx,20} = zone;

data{stateIdx,21} = housingUnits*data{stateIdx,19};

data{stateIdx,22}= housingUnits;
data{stateIdx,23}= selectedHeadroom;
% Define electrification levels (percent of total homes)
elec_levels = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]; % 20%, 40%, 60%, 80%, 100%



% Initialize array to store results
quantile_values = zeros(length(elec_levels), 1);

% Loop through different electrification levels
for i = 1:length(elec_levels)
    num_selected_ev = round(elec_levels(i) * n2); % Select required number of homes
    num_selected_home = round(elec_levels(i) * n1); % Select required number of homes

    % Randomly select indices
    selected_indices_home = randperm(n1, num_selected_home);
    selected_indices_ev = randperm(n2, num_selected_ev);

    
    % Compute total load for selected homes
    totalLoad_winter = quantile((Pwinter/s) + sum(p1baseWinter(:,selected_indices_home),2)/s + (sum(phBaseWinter(:,selected_indices_ev),2)/s) + (sum(p1basehpwhWinter(:,selected_indices_home),2)/s),0.99);
    totalLoad_summer = quantile((Psummer/s) + sum(p1baseSummer(:,selected_indices_home),2)/s + (sum(phBaseSummer(:,selected_indices_ev),2)/s) + (sum(p1basehpwhSummer(:,selected_indices_home),2)/s),0.99);
    
    % Compute the 99th percentile
   data{stateIdx,23+i} = max(totalLoad_winter,totalLoad_summer);

end





fprintf('State: %s, County: %s City: %s, \n', stateName, countyName, cityName);

end




% Define the data
headers = {
    'State', 'County', 'City', 'lat', 'lng', ...
    'R_detached', 'R_attached','Design Temp Heating (F)','Design Temp Cooling (F)', 'Heat pump size cooling Detached (kW)', ...
    'Heat pump size cooling aatched (kW)','Heat pump size heating detached (kW)','Heat pump size heating attached (kW)','Baseline todays peak  (MW)','Electrified Winter peak (MW)',...
    'Electrified Summer Peak (MW)','Change in electrified peak (MW)','Today - Future (MW)','Cost','zone',...
    'TotalCost','HousingUnits','TodaysHeadroom','10%_elec','20%_elec',...
    '30%_elec','40%_elec','50%_elec','60%_elec','70%_elec','80%_elec','90%_elec','100%'...
   ;
};

% Convert the cell array to a table
dataTable = cell2table(data(1:end, :), 'VariableNames', headers(1, :));

% Create a folder named 'combinedData' if it doesn't exist
folderName = 'Electrification';
if ~exist(folderName, 'dir')
    mkdir(folderName);
end

% Specify the file path including the folder
filePath = fullfile(folderName, sprintf('%s.xlsx', stateName));

%Save the table to an Excel file in the 'combinedData' folder
writetable(dataTable, filePath);



end




toc

