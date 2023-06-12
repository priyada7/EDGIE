%% introduction
% This script simulates the grid impacts of electric vehicles, heat pumps and heat
% pump water heaters on the grid.

clear
tic

%% identifiers that users would like to simulate
batteryDeg =  1; % battery degradation with outside temperature is implemented
matpower   =  1; % set 1 if thermal model of transformer and voltage regulations needs to be solved

rng(1)
%% weather data 
% timing (hour 0 is noon to make EV simulation easier)
ti = 0;                      % initial time, h
nDays=5;                     % no. of days
tf = nDays*24;               % final time, h
dt = 1;                      %  time step, h
K = tf/dt;                   % number of time steps
t = (0:dt:tf)';              % time span, h

%% download weather data
weatherTime = (datetime(2019,1,1,0,0,0):hours(dt):datetime(2020,1,1,0,0,0))';
fileName = 'cambridge-weather-2019.csv';      % load weather file
[thetaFull,solarFull] = importWeather(fileName,weatherTime);
tStart = datetime(2019,1,18,0,0,0);           % Start date for simulation
tEnd = tStart + hours(tf);
tt = (tStart:hours(dt):tEnd)';
theta = thetaFull(weatherTime>=tStart & weatherTime<tEnd);
solar = solarFull(weatherTime>=tStart & weatherTime<tEnd);


%% baseline electrical load data 
% dimensions
n1 =2;               % number of homes (= number of HPs)
carsPerHome = 2;     % number of cars per home
n2 = carsPerHome*n1; % number of EVs
perBushome = 30;     % factor which decides how many house loads will be on each bus of 33 bus network (normally 990 homes gives 30 homes on each bus in a 33 bus network)

% import raw data
fileName = 'MFRED-2019-NYC-Apartments-Electricity-Data.csv';  % load base load file
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
% aggregate baseline power of work neighborhood, kW
workAnnualElectricityIntensity = 17;                        % kWh/sqft/yr, from CBECS 2012 Table PBA4
workSqftPerPerson = 250;
workMeanPower = 2*n1*workSqftPerPerson*workAnnualElectricityIntensity/8760; % kW
workLoadFactor = 0.33;                                      % ratio of peak power to mean power
workPeakPower = workMeanPower/workLoadFactor; % kW
rho = 0.7;                                                  % ratio of daily mean power to daily peak power
pWorkBase = rho*workPeakPower + (1-rho)*workPeakPower*sin(2*pi*(t(1:K)-5)/24);
%% heat pump parameters

%Building data
floorArea = trirnd(200,250,1,n1);                       % conditioned floor area, m^2
ceilingHeight = trirnd(3,3.6,1,n1);                     % ceiling height, m
rhoc = 3.42e-4;                                         % volumetric thermal capacitance of air, kWh/C/m^3
C = trirnd(16,18,1,n1).*(rhoc*floorArea.*ceilingHeight); % thermal capacitance, kWh/C

thetaBase = f2c(50);                                    % base temperature for heading degree-hours, C
HDH = dt*sum(max(0,thetaBase-theta));                   % heating degree-hours in Boston at base temperature 50 F, C*h
RValue = gauss(550,550,1,n1);                                           % m^2*C/kW, value obtained from DC house validation
R = RValue./(4.*ceilingHeight.*sqrt(2.*floorArea));     % thermal resistance, C/kW

a1 = exp(-dt./(R.*C));                                  % discrete-time dynamics parameter
thetaLow = f2c(5);                                      % first temperature point for heat pump COP curve, C
thetaHigh = f2c(45);                                    % second temperature point, C
etaLow = trirnd(1.5,2,1,n1);                            % first COP point
etaHigh = etaLow + 1.5 + trirnd(0,0.5,1,n1);            % second COP point
eta1 = etaLow + (etaHigh-etaLow).*(theta-thetaLow)./(thetaHigh-thetaLow); % COP curve
Tset = repmat(f2c(trirnd(65,74,1,n1)),K+1,1);                             % heating temperature setpoint, C, tuned to cold-climate occupied values in RECS 2015 Table HC6.6
nSetback = 1;                                                             %round(n1/2); % number of units with night setbacks
for i=1:nSetback
    tDown = randi([22,23]);                             % setback time, h
    tUp = randi([7 10]);                                % recovery time, h
    setback = 5*trirnd(3,6,1,1)/9;                      % setback depth, C
    Tset(:,i) = Tset(:,i) - setback;                    % start with setback over all hours
    Tset(mod(t,24)>=tUp & mod(t,24)<=tDown,i) = Tset(mod(t,24)>=tUp & mod(t,24)<=tDown,i) + setback;
end

% qe = P/n1;                                          % exogenous thermal power from non-HP/EV electrical loads, kW
% qe = qe + 0.1*mean(qe)*randn(K,n1)/3;               % random noise
% qe = qe + trirnd(8,12,K,n1).*floorArea.*max(0,sin(2*pi*(t(1:K)-6)/24))/1000; % solar thermal power (8-12 W/m^2)


electricLoadThermalPower = P/n1; % thermal power from electrical loads other than heat pump, kW
bodyThermalPower = gauss(0.1,0.5,K,n1); % thermal power from body heat, kW
solarThermalPower = 0.03*gauss(0.9,1.1,K,n1).*floorArea.*solar; % thermal power from sunlight, kW 
qe = electricLoadThermalPower + bodyThermalPower + solarThermalPower; % exogenous thermal power, kW (8-12 W/m^2)


pMax1 = 4.5*ones(K,n1);                               % max electical capcity of heat pump
pMaxAux = 19.5*ones(K,n1);                            % max electical capcity of resistance backup element


fprintf('Mean exogenous thermal power intensity: %.3g W/m^2.\n',...
    mean(mean(qe)./floorArea)*1000)


%% heat pump water heater ; lumped model
% set timing and number of buildings
L = n1;                             % number of buildings
tfw = tf;                           % 365*24; % final time, h
Kw = K ;                            % number of time steps
startTime = datetime(2019,1,1,0,0,0);           % start time
endTime = datetime(2020,1,1,0,0,0) - hours(dt); % end time
timeSpan = (startTime:hours(dt):endTime)';      % time span as datetime

% generate water heater models
[Rw,aw,ww,thetaw,Tsetw,etaw,pMaxHP,pMaxR] = generateWaterHeaterModels('DHWEventGeneratorOutput.csv',L,timeSpan);

for i=1:length(timeSpan)
    if tStart == timeSpan(i)
        startIndex = i;
    end
    if tEnd == timeSpan(i)
        endIndex = i;
    end
end

%match the time step with weather data
etaw = etaw(startIndex:endIndex-1,:);
ww = ww(startIndex:endIndex-1,:);
thetaw = thetaw(startIndex:endIndex-1,:);


%% electric vehicle parameters
% electric vehicle batteries
rho = trirnd(0.97,0.99,1,n2);       % fraction of charge remaining after 24 h
tau = -24./log(rho);                % dissipation rate, 1/h
a2 = exp(-dt./tau);                 % discrete-time dynamics parameter
etac = trirnd(0.9,0.95,1,n2);       % charge efficiency
etad = etac;                        % discharge efficiency
pcMax = 240*trirnd(24,48,1,n2)/1e3; % charge capacity, kW (level 2: 240 V, 24-32 A)


pdMax = trirnd(30,40,1,n2);                 % discharge capacity, kW (100 mile/h at 0.35 kWh/mile: 35 kW)
eMax = 54 + trirnd(0,28,1,n2);              % energy capacity, kWh (Tesla Model 3, 54-82 kWh)
eMin = trirnd(0.15,0.25,1,n2).*eMax;        % user-specified minimum energy, kWh
e0 = eMin + trirnd(0,1,1,n2).*(eMax-eMin);  % initial energy, kWh
daysPerCharge = 5.5;                        % average number of days between full charges
p0 = (rand(1,n2)<=1/daysPerCharge).*pcMax;  % initial charging power, kW

% electric vehicle commutes
tw = round(trirnd(5,11,1,n2));      % commute to work time of day, h
th = round(trirnd(15,23,1,n2));     % commute to home time of day, h
if batteryDeg ==1
commuteDistance = trirnd(10,20,1,n2);% miles
                                     % source for drivng efficiency modeling https://doi.org/10.1016/j.trd.2019.07.025
drivingEfficiency = zeros(1,nDays);
for i=1:length(theta)
if theta(i,1) < 22
    drivingEfficiency(i) = 0.3392 - 0.005238*theta(i,1) - 0.0001078*theta(i,1)^2 + 1.04710e-5*theta(i,1)^3 + 3.955e-7*theta(i,1)^4-1.362e-8*theta(i,1)^5 - 3.109e-10*theta(i,1)^6;
else
    drivingEfficiency(i)  = 0.4211-0.01627*theta(i,1) + 0.0004229*theta(i,1)^2;
end
end
ec=commuteDistance.*trirnd(min(drivingEfficiency),max(drivingEfficiency),1,n2);
else
    ec = trirnd(3,4,1,n2);          % commute energy, kWh (12 miles, 0.3 miles/kWh)
end
atHome = zeros(K,n2);               % indicator that vehicle's at home
atWork = zeros(K,n2);               % indicator that car's at work
for i=1:n2
    atHome(:,i) = mod(t(1:K),24)<tw(i) | mod(t(1:K),24)>th(i);
    atWork(:,i) = mod(t(1:K),24)>tw(i) & mod(t(1:K),24)<th(i);
end


onRoad = ~atHome & ~atWork; % indicator that vehicle's on the road

%% baseline simulation
% heat pumps

Tbase2 = zeros(K+1,n1);      % indoor temperature, C
Tbase2(1,:) = Tset(1,:);
p1base = zeros(K,n1);        % heat pump electric power, kW
pHpbase2 = zeros(K,n1);      % heat pump electric power, kW
pAux = zeros(K,n1);          % resistance backup electric power, kW

for k=1:K  
%  % to compute power input from Heat Pump
[row1, col1]=size(p1base);
for i = 1 : row1
      Qck(i,:) = (((Tset(i+1,:) - a1.*Tbase2(i,:))./(1-a1) - theta(i))./R - qe(i,:)); 
      for j = 1 : col1
          if Qck(i,j) <= 0
             Qhp(i,j) =0;
             Qck(i,j) =0;
         elseif Qck(i,j) <= eta1(i,j)*pMax1(i,j) & Qck(i,j) >=0
          Qhp(i,j) = Qck(i,j);
          else
              Qhp(i,j) = eta1(i,j)*pMax1(i,j);
          end
          pHpbase2(i,j) = Qhp(i,j)/eta1(i,j);
      end
end
 % to compute power input from Aux
for i = 1 : row1
     Qck2(i,:) = (((Tset(i+1,:) - a1.*Tbase2(i,:))./(1-a1) - theta(i))./R - qe(i,:)); 
     for j = 1 : col1
          if Qck2(i,j) <=0
              Qck2(i,j) =0;
               pAux(i,j) = 0;
          elseif Qck2(i,j) <= eta1(i,j)*pMax1(i,j)
              pAux(i,j) = 0; 
         elseif Qck2(i,j)> eta1(i,j)*pMax1(i,j) & Qck2(i,j)- eta1(i,j)*pMax1(i,j) <= pMaxAux(i,j)
              pAux(i,j) = Qck2(i,j)-eta1(i,j)*pMax1(i,j);
          else
              pAux(i,j) = pMaxAux(i,j);
      end

      end
end
  Tbase2(k+1,:) = a1.*Tbase2(k,:) + (1-a1).*(theta(k) + R.*(eta1(k,:).*pHpbase2(k,:) +  pAux(k,:) + qe(k,:)));
   p1base = pHpbase2+pAux;
end
HPload   = sum(pHpbase2,2) + sum(pAux,2);

%% electric vehicles (charging only at home, discharging only to drive)
eBase = zeros(K+1,n2); % stored energy, kWh
eBase(1,:) = e0;
pcBase = zeros(K,n2); % home or work charge power, kW
pdBase = zeros(K,n2); % discharge power, kW
for k=1:K
    % initialize charging power with previous chargin power
    if k==1
        pck = p0;
    else
        pck = pcBase(k-1,:);
    end
    
    % start charging if energy dropped below minimum
    pck(eBase(k,:)<eMin) = pcMax(eBase(k,:)<eMin);
    
    % stop charging if driving
    pck(~atHome(k,:) & ~atWork(k,:)) = 0;
    
    % saturate charging power
    pcBase(k,:) = max(0,min(pck,(eMax-a2.*eBase(k,:))./(1-a2)./tau./etac));
    
    % discharge power if on the road
    pdBase(k,onRoad(k,:)) = etad(onRoad(k,:)).*ec(onRoad(k,:))/dt;
    
    % update stored energy
    eBase(k+1,:) = a2.*eBase(k,:) + (1-a2).*tau.*(etac.*pcBase(k,:) - pdBase(k,:)./etad);
end
phBase = atHome.*pcBase; % baseline home charge power, kW
pwBase = atWork.*pcBase; % baseline work charge power, kW



 %% water heaters


Qcw2=zeros(Kw,L);
Tsetw = repmat(f2c(trirnd(120,130,1,L)),Kw+1,1);
Tbasew = zeros(Kw+1,L); % water temperature, C
Tbasew(1,:) = Tsetw(1,:);


pAuxhpwh= zeros(Kw+1,L);
p1basew = zeros(Kw,L); % electric power, kW
p1basehpwh = zeros(Kw,L);
for k=1:K
   [row1, col1]=size(p1basew);
 for i = 1 : row1  
      Qcw2(i,:) = ((Tsetw(i+1,:) - aw.*Tbasew(i,:))./(1-aw) - thetaw(i))./Rw + ww(i,:); % power to exactly track setpoint, kW
   for j = 1 : col1
       if Qcw2(i,j) <=0
           Qcw2(i,j)=0;
            pAuxhpwh(i,j) = 0;
       elseif Qcw2(i,j) <= etaw(i,j)*pMaxHP(1,j)
       pAuxhpwh(i,j) = 0;
   elseif Qcw2(i,j) >etaw(i,j)*pMaxHP(1,j) & Qcw2(i,j)-etaw(i,j)*pMaxHP(1,j)<= pMaxR(1,j)
       pAuxhpwh(i,j) = Qcw2(i,j) - etaw(i,j)*pMaxHP(1,j);
   else
       pAuxhpwh(i,j) = pMaxR(1,j);
   end
   end
 end
  

  Tbasew(k+1,:) = aw.*Tbasew(k,:) + (1-aw).*(thetaw(k) + Rw.*( pAuxhpwh(k,:) - ww(k,:) ));
 p1basehpwh = pAuxhpwh;
end

HPWHload = sum(pAuxhpwh,2);

%% baseline plots
% temperature plot
lw = 2; % line width
fs = 24; % font size
tLim = [ti,tf]; % time axis limits, h
tTicks = ti:4:tf; % time axis ticks, h

% temperature plot
figure(1), clf
subplot(2,1,1), plot(t,Tbase2,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks)
ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Heat pumps','fontsize',fs)

% HP power plot
subplot(2,1,2), plot(t(1:K),p1base,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Electric','power (kW)'},'fontsize',fs)
xlabel('Hour','fontsize',fs)

% energy plot
figure(2), clf
subplot(4,1,1), plot(t,eBase./eMax,'linewidth',lw), grid on
hold on, plot(t,mean(eBase./eMax,2),'k','linewidth',2*lw)
xlim(tLim), xticks(tTicks)
ylabel({'State of','charge'},'fontsize',fs), ylim([0,1])
title('Electric vehicles','fontsize',fs)

% home charge plot
subplot(4,1,2), plot(t(1:K),phBase,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Home charge','power (kW)'},'fontsize',fs)

% work charge plot
subplot(4,1,3), plot(t(1:K),pwBase,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Work charge','power (kW)'},'fontsize',fs)

% discharge plot
subplot(4,1,4), plot(t(1:K),pdBase,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Discharge','power (kW)'},'fontsize',fs)
xlabel('Hour','fontsize',fs)

%% baseline aggregate power plot
% set scale of plot
if max(P)<1e3
    s = 1; % unit scalar, divide by s to convert from kW to plot units
    yMax = ceil(max(P/s + sum(phBase,2)/s+sum(p1base,2)/s)/500)*500; % y-axis limit
else
    s = 1e3; % for MW
    yMax = ceil(max(P/s + sum(phBase,2)/s+sum(p1base,2)/s)); % y-axis limit
end


figure(3), clf
fullStairs(tt,P/s,LineWidth=2); hold on
fullStairs(tt,P/s + HPWHload(1:end-1)/s,'m',LineWidth=2),
fullStairs(tt,P/s + HPWHload(1:end-1)/s + sum(phBase,2)/s ,'g',LineWidth=2),
fullStairs(tt,P/s + HPWHload(1:end-1)/s + sum(phBase,2)/s + HPload/s,'k',LineWidth=2)
yline(summerPeak/s,'r',LineWidth=2,LineStyle=":")
if s==1e3, ylabel('Aggregate power (MW)','fontsize',24), end
if s==1, ylabel('Aggregate power (kW)','fontsize',24), end
ylim([0 max(ylim)])
grid on
grid minor
legend('Base Electrified Load','Add HPWH ','Add EV','Add HP','Summer Peak',Location='best')



%% optimization
pie = 0.15; % energy price, $/kWh
pid = 50; % peak demand price, $/kW
dT = 2; % allowable temperature deviation from setpoint, C
pic = 0.5; % discomfort price, $/C
cvx_solver_settings('MIPGapAbs',1e-3,'MIPGap',1e-3,'NumericFocus',3)
cvx_begin 
     variables T(K+1,n1) Tw(K+1,n1) Qdotc(K,n1)  p2(K,n1)  e(K+1,n2) ph(K,n2) pw(K,n2) pd(K,n2)
     expressions p1(K,n1)  % intermediate variables (these are just placeholders for text expressions)
     p1 = Qdotc./eta1 + (1-1./eta1).*pos(Qdotc - eta1.*pMax1); % thermal -> electric power map (note assignment (=) vs. equality (==) operator)
   
     minimize( pid*max(P + sum(p1,2)+ sum(p2,2) + sum(ph - atHome.*pd,2))...        % home neighborhood peak
        + pid*max(0,max(pWorkBase + sum(pw - atWork.*pd,2)) - max(pWorkBase))...    % work neighborhood peak
        + pie*sum(P + sum(p1,2) + sum(p2,2) + sum(ph + pw,2))...                    % energy
        + pic*(sum(sum(abs(T-Tset))) -  sum(sum(abs(Tbase2-Tset)))) ...
        +  pic*(sum(sum(abs(Tw-Tsetw))) -  sum(sum(abs(Tbasew-Tsetw)))))            % discomfort
    % discomfort
 subject to                                                                         % constraints
%heat pump
 T(1,:) == Tset(1,:);                                                               % initial condition
 T(K+1,:) == Tset(K+1,:);                                                           % final condition: return to setpoint
 T(2:K+1,:) == repmat(a1,K,1).*T(1:K,:)...
     + (1-repmat(a1,K,1)).*(repmat(theta,1,n1) + repmat(R,K,1).*(Qdotc + qe));
 zeros(K,n1) <= Qdotc <= (eta1.*pMax1 + pMaxAux).*ones(K,n1)                        % heat pump electric power capacity limits
 abs(T - Tset) <= dT

 %heat-pump water heater
 Tw(1,:) == Tsetw(1,:);                                                             % initial condition
 Tw(K+1,:) == Tsetw(K+1,:);                                                         % final condition: return to setpoint
 Tw(2:K+1,:) == repmat(aw,K,1).*Tw(1:K,:)...
     + (1-repmat(aw,K,1)).*(repmat(thetaw,1,n1) + repmat(Rw,K,1).*(p2 - ww));
 abs(Tw - Tsetw) <= 10
 zeros(K,n1) <= p2 <= pMaxR.*ones(K,n1)


 % EV energy
    e(1,:) == e0                                                                    % initial condition
    e(K+1,:) == 0.3.*eMax(1,:)                                                    %eBase(K+1,:) % terminal condition: return to initial charge
    e(2:K+1,:) == repmat(a2,K,1).*e(1:K,:)...
        + (1-repmat(a2,K,1)).*repmat(tau,K,1).*...
        (repmat(etac,K,1).*(ph + pw) - pd./repmat(etad,K,1))
    repmat(eMin,K+1,1) <= e <= repmat(eMax,K+1,1)
    
    % EV power
    zeros(K,n2) <= ph <= repmat(pcMax,K,1)
    zeros(K,n2) <= pw <= repmat(pcMax,K,1)
    zeros(K,n2) <= pd <= repmat(pdMax,K,1)
    for i=1:n2
        pd(mod(t,24)==tw(i),i) == etad(i)*ec(i)/dt                                  % commute to work
        pd(mod(t,24)==th(i),i) == etad(i)*ec(i)/dt                                  % commute to home
    end
    ph(~atHome) == 0 % at-home charging
    pw(~atWork) == 0 % at-work charging
        
 
cvx_end
fprintf('CVX exit status: %s.\n',cvx_status)
% expand sparse matrices
e = full(e);
ph = full(ph);
pw = full(pw);
pd = full(pd);
T = full(T);
p1 = full(p1);
p2 = full(p2);
toc
%%
if cvx_status == 'Failed'
    return
end
%% optimized plots
% temperature plot
figure(1), clf
subplot(2,1,1), plot(t,T,'linewidth',lw), grid on
xlim(tLim), xticks(tTicks)
ylabel({'Indoor','temperature ($^\circ$C)'},'fontsize',fs)
title('Heat pumps','fontsize',fs)

% HP power plot
subplot(2,1,2), plot(t(1:K),p1,'linewidth',lw), grid on
hold on, plot(t(1:K),mean(p1,2),'k','linewidth',2*lw)
xlim(tLim), xticks(tTicks), ylim([0,max(ylim)])
ylabel({'Electric','power (kW)'},'fontsize',fs)
xlabel('Hour of day (0 = midnight)','fontsize',fs)

% individual energy
figure(2), clf
subplot(4,1,1), stairs(t,e./eMax,'linewidth',lw), grid on
hold on, stairs(t,mean(e./eMax,2),'k','linewidth',2*lw)
ylim([0,max(ylim)]), ylabel({'State of','charge'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)

% individual home charge power
subplot(4,1,2), plot(t(1:K),ph,'linewidth',lw), grid on
ylim([0,max(ylim)]), ylabel({'Home charge','power (kW)'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)

% individual work charge power
subplot(4,1,3), plot(t(1:K),pw,'linewidth',lw), grid on
ylim([0,max(ylim)]), ylabel({'Work charge','power (kW)'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)

% individual discharge power
subplot(4,1,4), plot(t(1:K),pd,'linewidth',lw), grid on
ylim([0,max(ylim)]), ylabel({'Discharge','power (kW)'},'fontsize',fs)
xlim([0,tf]), xticks(tTicks)
xlim([0,tf]), xticks(tTicks), xlabel('Hour','fontsize',fs)

%% optimized aggregate power plots
figure(3), clf 
subplot(2,1,1),fullStairs(t,P + sum(p1base,2) + sum(phBase,2),'m','linewidth',lw), grid on
hold on, fullStairs(t,P + sum(p1,2) + sum(ph - atHome.*pd,2),'k','linewidth',lw)
stairs(t,0*t+summerPeak,'b--','linewidth',lw)
ylim([0,max(ylim)])
ylabel({'Aggregate home','power (kW)'},'fontsize',fs)

subplot(2,1,2), fullStairs(t,pWorkBase,'m','linewidth',lw), grid on
hold on, fullStairs(t,pWorkBase+sum(pw - atWork.*pd,2),'k','linewidth',lw)
ylim([0,max(ylim)])
ylabel({'Aggregate work','power (kW)'},'fontsize',fs)

%% baseline aggregate home power plots, peak day only
% day selection
tLim = [0 24];
tShift = t - 2*24;
fs = 30;
resolution = 300; % plot printing resolution
toPrint = 0;

% plot
figure(4), clf
fullStairs(tShift,P/s,'b','linewidth',lw), grid on % base power
hold on, fullStairs(tShift,P/s + sum(phBase,2)/s,'k','linewidth',lw) % EV power
fullStairs(tShift,P/s + sum(p1base,2)/s,'r','linewidth',lw) % HP power
fullStairs(tShift,P/s + sum(phBase,2)/s+sum(p1base,2)/s,'m','linewidth',lw) % EV + HP power
fullStairs(tShift,0*P + summerPeak/s,'b--','linewidth',lw) % summer peak
xlim(tLim), xticks(tTicks), ylim([0,yMax])
if s==1e3, ylabel('Aggregate power (MW)','fontsize',fs), end
if s==1, ylabel('Aggregate power (kW)','fontsize',fs), end
xlabel('Hour of day (0 = midnight)','fontsize',fs)


%% optimized aggregate home power plots, peak day only
% plot
figure(5), clf
fullStairs(tShift,P + sum(p1base,2) + sum(phBase,2),'m','linewidth',lw), grid on
hold on, fullStairs(tShift,P + sum(p1,2) + sum(ph - atHome.*pd,2),'k','linewidth',lw)
stairs(tShift,0*tShift+summerPeak,'b--','linewidth',lw)
ylim([0,max(ylim)])
ylabel('Aggregate power (kW)','fontsize',fs)
xlim(tLim), xticks(tTicks), ylim([0,yMax])
if s==1e3, ylabel('Aggregate power (MW)','fontsize',fs), end
if s==1, ylabel('Aggregate power (kW)','fontsize',fs), end
xlabel('Hour of day (0 = midnight)','fontsize',fs)



figure(6), clf 
fullStairs(tt,P/s + sum(p1base,2)/s + sum(phBase,2)/s + sum(p1basehpwh(1:end-1,:),2)/s,'k','linewidth',lw), grid on
hold on, fullStairs(tt,P/s + sum(p1,2)/s + sum(ph - atHome.*pd,2)/s + sum(p2,2)/s ,'m--','linewidth',lw)
ylim([0,max(ylim)])
grid on
grid minor
if s==1e3, ylabel({'Aggregate',' power (MW)'},'fontsize',fs), end
if s==1, ylabel('Aggregate power (kW)','fontsize',fs), end
legend('Unoptimized','Optimized',Location="best")

%% Matpower load calculations using case33bw and transformer modeling
if matpower==1
    Unopt=phBase;
    Opt=ph-atHome.*pd;


    ind=1;
    BUnopt=zeros(nDays*24,n2);
    BOpt=zeros(nDays*24,n2);
    for i=1:(n2-1)
        BUnopt(:,i) =Unopt(:,ind)+Unopt(:,ind+1);
        BOpt(:,i)=Opt(:,ind)+Opt(:,ind+1);
        ind=ind+1;
    end

    j=(1:2:length(Unopt));
    CUnopt=zeros(nDays*24,n1);
    COpt=zeros(nDays*24,n1);
    for i=1:n1
        CUnopt(:,i) = BUnopt(:,j(i));
        COpt(:,i) = BOpt(:,j(i));
    end

    baseload =P'/n1;
    for i=1:n1
        baseload(i,:) = baseload(1,:);
    end
    loadprofileUnopt = (CUnopt + p1base + p1basehpwh(1:end-1,:))'+ baseload;
    loadprofileOpt = (COpt + p1 + p2 )'+ baseload;

    multiplier =1;
    matpowerloadUnopt= (multiplier*loadprofileUnopt);

    matpowerloadOpt= (multiplier*loadprofileOpt);



    beta = (1:perBushome:n1+perBushome);
    for i=1:length(beta)-1
        houseloadUnopt(i,:)= sum(matpowerloadUnopt(beta(:,i):beta(:,i+1)-1,:));
        houseloadOpt(i,:)= sum(matpowerloadOpt(beta(:,i):beta(:,i+1)-1,:));
    end

    mpcUnopt=loadcase('case33bw');
    mpcUnopt.gen(:,6)=1.05 ;  %substation voltage
    for i=1:120
        mpcUnopt.bus(:,3)=houseloadUnopt(:,i)/1e3;
        mpcUnopt.bus(:,4)=0.35.*houseloadUnopt(:,i)/1e3;
        resultsUnopt(:,i)=runpf(mpcUnopt);
        generationDataUnopt(:,i) = resultsUnopt(i).gen(:,2);
        loadDataUnopt = sum(houseloadUnopt)/1e3;
        LossDataUnopt(i,:) = generationDataUnopt(:,i) - loadDataUnopt(:,i);
        LossDataPercUnopt(i,:) = (generationDataUnopt(:,i)- loadDataUnopt(:,i))*100/generationDataUnopt(:,i);
    end

    mpcOpt=loadcase('case33bw');
    mpcOpt.gen(:,6)=1.05   %substation voltage
    for i=1:120
        mpcOpt.bus(:,3)=houseloadOpt(:,i)/1e3;
        mpcOpt.bus(:,4)=0.35.*houseloadOpt(:,i)/1e3;
        resultsOpt(:,i)=runpf(mpcOpt);
        generationDataOpt(:,i) = resultsOpt(i).gen(:,2);
        loadDataOpt = sum(houseloadOpt)/1e3;
        LossDataOpt(i,:) = generationDataOpt(:,i) - loadDataOpt(:,i);
        LossDataPercOpt(i,:) = (generationDataOpt(:,i)- loadDataOpt(:,i))*100/generationDataOpt(:,i);
    end



    figure(8)
    ylim([0.5 1.05]);
    [maxLoss, index1] = max(LossDataUnopt(:,1));
    plot(resultsUnopt(index1).bus(:,1),resultsUnopt(index1).bus(:,8));
    xlabel('Bus','fontsize',fs)
    ylabel('Voltage in p.u.','fontsize',fs)
   
%     figure(9)
%     yyaxis left
%     plot(t(2:end),LossDataUnopt,'LineWidth',4);
%     xlabel('Hour','fontsize',fs)
%     ylabel('Power Loss (MW)','fontsize',fs)
%     yyaxis right
%     plot(t(2:end),LossDataPercUnopt,'LineWidth',4);
%     ylabel('Power Loss %','fontsize',fs)

    figure(10)
    ylim([0.5 1.05]);
    [maxLoss, index2] = max(LossDataOpt(:,1))
    plot(resultsOpt(index2).bus(:,1),resultsOpt(index2).bus(:,8));
    xlabel('Bus','fontsize',fs)
    ylabel('Voltage in p.u.','fontsize',fs)
   
%     figure(11)
%     yyaxis left
%     plot(t(2:end),LossDataOpt,'LineWidth',4);
%     xlabel('Hour','fontsize',fs)
%     ylabel('Power Loss (MW)','fontsize',fs)
%     yyaxis right
%     plot(t(2:end),LossDataPercOpt,'LineWidth',4);
%     ylabel('Power Loss %','fontsize',fs)

    figure(12)
    plot(resultsUnopt(index1).bus(:,1),resultsUnopt(index1).bus(:,8),'LineWidth',2); hold on
    plot(resultsOpt(index2).bus(:,1),resultsOpt(index2).bus(:,8),'LineWidth',2);
    xlabel(" Bus ",'fontsize',fs);
    ylabel('Voltage in p.u.','fontsize',fs)
    legend("Unoptimized","Optimized")
    ylim([0 max(ylim)])
    hold off

%     figure(13)
%     plot(t(2:end),LossDataUnopt,'LineWidth',4);hold on
%     plot(t(2:end),LossDataOpt,'LineWidth',4);
%     xlabel('Hour','fontsize',fs)
%     ylabel('Power Loss (MW)','fontsize',fs)
%     legend("Unoptimized","Optimized")
%     toc
%     hold off

    Irated = 36600;


    IactualUnopt = sum(houseloadUnopt)/0.240;
    IactualOpt = sum(houseloadOpt)/0.240;
    Tamb=theta;
    IpuUnopt = (round(IactualUnopt/Irated,1))';
    IpuOpt = (round(IactualOpt/Irated,1))';




    [TOT1, HST1, time] = transformer_git(Tamb, IpuUnopt);
    figure(14)
    plot(time,TOT1,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Top oil temperature (C)','fontsize',fs)

    figure(15)
    plot(t(1:end-1),HST1,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Hot spot temperature (C)','fontsize',fs)


    [TOT2, HST2, time] = transformer_git(Tamb, IpuOpt);
    figure(16)
    plot(time,TOT2,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Top oil temperature (C)','fontsize',fs)

    figure(17)
    plot(t(1:end-1),HST2,'LineWidth',3);
    xlabel('Hour','fontsize',fs)
    ylabel('Hot spot temperature (C)','fontsize',fs)

for i=1:length(HST2)
lossUnoptT(i,1) = exp(15000/383-15000/(HST1(i,1)+273));
end


fprintf('Loss of life UNOPT: %.5g h.\n',sum(lossUnoptT))

for i=1:length(HST2)
lossOptT(i,1) = exp(15000/383-15000/(HST2(i,1)+273));
end
fprintf('Loss of life OPT: %.5g h.\n',sum(lossOptT))


end
